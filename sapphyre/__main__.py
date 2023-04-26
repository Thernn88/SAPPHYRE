#!/usr/bin/env python
"""
Order:
    1. Prepare
    2. Diamond
    3. Reporter
Post-processing:
    5. align
    6. pal2nal
    7. FlexCull (optional)
    8. Outlier
    9. Merge
    10. Combine
"""
import argparse


class CaseInsensitiveArgumentParser(argparse.ArgumentParser):
    def _parse_known_args(self, arg_strings, *args, **kwargs):
        # Iterate through the subparsers to find the command
        lower_arg_string = list(map(str.lower, arg_strings))
        for action in self._actions:
            if isinstance(action, argparse._SubParsersAction):
                subparsers = action.choices
                for subparser in subparsers:
                    # Check if the command matches the argument string (case-insensitive)
                    if subparser.lower() in lower_arg_string:
                        # Replace the argument string with the command (case-sensitive)
                        arg_strings[
                            lower_arg_string.index(subparser.lower())
                        ] = subparser
                        break
                break
        return super()._parse_known_args(arg_strings, *args, **kwargs)


def subcmd_prepare(subparsers):
    par = subparsers.add_parser(
        "Prepare",
        help="Loads NT input files (.fa, .fas or .fasta) into a rocksdb database. "
        "Unique NT sequences stored with a duplicate count stored for later use. "
        "Performs six-fold translation of base NT sequences into AA translation. "
        "Unique AA translation stored with duplicate counts stored for later use.",
    )
    par.add_argument("INPUT", help="Path to directory of Input folder", action="store")
    par.add_argument(
        "-clear",
        "--clear_database",
        action="store_true",
        help="Overwrite existing rocksdb database.",
    )
    par.add_argument(
        "-ml",
        "--minimum_sequence_length",
        default=90,
        type=int,
        help="Minimum input sequence length.",
    )
    par.add_argument(
        "-cs",
        "--chunk_size",
        default=500000,
        type=int,
        help="Sequences per chunk in DB.",
    )
    par.add_argument(
        "-k",
        "--keep_prepared",
        action="store_true",
        help="Writes the prepared input fasta into the output taxa directory.",
    )
    par.set_defaults(func=prepare, formathelp=par.format_help)


def prepare(args):
    from . import prepare

    if not prepare.main(args):
        print()
        print(args.formathelp())


def subcmd_diamond(subparsers):
    par = subparsers.add_parser(
        "Diamond",
        help="To be written.",
    )
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument(
        "-ovw", "--overwrite", action="store_true", help="Overwrite existing files."
    )
    par.add_argument("-d", "--debug", action="store_true", help="Enable debug out.")

    # Should be based on Taxa name. Check TODO.
    par.add_argument(
        "-strict",
        "--strict-search-mode",
        action="store_true",
        help="Only allow sequences that hit on every present taxa.",
    )
    par.add_argument(
        "-s",
        "--sensitivity",
        choices=["very", "ultra"],
        default="very",
        help="Diamond sensitivty.",
    )
    par.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=0.000005,
        help="Diamond blast evalue threshold.",
    )
    par.add_argument(
        "-t",
        "--top",
        type=int,
        default=10,
        help="Diamond top percentage cull.",
    )
    par.add_argument(
        "-tr",
        "--top-ref",
        type=float,
        default=0.1,
        help="Dynamically adjusts %% of hits a reference is less than our top 5 and still a good ref.",
    )
    par.add_argument(
        "-ip",
        "--internal-percent",
        type=float,
        default=0.3,
        help="Percentage of overlap required to constitute an internal overlap kick.",
    )
    par.set_defaults(func=diamond, formathelp=par.format_help)


def diamond(args):
    from . import diamond

    if not diamond.main(args):
        print()
        print(args.formathelp())


def subcmd_reporter(subparsers):
    par = subparsers.add_parser(
        "reporter",
        help="Trims mapped sequence to mapped region." "Produces aa and nt output.",
    )
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument(
        "-bp",
        "--minimum_bp",
        type=int,
        default=20,
        help="Amount of bp required after trim.",
    )
    par.add_argument(
        "-m",
        "--matches",
        type=int,
        default=7,
        help="Amount of matches for dynamic pairwise aligned edge trim.",
    )
    par.add_argument(
        "--gene_list_file",
        type=str,
        default=None,
        help="Path to a txt file containing target genes for processing. Processes all genes in the input folder if not specified.",
    )
    par.add_argument(
        "-bm",
        "--blosum_mode",
        choices=["exact", "strict", "lax"],
        default="strict",
        help="Trim distance mode.",
    )
    par.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="Enable debug",
    )
    par.add_argument(
        "-clear",
        "--clear_output",
        action="store_false",
        default=True,
        help="Clear output folder before running.",
    )
    par.set_defaults(func=reporter, formathelp=par.format_help)


def reporter(args):
    from . import reporter

    mainargs = reporter.MainArgs(
        args.verbose,
        args.processes,
        args.debug,
        args.INPUT,
        args.orthoset_input,
        args.orthoset,
        args.compress,
        args.matches,
        args.blosum_mode,
        args.minimum_bp,
        args.gene_list_file,
        args.clear_output,
    )
    if not reporter.main(mainargs):
        print(args.formathelp())


def subcmd_outlier(subparsers):
    par = subparsers.add_parser(
        "outlier",
        help="Calculates a Blosum62 distance matrix which are used to remove outlier "
        "sequences above a threshold.",
    )
    par.add_argument("INPUT", help="Path to taxa", action="extend", nargs="+")
    par.add_argument("-o", "--output", default="outlier", help="Output folder")
    par.add_argument(
        "-t",
        "--threshold",
        type=int,
        default=100,
        help="Greater than reference mean to be counted as an outlier. Default is 2x.",
    )
    par.add_argument(
        "--no-references",
        action="store_true",
        help="Disable output of reference sequences",
    )
    par.add_argument(
        "-ccp",
        "--col-cull-percent",
        type=float,
        default=0.33,
        help="Minimum percentage of data required not to cull column.",
    )
    par.add_argument(
        "-refcp",
        "--ref-gap-percent",
        type=float,
        default=0.5,
        help="Minimum percent of non-gap characters in a reference sequence after cull to pass.",
    )
    par.add_argument(
        "-refmp",
        "--ref-min-percent",
        type=float,
        default=0.5,
        help="Minimum percent of references required after kick",
    )
    par.add_argument(
        "-imp",
        "--index_group_min_bp",
        type=int,
        default=20,
        help="Minimum bp for index group after column cull.",
    )
    par.add_argument(
        "-ict",
        "--internal_consensus_threshold",
        type=float,
        default=0.85,
        help="Consensus threshold for internal summary.",
    )
    par.add_argument(
        "-ikt",
        "--internal_kick_threshold",
        type=int,
        default=4,
        help="Amount of mismatches required to constitute an internal kick.",
    )
    par.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="Log outliers to csv files",
    )
    par.set_defaults(func=outlier, formathelp=par.format_help)


def outlier(args):
    from . import outlier

    if not outlier.main(args):
        print()
        print(args.formathelp)


def subcmd_Merge(subparsers):
    par = subparsers.add_parser(
        "Merge",
        help="Reference-guided De-novo Assembly Algorithm which merges overlapping reads "
        "into contiguous segments (Contigs).",
    )
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument(
        "-aa",
        "--aa_input",
        type=str,
        default="outlier/aa",
        help="Path to directory of AA folder",
    )
    par.add_argument(
        "-nt",
        "--nt_input",
        type=str,
        default="outlier/nt",
        help="Path to directory of NT folder",
    )
    par.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable debug. When enabled displays each component of merged headers.",
    )
    par.add_argument(
        "-io",
        "--ignore_overlap_chunks",
        action="store_true",
        default=False,
        help="Ignore overlapping chunks and merge all candidates for a reference taxon.",
    )
    par.add_argument(
        "-m",
        "--majority",
        type=float,
        default=0.55,
        help="Percentage for majority ruling.",
    )
    par.add_argument(
        "-mc",
        "--majority_count",
        type=int,
        default=4,
        help="Percentage for majority ruling.",
    )
    par.set_defaults(func=Merge, formathelp=par.format_help)


def Merge(args):
    from . import merge

    if not merge.main(args):
        print()
        print(args.formathelp())


def subcmd_Combine(subparsers):
    dsc = (
        "Initial Dataset Construction: Merges AA and NT data files across "
        "multiple taxa into a single AA and NT pair per gene."
    )
    par = subparsers.add_parser("Combine", help=dsc, description=dsc)
    par.add_argument("INPUT", help="Paths of directories.", action="extend", nargs="+")
    par.add_argument(
        "-t",
        "--prepend-directory",
        action="store",
        type=str,
        dest="DIRECTORY",
        help="Prepend DIRECTORY to the list of INPUT.",
    )
    par.add_argument(
        "-out",
        "--output-directory",
        type=str,
        required=True,
        help="Target directory for merged output. (Example: MergedGenes)",
    )
    par.set_defaults(func=Combine, formathelp=par.format_help)


def Combine(args):
    from . import combine

    if not combine.main(args):
        print()
        print(args.formathelp())


def subcmd_align(subparsers):
    par = subparsers.add_parser(
        "align", help="Aligns AA sequences against existing reference alignment."
    )
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="Output intermediate files for debug",
    )
    par.set_defaults(func=align, formathelp=par.format_help)


def align(args):
    from . import align

    if not align.main(args):
        print()
        print(args.formathelp())


def subcmd_pal2nal(subparsers):
    par = subparsers.add_parser(
        "Pal2Nal",
        help="Mirrors Amino Acid Alignment to Nucleotide data using specified NCBI "
        "table. Also performs basic error checking on data.",
    )
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument("-t", "--table", type=int, default=1, help="Table ID.")
    par.set_defaults(func=pal2nal, formathelp=par.format_help)


def pal2nal(args):
    from . import pal2nal

    if not pal2nal.main(args):
        print()
        print(args.formathelp())


def subcmd_flexcull(subparsers):
    par = subparsers.add_parser(
        "FlexCull",
        help="Adaptive End Trimming algorithm that trims the start and end of "
        "candidate reads to remove introns and/or other LQ bases.",
    )
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument(
        "-o", "--output", type=str, default="trimmed", help="Output Directory."
    )
    par.add_argument(
        "-aa", "--amino-acid", type=str, default="align", help="AA Folder Name."
    )
    par.add_argument(
        "-nt", "--nucleotide", type=str, default="nt_aligned", help="NT Folder Name."
    )  #
    par.add_argument(
        "-bm",
        "--blosum_mode",
        choices=["exact", "strict", "lax"],
        default="exact",
        help="Blosum strictness setting.",
    )
    par.add_argument(
        "-mm",
        "--mismatches",
        type=int,
        default=0,
        help="Amount mismatches allowed per trim.",
    )
    par.add_argument(
        "-m",
        "--matches",
        type=int,
        default=5,
        help="Amount of base pairs that have to match reference.",
    )  #
    par.add_argument(
        "-cc",
        "--column_cull",
        type=float,
        default=0.1,
        help="Percentage of reference columns that must contain data.",
    )
    par.add_argument(
        "-gt",
        "--gap_threshold",
        type=float,
        default=0.50,
        help="Percentage of references that must contain a gap to allow a match to continue.",
    )
    par.add_argument(
        "-bp", "--base-pair", type=int, default=20, help="Minimum bp after cull."
    )
    par.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable debug. When enabled Output log of culls.",
    )
    par.set_defaults(func=flexcull, formathelp=par.format_help)


def flexcull(args):
    from . import flexcull

    flexargs = flexcull.MainArgs(
        args.verbose,
        args.processes,
        args.debug,
        args.INPUT,
        args.output,
        args.amino_acid,
        args.nucleotide,
        args.matches,
        args.base_pair,
        args.compress,
        args.gap_threshold,
        args.mismatches,
        args.column_cull,
        args.blosum_mode,
    )
    if not flexcull.main(flexargs):
        print()
        print(args.formathelp())


def finalize(args):
    from . import finalize

    if not finalize.main(args):
        print(args.formathelp())


def subcmd_finalize(subparsers):
    par = subparsers.add_parser(
        "Finalize",
        help="Contains a handful of useful functions to finalize a dataset. "
        "Kick taxa: Exclude taxa listed in the --kick_file txt file. "
        "Kick columns: Removes columns where there is greater than or equal to "
        "a supplied percentage of gap characters over data characters. Also "
        "kicks sequences with less than a supplied integer of data characters. "
        "supplied percentage: -kp, supplied integer: -bp"
        "Stop codon: Replaces stop characters with X for AA and N for NT. "
        "Rename taxon: Uses the supplied --names_csv file to replace reference taxon "
        "names using the taxa id. Renames all ids in col 0 of the csv to col 1"
        "Sort: Sorts genes based on presence in the supplied target txt file. "
        "Concat: Merges all the final sequences together into a single fasta file.",
    )

    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument(
        "-k", "--kick_file", type=str, default="TaxaKick.txt", help="Percent"
    )
    par.add_argument(
        "-t",
        "--target_file",
        type=str,
        default="TARGET_GENES_SYRPHID.txt",
        help="Percent",
    )
    par.add_argument("-n", "--names_csv", type=str, default="names.csv", help="Percent")
    par.add_argument(
        "-kp",
        "--kick_percentage",
        type=float,
        default=0.85,
        help="Adjustable value for kick columns. Float value of minimum percentage of non-gap characters",
    )
    par.add_argument(
        "-mb",
        "--minimum_bp",
        type=int,
        default=30,
        help="Adjustable value for kick columns. Integer value of minimum base pairs",
    )
    par.add_argument(
        "-s",
        "--sort",
        action="store_true",
        help="Sort taxa based on target file provided.",
    )
    par.add_argument(
        "-kt",
        "--kick_taxa",
        action="store_true",
        help="Kick taxa present in kick file provided.",
    )
    par.add_argument(
        "-kc",
        "--kick_columns",
        action="store_true",
        help="Kick columns based on amount of gap characters.",
    )
    par.add_argument(
        "-stop",
        "--stopcodon",
        action="store_true",
        help="Replace Stop Codons with X for AA and N for NT.",
    )
    par.add_argument(
        "-r",
        "--rename",
        action="store_true",
        help="Rename reference taxa using names csv file provided.",
    )
    par.add_argument(
        "-con",
        "--concat",
        action="store_true",
        help="Merge resulting on target genes into a final fasta file.",
    )
    par.set_defaults(func=finalize, formathelp=par.format_help)


def subcmd_download(sp):
    par = sp.add_parser(
        "download", help="Download fastq files from www.ncbi.nlm.nih.gov"
    )
    par.add_argument(
        "INPUT",
        help="Path to the CSV file input",
    )
    par.add_argument(
        "-b",
        "--bin",
        help="Path to SRA Toolkit. Will try system's PATH if not used.",
        required=False,
    )
    par.set_defaults(func=download, formathelp=par.format_help)


def download(argsobj):
    from . import download

    if not download.main(argsobj):
        print()
        print(argsobj.formathelp())


def subcmd_archiver(sp):
    par = sp.add_parser("Archiver", help="Recursive archiver/unarchiver")  # TODO add me
    par.add_argument(
        "INPUT",
        help="Path to input directory.",
        action="extend",
        nargs="+",
    )
    par.add_argument(
        "-sd",
        "--specific_directories",
        help="Directories to archive/unarchive.",
        action="extend",
        nargs="+",
    )
    par.add_argument(
        "-u",
        "--unarchive",
        default=False,
        action="store_true",
        help="Unarchive directories.",
    )
    par.set_defaults(func=archiver, formathelp=par.format_help)


def archiver(argsobj):
    from . import archiver

    if not archiver.main(argsobj):
        print()
        print(argsobj.formathelp())


def subcmd_makeref(sp):
    par = sp.add_parser("Makeref", help="Reference set generation tool")  # TODO add me
    par.add_argument(
        "INPUT",
        type=str,
        help="Input file for current set (Either .fa or .sqlite)",
    )
    par.add_argument(
        "-m",
        "--align_method",
        choices=["clustal", "mafft"],
        default="clustal",
        help="What alignment method to use.",
    )
    par.add_argument(
        "-s", "--set", type=str, help="Name of the current set being produced/modified."
    )
    par.add_argument(
        "-od",
        "--orthoset_dir",
        type=str,
        default="orthosets",
        help="Path to the Orthosets dir.",
    )
    par.add_argument(
        "-ovw", "--overwrite", action="store_true", help="Overwrite existing files."
    )
    par.set_defaults(func=makeref, formathelp=par.format_help)


def makeref(argsobj):
    from . import makeref

    if not makeref.main(argsobj):
        print()
        print(argsobj.formathelp())


def subcmd_wrap_final(sp):
    par = sp.add_parser("reconcile", help="Wrapper for Second Run")  # TODO add me
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument("-t", "--table", type=int, default=1, help="Table ID.")
    par.add_argument(
        "-aa",
        "--aa_input",
        type=str,
        default="align",
        help="Path to directory of AA folder",
    )
    par.add_argument(
        "-nt",
        "--nt_input",
        type=str,
        default="nt_aligned",
        help="Path to directory of NT folder",
    )
    par.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable debug. When enabled displays each component of merged headers.",
    )
    par.add_argument(
        "-io",
        "--ignore_overlap_chunks",
        action="store_false",
        default=True,
        help="Ignore overlapping chunks and merge all candidates for a reference taxon.",
    )
    par.add_argument(
        "-m",
        "--majority",
        type=float,
        default=0.55,
        help="Percentage for majority ruling.",
    )
    par.add_argument(
        "-mc",
        "--majority_count",
        type=int,
        default=4,
        help="Sites required for majority ruling.",
    )
    par.set_defaults(func=wrap_final, formathelp=par.format_help)


def wrap_final(argsobj):
    from . import merge, pal2nal, align

    print("Triggering Align")
    if not align.main(argsobj):
        print()
        print(argsobj.formathelp())
    print("Triggering Pal2Nal")
    if not pal2nal.main(argsobj):
        print()
        print(argsobj.formathelp())
    print("Triggering Merge")
    if not merge.main(argsobj):
        print()
        print(argsobj.formathelp())


if __name__ == "__main__":
    parser = CaseInsensitiveArgumentParser(
        prog="sapphyre",
        # TODO write me
        description="Order: Prepare, Hmmsearch, BlastPal, Reporter, "
        "align, pal2nal, FlexCull (optional), outlier, Merge, Combine",
        epilog="sapphyre  Copyright (C) 2022  Sapphyre Team\n"
        "License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.\n"
        "This program comes with ABSOLUTELY NO WARRANTY.\n"
        "This is free software, and you are welcome to redistribute it "
        "under certain conditions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--compress",
        action="store_true",
        default=False,
        help="Output fasta files as compressed files using gzip",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbosity level. Repeat for increased verbosity.",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="Number of threads used to call processes.",
    )
    parser.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="orthosets",
        help="Path to directory of Orthosets folder.",
    )
    parser.add_argument(
        "-os",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Current Orthoset to be used.",
    )
    subparsers = parser.add_subparsers()
    # The order in which those functions are called define the order in which
    # the subcommands will be displayed.
    subcmd_prepare(subparsers)
    subcmd_diamond(subparsers)
    subcmd_reporter(subparsers)

    # align > pal2nal > FlexCull (optional) > outlier > Merge
    subcmd_align(subparsers)
    subcmd_pal2nal(subparsers)
    subcmd_flexcull(subparsers)
    subcmd_outlier(subparsers)
    subcmd_Merge(subparsers)
    subcmd_Combine(subparsers)
    subcmd_download(subparsers)
    subcmd_archiver(subparsers)

    # Finalize
    subcmd_finalize(subparsers)

    # Makeref
    subcmd_makeref(subparsers)

    # Final wrapper
    subcmd_wrap_final(subparsers)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        parser.exit()
    args.func(args)
