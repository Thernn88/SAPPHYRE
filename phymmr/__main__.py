#!/usr/bin/env python
"""
Order:
    1. Prepare
    2. Hmmsearch
    3. BlastPal
    4. Reporter
Post-processing:
    5. mafft
    6. pal2nal
    7. FlexCull (optional)
    8. OutlierCheck
    9. MergeOverlap
    10. MergeGenes
"""
import argparse


def subcmd_prepare(subparsers):
    par = subparsers.add_parser(
        "Prepare",
        help="Loads NT input files (.fa, .fas or .fasta) into a rocksdb database. "
        "Unique NT sequences stored with a duplicate count stored for later use. "
        "Performs six-fold translation of base NT sequences into AA translation. "
        "Unique AA translation stored with duplicate counts stored for later use.",
    )
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="store"
    )
    par.add_argument(
        "-c",
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
        "-sl",
        "--sequences_per_level",
        default=100000,
        type=int,
        help="Amount of sequences to store per database entry.",
    )
    par.add_argument(
        "-k",
        "--keep_prepared",
        action="store_true",
        help="Writes the prepared input fasta into the output taxa directory.",
    )
    par.set_defaults(
        func=prepare, formathelp=par.format_help
    )


def prepare(args):
    from . import prepare

    if not prepare.main(args):
        print()
        print(args.formathelp())


def subcmd_hmmsearch(subparsers):
    par = subparsers.add_parser(
        "Hmmsearch",
        help="Queries protein translations against profile HMMs using the HMMER3 "
        "external. Filters HMMER3 output using 3 custom filters: MultiGene, "
        "InternalMulti & InternalSubpar to remove LQ hits and prevent sequence reuse. "
        "Loads hits into RocksDB after filters finish.",
    )
    par.add_argument(
        "INPUT", help="Path to input directory.", action="extend", nargs="+"
    )
    par.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="orthosets",
        help="Path to directory of Orthosets folder",
    )
    par.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    par.add_argument(
        "-ovw",
        "--overwrite",
        default=False,
        action="store_true",
        help="Remake domtbl files even if previous file exists.",
    )
    par.add_argument(
        "-s", "--score", type=float, default=40, help="Score threshold. Defaults to 40"
    )
    par.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=0,
        help="Evalue threshold. Defaults to 0",
    )
    par.add_argument(
        "--excluded-list",
        default=False,
        help="File containing names of genes to be excluded",
    )
    par.add_argument(
        "--wanted-list", default=False, help="File containing list of wanted genes"
    )
    par.add_argument(
        "--remake-protfile",
        default=False,
        action="store_true",
        help="Force creation of a new protfile even if one already exists.",
    )
    par.add_argument(
        "-sdm",
        "--score_diff_multi",
        type=float,
        default=1.05,
        help="Multi-gene Score Difference Adjustment",
    )
    par.add_argument(
        "-mom",
        "--min_overlap_multi",
        type=float,
        default=0.3,
        help="Multi-gene Minimum Overlap Adjustment",
    )
    par.add_argument(
        "-momi",
        "--minimum_overlap_internal_multi",
        type=float,
        default=0.5,
        help="Internal Multi Minimum Overlap Adjustment",
    )
    par.add_argument(
        "-sdi",
        "--score_diff_internal",
        type=float,
        default=1.5,
        help="Internal Score Difference Adjustmen",
    )
    par.add_argument(
        "-moi",
        "--min_overlap_internal",
        type=float,
        default=0.9,
        help="Internal Minimum Overlap Adjustment",
    )
    par.add_argument(
        "-m",
        "--max_hmm_batch_size",
        default=250000,
        type=int,
        help="Max hits per hmmsearch batch in db. Default: 500 thousand.",
    )
    par.add_argument(
        "-d", "--debug", type=int, default=0, help="Output debug logs."
    )
    par.set_defaults(
        func=hmmsearch, formathelp=par.format_help
    )


def hmmsearch(args):
    from . import hmmsearch
    if not hmmsearch.main(args):
        print(args.formathelp())


def subcmd_blastpal(subparsers):
    par = subparsers.add_parser(
        "BlastPal",
        help="Blasts Hmmsearch hits using NCBI-Blast against reference sequences. "
        "Loads resulting output into RocksDB",
    )

    par.add_argument(
        "INPUT",
        action="extend",
        nargs="+",
        help="Path to directory of Input folder",
    )
    par.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="orthosets",
        help="Path to directory of Orthosets folder",
    )
    par.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    par.add_argument(
        "-bs",
        "--blast_minimum_score",
        type=float,
        default=40.0,
        help="Minimum score filter in blast.",
    )
    par.add_argument(
        "-be",
        "--blast_minimum_evalue",
        type=float,
        default=0.00001,
        help="Minimum evalue filter in blast.",
    )
    par.add_argument(
        "-ovw",
        "--overwrite",
        action="store_true",
        help="Overwrite existing blast results.",
    )
    par.set_defaults(
        func=blastpal, formathelp=par.format_help
    )


def blastpal(args):
    from . import blastpal
    if not blastpal.main(args):
        print(args.formathelp())


def subcmd_reporter(subparsers):
    par = subparsers.add_parser(
        "Reporter",
        help="Checks Blast results to ensure a hit is reciprocal. Queries a sequence "
        "using exonerate to align it against a target reference and trim it to mapped "
        "region. Produces aa and nt output.",
    )
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="orthosets",
        help="Path to directory of Orthosets folder",
    )
    par.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    par.add_argument(
        "-ml", "--min_length", type=int, default=30, help="Minimum Transcript Length"
    )
    par.add_argument(
        "-ms", "--min_score", type=float, default=40, help="Minimum Hit Domain Score"
    )
    par.add_argument("-d", "--debug", type=int, default=0, help="Verbose debug.")
    par.set_defaults(
        func=reporter, formathelp=par.format_help
    )


def reporter(args):
    from . import reporter
    mainargs = reporter.MainArgs(
        args.verbose,
        args.processes,
        args.debug,
        args.INPUT,
        args.orthoset_input,
        args.orthoset,
        args.min_length,
        args.min_score
    )
    if not reporter.main(mainargs):
        print(args.formathelp())


def subcmd_outliercheck(subparsers):
    par = subparsers.add_parser(
        "OutlierCheck",
        help="Calculates a Blosum62 distance matrix which are used to remove outlier "
        "sequences above a threshold.",
    )
    par.add_argument(
        "INPUT", help="Path to taxa", action="extend", nargs="+"
    )
    par.add_argument(
        "-o", "--output", default="outlier", help="Output folder"
    )
    par.add_argument(
        "-t",
        "--threshold",
        type=int,
        default=50,
        help="Greater than reference mean to be counted as an outlier. Default is 2x.",
    )
    par.add_argument(
        "--no-references",
        action="store_true",
        help="Disable output of reference sequences",
    )
    par.add_argument(
        "-s",
        "--sort",
        choices=["cluster", "original"],
        default="original",
        help="Sort candidate output by cluster and taxa, or preserver original order.",
    )
    par.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="Log outliers to csv files",
    )
    par.set_defaults(
        func=outliercheck, formathelp=par.format_help
    )


def outliercheck(args):
    from . import outliercheck
    if not outliercheck.main(args):
        print()
        print(args.formathelp)


def subcmd_mergeoverlap(subparsers):
    par = subparsers.add_parser(
        "MergeOverlap",
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
        default="aa",
        help="Path to directory of AA folder",
    )
    par.add_argument(
        "-nt",
        "--nt_input",
        type=str,
        default="nt",
        help="Path to directory of NT folder",
    )
    par.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable debug. When enabled displays each component of merged headers.",
    )
    par.add_argument(
        "-c",
        "--comparison",
        type=str,
        default="Drosophila_melanogaster",
        help="Fallback Comparison Taxa. Sequence in which Sequence A and Sequence B is "
        "compared to in split calculation.",
    )
    par.add_argument(
        "-fo",
        "--force_overlap_chunks",
        action="store_true", 
        default=False,
        help="Only merge overlapping chunks within each reference taxon candidate merge."
    )
    par.add_argument(
        "-m",
        "--majority",
        type=float,
        default=0.66,
        help="Percentage for majority ruling.",
    )
    par.add_argument(
        "-mc",
        "--majority_count",
        type=int,
        default=4,
        help="Percentage for majority ruling.",
    )
    par.set_defaults(
        func=mergeoverlap, formathelp=par.format_help
    )


def mergeoverlap(args):
    from . import merge_overlap
    if not merge_overlap.main(args):
        print()
        print(args.formathelp())


def subcmd_mergegenes(subparsers):
    dsc = "Initial Dataset Construction: Merges AA and NT data files across "\
          "multiple taxa into a single AA and NT pair per gene."
    par = subparsers.add_parser(
        "MergeGenes", help=dsc, description=dsc
    )
    par.add_argument(
        "INPUT", help="Paths of directories.", action="extend", nargs="+"
    )
    par.add_argument(
        "-t", "--prepend-directory", action="store", type=str, dest="DIRECTORY",
        help="Prepend DIRECTORY to the list of INPUT."
    )
    par.add_argument(
        "-o", "--output-directory", type=str, required=True,
        help="Target directory for merged output. (Example: MergedGenes)"
    )
    par.set_defaults(
        func=mergegenes, formathelp=par.format_help
    )


def mergegenes(args):
    from . import merge_genes
    if not merge_genes.main(args):
        print()
        print(args.formathelp())


def subcmd_mafft(subparsers):
    par = subparsers.add_parser(
        "mafft", help="Aligns AA sequences against existing reference alignment."
    )
    par.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    par.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="orthosets",
        help="Path to directory of Orthosets folder",
    )
    par.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    par.add_argument(
        "-l", "--linsi", action="store_true", default=False,
        help="Enable the use of mafft-linsi."
    )
    par.set_defaults(func=mafft, formathelp=par.format_help)


def mafft(args):
    from . import mafft
    if not mafft.main(args):
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
        "-aa", "--amino-acid", type=str, default="mafft", help="AA Folder Name."
    )
    par.add_argument(
        "-nt", "--nucleotide", type=str, default="nt_aligned", help="NT Folder Name."
    )
    par.add_argument(
        "-m",
        "--matches",
        type=int,
        default=3,
        help="Amount of nucleotides that have to match reference.",
    )
    par.add_argument(
        "-bp", "--base-pair", type=int, default=15, help="Minimum bp after cull."
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
        args.base_pair
    )
    if not flexcull.main(flexargs):
        print()
        print(args.formathelp())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="phymmr",
        # TODO write me
        description="Order: Prepare, Hmmsearch, BlastPal, Reporter, "
        "mafft, pal2nal, FlexCull (optional), OutlierCheck, MergeOverlap, MergeGenes",
        epilog="phymmr  Copyright (C) 2022  PhyMMR Team\n"
        "License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>.\n"
        "This program comes with ABSOLUTELY NO WARRANTY.\n"
        "This is free software, and you are welcome to redistribute it "
        "under certain conditions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
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

    subparsers = parser.add_subparsers()
    # The order in which those functions are called define the order in which
    # the subcommands will be displayed.
    subcmd_prepare(subparsers)
    subcmd_hmmsearch(subparsers)
    subcmd_blastpal(subparsers)
    subcmd_reporter(subparsers)

    # mafft > pal2nal > FlexCull (optional) > OutlierCheck > MergeOverlap
    subcmd_mafft(subparsers)
    subcmd_pal2nal(subparsers)
    subcmd_flexcull(subparsers)
    subcmd_outliercheck(subparsers)
    subcmd_mergeoverlap(subparsers)
    subcmd_mergegenes(subparsers)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        parser.exit()
    args.func(args)
