#!/usr/bin/env python
"""
Order:
    1. PrepareDB
    2. HmmsearchDB
    3. BlastPalDB
    4. ReporterDB
Post-processing:
    5. mafft
    6. pal2nal
    7. FlexCull (optional)
    8. OutlierCheck
    9. MergeOverlap
    10. MergeGenes
"""
import argparse


def subcmd_preparedb(subparsers):
    parser_preparedb = subparsers.add_parser(
        "PrepareDB",
        help="Loads NT input files (.fa, .fas or .fasta) into a rocksdb database. "
        "Unique NT sequences stored with a duplicate count stored for later use. "
        "Performs six-fold translation of base NT sequences into AA translation. "
        "Unique AA translation stored with duplicate counts stored for later use.",
    )
    parser_preparedb.add_argument(
        "INPUT", help="Path to directory of Input folder", action="store"
    )
    parser_preparedb.add_argument(
        "-c",
        "--clear_database",
        action="store_true",
        help="Overwrite existing rocksdb database.",
    )
    parser_preparedb.add_argument(
        "-ml",
        "--minimum_sequence_length",
        default=90,
        type=int,
        help="Minimum input sequence length.",
    )
    parser_preparedb.add_argument(
        "-sl",
        "--sequences_per_level",
        default=100000,
        type=int,
        help="Amount of sequences to store per database entry.",
    )
    parser_preparedb.add_argument(
        "-k",
        "--keep_prepared",
        action="store_true",
        help="Writes the prepared input fasta into the output taxa directory.",
    )
    parser_preparedb.set_defaults(
        func=preparedb, formathelp=parser_preparedb.format_help
    )


def preparedb(args):
    from . import preparedb

    if not preparedb.main(args):
        print()
        print(args.formathelp())


def subcmd_hmmsearchdb(subparsers):
    parser_hmmsearchdb = subparsers.add_parser(
        "HmmsearchDB",
        help="Queries protein translations against profile HMMs using the HMMER3 "
        "external. Filters HMMER3 output using 3 custom filters: MultiGene, "
        "InternalMulti & InternalSubpar to remove LQ hits and prevent sequence reuse. "
        "Loads hits into RocksDB after filters finish.",
    )
    parser_hmmsearchdb.add_argument(
        "INPUT", help="Path to input directory.", action="extend", nargs="+"
    )
    parser_hmmsearchdb.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="PhyMMR/orthosets",
        help="Path to directory of Orthosets folder",
    )
    parser_hmmsearchdb.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    parser_hmmsearchdb.add_argument(
        "-ovw",
        "--overwrite",
        default=False,
        action="store_true",
        help="Remake domtbl files even if previous file exists.",
    )
    parser_hmmsearchdb.add_argument(
        "-s", "--score", type=float, default=40, help="Score threshold. Defaults to 40"
    )
    parser_hmmsearchdb.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=0,
        help="Evalue threshold. Defaults to 0",
    )
    parser_hmmsearchdb.add_argument(
        "--excluded-list",
        default=False,
        help="File containing names of genes to be excluded",
    )
    parser_hmmsearchdb.add_argument(
        "--wanted-list", default=False, help="File containing list of wanted genes"
    )
    parser_hmmsearchdb.add_argument(
        "--remake-protfile",
        default=False,
        action="store_true",
        help="Force creation of a new protfile even if one already exists.",
    )
    parser_hmmsearchdb.add_argument(
        "-sdm",
        "--score_diff_multi",
        type=float,
        default=1.05,
        help="Multi-gene Score Difference Adjustment",
    )
    parser_hmmsearchdb.add_argument(
        "-mom",
        "--min_overlap_multi",
        type=float,
        default=0.3,
        help="Multi-gene Minimum Overlap Adjustment",
    )
    parser_hmmsearchdb.add_argument(
        "-momi",
        "--minimum_overlap_internal_multi",
        type=float,
        default=0.5,
        help="Internal Multi Minimum Overlap Adjustment",
    )
    parser_hmmsearchdb.add_argument(
        "-sdi",
        "--score_diff_internal",
        type=float,
        default=1.5,
        help="Internal Score Difference Adjustmen",
    )
    parser_hmmsearchdb.add_argument(
        "-moi",
        "--min_overlap_internal",
        type=float,
        default=0.9,
        help="Internal Minimum Overlap Adjustment",
    )
    parser_hmmsearchdb.add_argument(
        "-m",
        "--max_hmm_batch_size",
        default=250000,
        type=int,
        help="Max hits per hmmsearch batch in db. Default: 500 thousand.",
    )
    parser_hmmsearchdb.add_argument(
        "-d", "--debug", type=int, default=0, help="Output debug logs."
    )
    parser_hmmsearchdb.set_defaults(
        func=hmmsearchdb, formathelp=parser_hmmsearchdb.format_help
    )


def hmmsearchdb(args):
    from . import hmmsearchdb
    if not hmmsearchdb.main(args):
        print(args.formathelp())


def subcmd_blastpaldb(subparsers):
    parser_blastpal = subparsers.add_parser(
        "BlastPalDB",
        help="Blasts Hmmsearch hits using NCBI-Blast against reference sequences. "
        "Loads resulting output into RocksDB",
    )

    parser_blastpal.add_argument(
        "INPUT",
        action="extend",
        nargs="+",
        help="Path to directory of Input folder",
    )
    parser_blastpal.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="PhyMMR/orthosets",
        help="Path to directory of Orthosets folder",
    )
    parser_blastpal.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    parser_blastpal.add_argument(
        "-bs",
        "--blast_minimum_score",
        type=float,
        default=40.0,
        help="Minimum score filter in blast.",
    )
    parser_blastpal.add_argument(
        "-be",
        "--blast_minimum_evalue",
        type=float,
        default=0.00001,
        help="Minimum evalue filter in blast.",
    )
    parser_blastpal.add_argument(
        "-ovw",
        "--overwrite",
        action="store_true",
        help="Overwrite existing blast results.",
    )
    parser_blastpal.set_defaults(
        func=blastpaldb, formathelp=parser_blastpal.format_help
    )


def blastpaldb(args):
    from . import blastpaldb
    if not blastpaldb.main(args):
        print(args.formathelp())


def subcmd_reporterdb(subparsers):
    parser_reporterdb = subparsers.add_parser(
        "ReporterDB",
        help="Checks Blast results to ensure a hit is reciprocal. Queries a sequence "
        "using exonerate to align it against a target reference and trim it to mapped "
        "region. Produces aa and nt output.",
    )
    parser_reporterdb.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    parser_reporterdb.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="PhyMMR/orthosets",
        help="Path to directory of Orthosets folder",
    )
    parser_reporterdb.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    parser_reporterdb.add_argument(
        "-ml", "--min_length", type=int, default=30, help="Minimum Transcript Length"
    )
    parser_reporterdb.add_argument(
        "-ms", "--min_score", type=float, default=40, help="Minimum Hit Domain Score"
    )
    parser.add_argument("-d", "--debug", type=int, default=0, help="Verbose debug.")
    parser_reporterdb.set_defaults(
        func=reporterdb, formathelp=parser_reporterdb.format_help
    )


def reporterdb(args):
    from . import reporterdb
    mainargs = reporterdb.MainArgs(
        args.verbose,
        args.processes,
        args.debug,
        args.INPUT,
        args.orthoset_input,
        args.orthoset,
        args.min_length,
        args.min_score
    )
    if not reporterdb.main(mainargs):
        print(args.formathelp())


def subcmd_outliercheck(subparsers):
    parser_outliercheck = subparsers.add_parser(
        "OutlierCheck",
        help="Calculates a Blosum62 distance matrix which are used to remove outlier "
        "sequences above a threshold.",
    )
    parser_outliercheck.add_argument(
        "INPUT", help="Path to taxa", action="extend", nargs="+"
    )
    parser_outliercheck.add_argument(
        "-o", "--output", default="outlier", help="Output folder"
    )
    parser_outliercheck.add_argument(
        "-t",
        "--threshold",
        type=int,
        default=50,
        help="Greater than reference mean to be counted as an outlier. Default is 2x.",
    )
    parser_outliercheck.add_argument(
        "--no-references",
        action="store_true",
        help="Disable output of reference sequences",
    )
    parser_outliercheck.add_argument(
        "-s",
        "--sort",
        choices=["cluster", "original"],
        default="original",
        help="Sort candidate output by cluster and taxa, or preserver original order.",
    )
    parser_outliercheck.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="Log outliers to csv files",
    )
    parser_outliercheck.set_defaults(
        func=outliercheck, formathelp=parser_outliercheck.format_help
    )


def outliercheck(args):
    from . import outliercheck
    if not outliercheck.main(args):
        print()
        print(args.formathelp)


def subcmd_mergeoverlap(subparsers):
    parser_mergeoverlap = subparsers.add_parser(
        "MergeOverlap",
        help="Reference-guided De-novo Assembly Algorithm which merges overlapping reads "
        "into contiguous segments (Contigs).",
    )
    parser_mergeoverlap.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    parser_mergeoverlap.add_argument(
        "-aa",
        "--aa_input",
        type=str,
        default="aa",
        help="Path to directory of AA folder",
    )
    parser_mergeoverlap.add_argument(
        "-nt",
        "--nt_input",
        type=str,
        default="nt",
        help="Path to directory of NT folder",
    )
    parser_mergeoverlap.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable debug. When enabled displays each component of merged headers.",
    )
    parser_mergeoverlap.add_argument(
        "-c",
        "--comparison",
        type=str,
        default="Drosophila_melanogaster",
        help="Fallback Comparison Taxa. Sequence in which Sequence A and Sequence B is "
        "compared to in split calculation.",
    )
    parser_mergeoverlap.add_argument(
        "-m",
        "--majority",
        type=float,
        default=0.66,
        help="Percentage for majority ruling.",
    )
    parser_mergeoverlap.add_argument(
        "-mc",
        "--majority_count",
        type=int,
        default=4,
        help="Percentage for majority ruling.",
    )
    parser_mergeoverlap.set_defaults(
        func=mergeoverlap, formathelp=parser_mergeoverlap.format_help
    )


def mergeoverlap(args):
    from . import merge_overlap
    if not merge_overlap.main(args):
        print()
        print(args.formathelp())


def subcmd_mergegenes(subparsers):
    dsc = "Initial Dataset Construction: Merges AA and NT data files across "\
          "multiple taxa into a single AA and NT pair per gene."
    parser_mergegenes = subparsers.add_parser(
        "MergeGenes", help=dsc, description=dsc
    )
    parser_mergegenes.add_argument(
        "INPUT", help="Paths of directories.", action="extend", nargs="+"
    )
    parser_mergegenes.add_argument(
        "-t", "--prepend-directory", action="store", type=str, dest="DIRECTORY",
        help="Prepend DIRECTORY to the list of INPUT."
    )
    parser_mergegenes.add_argument(
        "-o", "--output-directory", type=str, required=True,
        help="Target directory for merged output. (Example: MergedGenes)"
    )
    parser_mergegenes.set_defaults(
        func=mergegenes, formathelp=parser_mergegenes.format_help
    )


def mergegenes(args):
    from . import merge_genes
    if not merge_genes.main(args):
        print()
        print(args.formathelp())


def subcmd_mafft(subparsers):
    parser_mafft = subparsers.add_parser(
        "mafft", help="Aligns AA sequences against existing reference alignment."
    )
    parser_mafft.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    parser_mafft.add_argument(
        "-oi",
        "--orthoset_input",
        type=str,
        default="PhyMMR/orthosets",
        help="Path to directory of Orthosets folder",
    )
    parser_mafft.add_argument(
        "-o",
        "--orthoset",
        type=str,
        default="Ortholog_set_Mecopterida_v4",
        help="Orthoset",
    )
    parser_mafft.set_defaults(func=mafft, formathelp=parser_mafft.format_help)


def mafft(args):
    from . import mafft
    if not mafft.main(args):
        print()
        print(args.formathelp())


def subcmd_pal2nal(subparsers):
    parser_pal2nal = subparsers.add_parser(
        "Pal2Nal",
        help="Mirrors Amino Acid Alignment to Nucleotide data using specified NCBI "
        "table. Also performs basic error checking on data.",
    )
    parser_pal2nal.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    parser_pal2nal.add_argument("-t", "--table", type=int, default=1, help="Table ID.")
    parser_pal2nal.set_defaults(func=pal2nal, formathelp=parser_pal2nal.format_help)


def pal2nal(args):
    from . import pal2nal
    if not pal2nal.main(args):
        print()
        print(args.formathelp())


def subcmd_flexcull(subparsers):
    parser_flexcull = subparsers.add_parser(
        "FlexCull",
        help="Adaptive End Trimming algorithm that trims the start and end of "
        "candidate reads to remove introns and/or other LQ bases.",
    )
    parser_flexcull.add_argument(
        "INPUT", help="Path to directory of Input folder", action="extend", nargs="+"
    )
    parser_flexcull.add_argument(
        "-o", "--output", type=str, default="trimmed", help="Output Directory."
    )
    parser_flexcull.add_argument(
        "-aa", "--amino-acid", type=str, default="mafft", help="AA Folder Name."
    )
    parser_flexcull.add_argument(
        "-nt", "--nucleotide", type=str, default="nt_aligned", help="NT Folder Name."
    )
    parser_flexcull.add_argument(
        "-m",
        "--matches",
        type=int,
        default=3,
        help="Amount of nucleotides that have to match reference.",
    )
    parser_flexcull.add_argument(
        "-bp", "--base-pair", type=int, default=15, help="Minimum bp after cull."
    )
    parser_flexcull.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable debug. When enabled Output log of culls.",
    )
    parser_flexcull.set_defaults(func=flexcull, formathelp=parser_flexcull.format_help)


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
        description="Order: PrepareDB, HmmsearchDB, BlastPalDB, ReporterDB, "
        "mafft, pal2nal, FlexCull (optional), OutlierCheck, MergeOverlap, MergeGenes",
        epilog="I am epilog",  # TODO write me
        # formatter_class=argparse.RawDescriptionHelpFormatter,
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
    subcmd_preparedb(subparsers)
    subcmd_hmmsearchdb(subparsers)
    subcmd_blastpaldb(subparsers)
    subcmd_reporterdb(subparsers)

    # mafft > pal2nal > FlexCull (optional) > OutlierCheck > MergeOverlap
    subcmd_mafft(subparsers)
    subcmd_pal2nal(subparsers)
    subcmd_flexcull(subparsers)
    subcmd_outliercheck(subparsers)
    subcmd_mergeoverlap(subparsers)
    subcmd_mergegenes(subparsers)

    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:
        parser.print_help()
