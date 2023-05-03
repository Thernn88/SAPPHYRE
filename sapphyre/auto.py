import os
import json
import argparse
from glob import glob


def main(args):
    default_config = {
        "Prepare": {
            "clear_database": True,
            "minimum_sequence_length": 90,
            "chunk_size": 500000,
            "keep_prepared": False,
        },
        "Diamond": {
            "overwrite": False,
            "debug": False,
            "strict_search_mode": True,
            "sensitivity": "very",
            "evalue": 0.000005,
            "top": 10,
            "top_ref": 0.1,
            "internal_percent": 0.3,
        },
        "Reporter": {
            "minimum_bp": 20,
            "matches": 7,
            "gene_list_file": None,
            "blosum_mode": "strict",
            "debug": False,
            "clear_output": True,
        },
        "Align": {"debug": False},
        "Pal2Nal": {"table": 1},
        "FlexCull": {
            "output": "trimmed",
            "amino_acid": "align",
            "nucleotide": "nt_aligned",
            "blosum_strictness": "exact",
            "mismatches": 1,
            "matches": 5,
            "column_cull": 0.1,
            "gap_threshold": 0.5,
            "base_pair": 20,
            "debug": False,
        },
        "Outlier": {
            "output": "outlier",
            "threshold": 100,
            "no_references": False,
            "col_cull_percent": 0.33,
            "ref_gap_percent": 0.5,
            "ref_min_percent": 0.5,
            "index_group_min_bp": 20,
            "internal_consensus_threshold": 0.85,
            "internal_kick_threshold": 4,
            "debug": False,
        },
        "Merge": {
            "aa_input": "outlier/aa",
            "nt_input": "outlier/nt",
            "debug": False,
            "ignore_overlap_chunks": False,
            "majority": 0.55,
            "majority_count": 4,
        },
    }

    config = default_config
    if args.config and not os.path.exists(args.config):
        print("Config file does not exist.")
        choice = input("Generate template or run default? (g/r): ")
        if choice.lower() == "g":
            print("Generating config file...")
            with open(args.config, "w") as f:
                json.dump(default_config, f, indent=4)
            return True
    elif args.config:
        with open(args.config, "r") as f:
            config = json.load(f)

    global_args = vars(args)
    global_args.pop("config")

    for script, sargs in config.items():
        print(f"Executing: {script}")
        this_args = global_args.copy()
        this_args.update(sargs)
        this_args = argparse.Namespace(**this_args)

        if script == "Prepare":
            from . import prepare

            this_args.INPUT = args.INPUT

            if not prepare.main(this_args):
                print("Error in Prepare.")
        elif script == "Diamond":
            from . import diamond

            this_args.INPUT = list(glob(os.path.join("datasets", os.path.split(args.INPUT)[-1], "*.fa")))

            if not diamond.main(this_args):
                print("Error in Diamond.")
        elif script == "Reporter":
            from . import reporter

            this_args.INPUT = list(glob(os.path.join("datasets", os.path.split(args.INPUT)[-1], "*.fa")))

            if not reporter.main(this_args):
                print("Error in Reporter.")
        elif script == "Align":
            from . import align

            this_args.INPUT = list(glob(os.path.join("datasets", os.path.split(args.INPUT)[-1], "*.fa")))

            if not align.main(this_args):
                print("Error in Align.")
        elif script == "Pal2Nal":
            from . import pal2nal

            this_args.INPUT = list(glob(os.path.join("datasets", os.path.split(args.INPUT)[-1], "*.fa")))

            if not pal2nal.main(this_args):
                print("Error in Pal2Nal.")
        elif script == "FlexCull":
            from . import flexcull

            this_args.INPUT = list(glob(os.path.join("datasets", os.path.split(args.INPUT)[-1], "*.fa")))

            if not flexcull.main(this_args):
                print("Error in FlexCull.")
        elif script == "Outlier":
            from . import outlier

            this_args.INPUT = list(glob(os.path.join("datasets", os.path.split(args.INPUT)[-1], "*.fa")))

            if not outlier.main(this_args):
                print("Error in Outlier.")
        elif script == "Merge":
            from . import merge

            this_args.INPUT = list(glob(os.path.join("datasets", os.path.split(args.INPUT)[-1], "*.fa")))

            if not merge.main(this_args):
                print("Error in Merge.")

    return True
