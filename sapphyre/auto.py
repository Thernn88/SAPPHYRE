import argparse
import json
import os
from glob import glob
from .timekeeper import KeeperMode, TimeKeeper

def main(args):
    default_config = {
        "prepare": {
            "clear_database": True,
            "minimum_sequence_length": 90,
            "chunk_size": 500000,
            "keep_prepared": False,
        },
        "diamond": {
            "overwrite": False,
            "debug": False,
            "strict_search_mode": True,
            "sensitivity": "very",
            "evalue": 6,
            "top": 10,
            "top_ref": 0.1,
            "internal_percent": 0.3,
        },
        "reporter": {
            "minimum_bp": 20,
            "matches": 7,
            "gene_list_file": None,
            "blosum_mode": "strict",
            "debug": False,
            "clear_output": True,
        },
        "align": {
            "debug": False,
            "second_run": False,
            "add_fragments": True,
        },
        "pal2nal": {"table": 1},
        "flexcull": {
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
        "outlier": {
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
        "merge": {
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
        with open(args.config) as f:
            config = json.load(f)

    start_script = args.start.title()

    global_args = vars(args)
    global_args.pop("config")
    global_args.pop("start")

    scripts = list(config.keys())
    if start_script.lower() not in scripts:
        options = ", ".join(scripts)
        print(
            f"{start_script} is not a valid start script. Please choose from: ({options})",
        )
        return False
    scripts = scripts[scripts.index(start_script.lower()) :]

    time = TimeKeeper(KeeperMode.DIRECT)
    for script in scripts:
        sargs = config[script]
        print(f"Executing: {script.title()}")
        this_args = global_args.copy()
        this_args.update(sargs)
        this_args = argparse.Namespace(**this_args)
        this_args.INPUT = (
            args.INPUT
            if script == "prepare"
            else sorted(
                glob(
                        os.path.join("datasets", os.path.split(args.INPUT)[-1], "*.fa"),
                    ),
                )
        )

        if script == "prepare":
            from . import prepare

            if not prepare.main(this_args):
                print("Error in Prepare.")
        elif script == "diamond":
            from . import diamond

            if not diamond.main(this_args):
                print("Error in Diamond.")
        elif script == "reporter":
            from . import reporter

            if not reporter.main(this_args):
                print("Error in Reporter.")
        elif script == "Align":
            from . import align

            if not align.main(this_args):
                print("Error in Align.")
        elif script == "pal2nal":
            from . import pal2nal

            if not pal2nal.main(this_args):
                print("Error in Pal2Nal.")
        elif script == "flexcull":
            from . import flexcull

            if not flexcull.main(this_args):
                print("Error in FlexCull.")
        elif script == "outlier":
            from . import outlier

            if not outlier.main(this_args):
                print("Error in Outlier.")
        elif script == "merge":
            from . import merge

            if not merge.main(this_args):
                print("Error in Merge.")

    time = time.differential()
    print(f"Took {time:.2f}s overall.")

    return True
