import argparse
import json
import os
from glob import glob
import gc
from .__main__ import (
    align_args,
    diamond_args,
    flexcull_args,
    # motif_args,
    merge_args,
    outlier_args,
    pal2nal_args,
    prepare_args,
    reporter_args,
    hmmsearch_args
)
from .timekeeper import KeeperMode, TimeKeeper


def get_args(arg_function):
    parser = argparse.ArgumentParser()
    arg_function(parser)
    variables = {}
    for action in parser._actions:
        if action.dest != "help":
            variables[action.dest] = action.default

    return variables


def main(args):
    gfm = args.gene_finding_mode == 2
    default_config = {
        "prepare": get_args(prepare_args),
        "diamond": get_args(diamond_args),
        "hmmsearch": get_args(hmmsearch_args),
        "reporter": get_args(reporter_args),
        "align": get_args(align_args),
        "pal2nal": get_args(pal2nal_args),
        "flexcull": get_args(flexcull_args),
        # "motif": get_args(motif_args),
        "outlier": get_args(outlier_args),
        "merge": get_args(merge_args),
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
        if gfm and (script == "merge" or script == "flexcull"):
            continue
        sargs = config[script]
        if script != "motif":
            print(f"\nExecuting: {script.title()}")
        this_args = global_args.copy()
        this_args.update(sargs)
        
        if args.assembly:
            if 'assembly' in this_args:
                this_args['assembly'] = True
        
        if args.overwrite:
            if 'overwrite' in this_args:
                this_args['overwrite'] = True

        this_args = argparse.Namespace(**this_args)
        for input_path in this_args.INPUT:
            if not os.path.isdir(input_path) or not os.path.exists(input_path):
                continue
            
            if args.solo:
                this_args.INPUT = [input_path]
            elif script == "prepare":
                this_args.INPUT = input_path
            else:
                this_args.INPUT = sorted(
                    glob(
                        os.path.join(args.in_folder, os.path.split(input_path)[-1], "*.fa"),
                    ),
                )
                
            if script == "prepare":
                from . import prepare

                if not prepare.main(this_args):
                    print("Error in Prepare.")
                gc.collect()
            elif script == "diamond":
                from . import diamond

                if not diamond.main(this_args):
                    print("Error in Diamond.")
                gc.collect()
            elif script == "hmmsearch":
                from . import hmmsearch

                if not hmmsearch.main(this_args):
                    print("Error in Hmmsearch.")
                gc.collect()
            elif script == "reporter":
                from . import reporter

                if not reporter.main(this_args):
                    print("Error in Reporter.")
                gc.collect()
            elif script == "align":
                from . import align

                if not align.main(this_args):
                    print("Error in Align.")
                gc.collect()
            elif script == "pal2nal":
                from . import pal2nal

                if not pal2nal.main(this_args):
                    print("Error in Pal2Nal.")
                gc.collect()
            elif script == "flexcull":
                from . import flexcull

                if not flexcull.main(this_args):
                    print("Error in FlexCull.")
                gc.collect()
            # elif script == "motif":
            #     if args.skip_motif:
            #         continue
            #     from . import motif

            #     if not motif.main(this_args):
            #         print("Error in Motif.")
            #     gc.collect()
            elif script == "outlier":
                from . import outlier

                if not outlier.main(this_args):
                    print("Error in Outlier.")
                gc.collect()
            elif script == "merge":
                from . import consensus_merge

                if not consensus_merge.main(this_args):
                    print("Error in Merge.")
                gc.collect()

    time = time.differential()
    print(f"Took {time:.2f}s overall.")

    return True
