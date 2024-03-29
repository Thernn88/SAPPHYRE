import argparse
import json
import os
from glob import glob

from .__main__ import (
    align_args,
    diamond_args,
    flexcull_args,
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
    default_config = {
        "prepare": get_args(prepare_args),
        "diamond": get_args(diamond_args),
        "hmmsearch": get_args(hmmsearch_args),
        "reporter": get_args(reporter_args),
        "align": get_args(align_args),
        "pal2nal": get_args(pal2nal_args),
        "flexcull": get_args(flexcull_args),
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
        sargs = config[script]
        print(f"\nExecuting: {script.title()}")
        this_args = global_args.copy()
        this_args.update(sargs)
        
        if args.overwrite:
            if 'overwrite' in this_args:
                this_args['overwrite'] = True

        this_args = argparse.Namespace(**this_args)
        if args.solo:
            this_args.INPUT = [args.INPUT]
        elif script == "prepare":
            this_args.INPUT = args.INPUT
        else:
            this_args.INPUT = sorted(
                glob(
                    os.path.join(args.in_folder, os.path.split(args.INPUT)[-1], "*.fa"),
                ),
            )

        if script == "prepare":
            from . import prepare

            if not prepare.main(this_args):
                print("Error in Prepare.")
        elif script == "diamond":
            from . import diamond

            if not diamond.main(this_args):
                print("Error in Diamond.")
        elif script == "hmmsearch":
            from . import hmmsearch

            if not hmmsearch.main(this_args):
                print("Error in Hmmsearch.")
        elif script == "reporter":
            from . import reporter

            if not reporter.main(this_args):
                print("Error in Reporter.")
        elif script == "align":
            from . import align

            if not align.main(this_args):
                print("Error in Align.")
        elif script == "pal2nal":
            from . import pal2nal

            if not pal2nal.main(this_args):
                print("Error in Pal2Nal.")
        elif script == "flexcull":
            if args.map:
                continue
            from . import flexcull

            if not flexcull.main(this_args):
                print("Error in FlexCull.")
        elif script == "outlier":
            from . import outlier

            if not outlier.main(this_args):
                print("Error in Outlier.")
        elif script == "merge":
            from . import consensus_merge

            if not consensus_merge.main(this_args):
                print("Error in Merge.")

    time = time.differential()
    print(f"Took {time:.2f}s overall.")

    return True
