import argparse
from . import blosum, collapser, excise, internal
import os
from .timekeeper import TimeKeeper, KeeperMode
from .utils import printv

def main(argsobj):
    timer = TimeKeeper(KeeperMode.DIRECT)
    to_move = []
    for folder in argsobj.INPUT:
        if not os.path.exists(folder):
            printv(
                "ERROR: All folders passed as argument must exists.", argsobj.verbose, 0
            )
            return

        printv(f"Processing: {folder}", argsobj.verbose)
        printv("Blosum62 Outlier Removal.", argsobj.verbose)
        this_args = vars(argsobj)
        this_args["INPUT"] = folder
        this_args = argparse.Namespace(**this_args)

        if not blosum.main(this_args):
            print()
            print(argsobj.formathelp())
            return

        printv("Checking for severe contamination.", argsobj.verbose)
        module_return_tuple = excise.main(this_args, False)
        if not module_return_tuple:
            print()
            print(argsobj.formathelp())

        is_flagged = module_return_tuple[1]
        if is_flagged:
            bad_folder = os.path.join(this_args.move_fails, os.path.basename(folder))
            to_move.append((folder, bad_folder))

        printv("Simple Assembly To Ensure Consistency.", argsobj.verbose)
        if not collapser.main(this_args):
            print()
            print(argsobj.formathelp())
            return

        printv("Detecting and Removing Ambiguous Regions.", argsobj.verbose)
        module_return_tuple = excise.main(this_args, True)
        if not module_return_tuple:
            print()
            print(argsobj.formathelp())

        printv("Removing Gross Consensus Disagreements.", argsobj.verbose)
        if not internal.main(this_args):
            print()
            print(argsobj.format)
    if to_move:
        printv("Moving Pre-flagged Folders.", argsobj.verbose)

        excise.move_flagged(to_move, this_args.processes)
        printv(f"Took {timer.differential():.2f} seconds overall.", argsobj.verbose)

    return True

if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre OutlierCheck"
    raise RuntimeError(
        MSG,
    )