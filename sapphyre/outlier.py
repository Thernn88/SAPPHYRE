import argparse
import os
from . import blosum, collapser, excise, internal
from .utils import printv
from .timekeeper import TimeKeeper, KeeperMode

def main(argsobj):
    timer = TimeKeeper(KeeperMode.DIRECT)
    to_move = []
    debug = 0 if argsobj.debug is None else argsobj.debug
    if debug > 1:
        printv("Debug Mode Enabled, Skipping final Excise.", argsobj.verbose, 0)
    for folder in argsobj.INPUT:
        if not os.path.exists(folder):
            printv(
                "ERROR: All folders passed as argument must exists.", argsobj.verbose, 0
            )
            return

        printv(f"Processing: {folder}", argsobj.verbose)
        printv("Blosum62 Outlier Removal.", argsobj.verbose)
        is_assembly = False
        this_args = vars(argsobj)
        this_args["INPUT"] = folder
        this_args = argparse.Namespace(**this_args)

        blosum_passed, is_assembly = blosum.main(this_args)
        if not blosum_passed:
            print()
            print(argsobj.formathelp())
            return        

        if not argsobj.no_excise:
            printv("Checking for severe contamination.", argsobj.verbose)
            module_return_tuple = excise.main(this_args, False)
            if not module_return_tuple:
                print()
                print(argsobj.formathelp())

            is_flagged = module_return_tuple[1]
            if is_flagged:
                bad_folder = os.path.join(this_args.move_fails, os.path.basename(folder))
                to_move.append((folder, bad_folder))

        if not is_assembly:
            printv("Removing Gross Consensus Disagreements.", argsobj.verbose)
            if not internal.main(this_args, False):
                print()
                print(argsobj.format)

        printv("Simple Assembly To Ensure Consistency.", argsobj.verbose)
        if not collapser.main(this_args):
            print()
            print(argsobj.formathelp())
            return

        if debug > 1:
            continue

        if not argsobj.no_excise:
            printv("Detecting and Removing Ambiguous Regions.", argsobj.verbose)
            module_return_tuple = excise.main(this_args, True)
            if not module_return_tuple:
                print()
                print(argsobj.formathelp())

        if is_assembly:
            printv("Removing Gross Consensus Disagreements.", argsobj.verbose)
            if not internal.main(this_args, True):
                print()
                print(argsobj.format)
    # if to_move:
    #     printv("Moving Pre-flagged Folders.", argsobj.verbose)

    #     excise.move_flagged(to_move, this_args.processes)
    printv(f"Took {timer.differential():.2f} seconds overall.", argsobj.verbose)

    return True

if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre OutlierCheck"
    raise RuntimeError(
        MSG,
    )
