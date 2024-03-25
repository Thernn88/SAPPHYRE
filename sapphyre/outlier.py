import argparse
import os
from shutil import rmtree

from . import blosum, excise, hmmfilter, internal, collapser
from .timekeeper import KeeperMode, TimeKeeper
from .utils import printv


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
        rmtree(os.path.join(folder, "outlier"), ignore_errors=True)
        printv("Blosum62 Outlier Removal.", argsobj.verbose)
        this_args = vars(argsobj)
        this_args["INPUT"] = folder
        this_args = argparse.Namespace(**this_args)

        blosum_passed, is_assembly, is_genome = blosum.main(this_args)
        if not blosum_passed:
            print()
            print(argsobj.formathelp())
            return
        from_folder = "blosum"

        if argsobj.map:
            continue

        if is_genome:
            printv("Simple Assembly To Ensure Consistency.", argsobj.verbose)
            if not collapser.main(this_args, from_folder):
                print()
                print(argsobj.formathelp())
                return
            continue

        printv("Filtering Using Hmmsearch Scores.", argsobj.verbose)
        if not hmmfilter.main(this_args, from_folder):
            print()
            print(argsobj.formathelp())
            return
        from_folder = "hmmfilter"

        if debug > 1 or this_args.add_hmmfilter_dupes:
            continue

        printv("Removing Gross Consensus Disagreements.", argsobj.verbose)
        if not internal.main(this_args, True, from_folder):
            print()
            print(argsobj.format)
        from_folder = "internal"

        if not argsobj.no_excise:
            printv("Detecting and Removing Ambiguous Regions.", argsobj.verbose)
            module_return_tuple = excise.main(this_args, True, from_folder)
            if not module_return_tuple:
                print()
                print(argsobj.formathelp())
            from_folder = "excise"

    printv(f"Took {timer.differential():.2f} seconds overall.", argsobj.verbose)

    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre OutlierCheck"
    raise RuntimeError(
        MSG,
    )
