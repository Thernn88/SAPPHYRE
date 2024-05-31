import argparse
import os
from shutil import rmtree
from wrap_rocks import RocksDB
from . import blosum, excise, hmmfilter, internal, cluster_consensus
from .timekeeper import KeeperMode, TimeKeeper
from .utils import printv


def get_genome_assembly(folder):
    rocks_db_path = os.path.join(folder, "rocksdb", "sequences", "nt")
    rocksdb_db = RocksDB(str(rocks_db_path))
    is_assembly = rocksdb_db.get("get:isassembly")
    is_assembly = is_assembly == "True"
    is_genome = rocksdb_db.get("get:isgenome")
    is_genome = is_genome == "True"
    del rocksdb_db
    return is_genome, is_assembly

def main(argsobj):
    timer = TimeKeeper(KeeperMode.DIRECT)
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
        
        this_args = vars(argsobj)
        this_args["INPUT"] = folder
        this_args = argparse.Namespace(**this_args)

        is_genome, is_assembly = get_genome_assembly(folder)

        from_folder = "trimmed"

        printv("Blosum62 Outlier Removal.", argsobj.verbose)
        blosum_passed = blosum.main(this_args, is_genome, from_folder)
        if not blosum_passed:
            print()
            print(argsobj.formathelp())
            return
        from_folder = "blosum"

        if argsobj.map:
            continue

        
        # printv("Simple Assembly To Ensure Consistency.", argsobj.verbose)
        # if not collapser.main(this_args, from_folder):
        #     print()
        #     print(argsobj.formathelp())
        #     return
        # continue

        if is_assembly:
            printv("Filtering Using Hmmsearch Scores.", argsobj.verbose)
            if not hmmfilter.main(this_args, from_folder):
                print()
                print(argsobj.formathelp())
                return
            from_folder = "hmmfilter"    

        if not is_genome:
            printv("Detecting and Removing Ambiguous Regions.", argsobj.verbose)
            excise_passed = excise.main(this_args, from_folder, is_genome, is_assembly or is_genome)
            if not excise_passed:
                print()
                print(argsobj.formathelp())
            from_folder = "excise"

        if not argsobj.gene_finding_mode and is_genome:
            printv("Filtering Clusters.", argsobj.verbose)
            if not cluster_consensus.main(this_args, from_folder):
                print()
                print(argsobj.formathelp())
                return
            from_folder = "clusters"

        if debug > 1 or this_args.add_hmmfilter_dupes:
            continue
        if not (is_assembly or is_genome):
            printv("Removing Gross Consensus Disagreements.", argsobj.verbose)
            if not internal.main(this_args, from_folder):
                print()
                print(argsobj.format)
            from_folder = "internal"
            
        if is_genome:
            printv("Detecting and Removing Ambiguous Regions.", argsobj.verbose)
            excise_passed = excise.main(this_args, from_folder, is_genome, is_assembly or is_genome)
            if not excise_passed:
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
