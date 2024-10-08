import argparse
import os
from shutil import rmtree
from sapphyre import excise, recover
from wrap_rocks import RocksDB
from . import blosum, genome_splice, hmmfilter, internal, cluster_consensus, read_assembly
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
        
        if argsobj.solo:
            script_name = argsobj.solo.lower()
            default_path = {
                "blosum": ("trimmed", blosum.main),
                "hmmfilter": ("blosum", hmmfilter.main),
                "clusters": ("blosum", cluster_consensus.main),
                "splice": ("clusters", genome_splice.main),
                "assembly": ("blosum", read_assembly.main),
                "internal": ("excise", internal.main),
                "recover": ("internal", recover.main),
            }
            from_folder, script = default_path[script_name]
            
            if script_name == "blosum":
                script(this_args, is_genome, from_folder)
            else:
                script(this_args, from_folder)
            
        else:
            printv(f"Processing: {folder}", argsobj.verbose)
            rmtree(os.path.join(folder, "outlier"), ignore_errors=True)
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


            if is_assembly:
                printv("Filtering Using Hmmsearch Scores.", argsobj.verbose)
                if not hmmfilter.main(this_args, from_folder):
                    print()
                    print(argsobj.formathelp())
                    return
                from_folder = "hmmfilter"    

            elif is_genome:
                if not argsobj.gene_finding_mode:
                    printv("Filtering Clusters.", argsobj.verbose)
                    if not cluster_consensus.main(this_args, from_folder):
                        print()
                        print(argsobj.formathelp())
                        return
                    from_folder = "clusters"
                    
                printv("Detecting and Splicing Ambiguous Regions in Clusters.", argsobj.verbose)
                excise_passed = genome_splice.main(this_args, from_folder)
                if not excise_passed:
                    print()
                    print(argsobj.formathelp())
                from_folder = "excise"
            else:   
                printv("Detecting and Removing Ambiguous Regions.", argsobj.verbose)
                excise_passed = read_assembly.main(this_args, from_folder)
                if not excise_passed:
                    print()
                    print(argsobj.formathelp())
                from_folder = "excise"
            
                printv("Removing Gross Consensus Disagreements.", argsobj.verbose)
                if not internal.main(this_args, from_folder):
                    print()
                    print(argsobj.format)
                from_folder = "internal"

                printv("Recovering lost reads.", argsobj.verbose)
                if not recover.main(this_args, from_folder):
                    print()
                    print(argsobj.format)
                from_folder = "recovered"

    printv(f"Took {timer.differential():.2f} seconds overall.", argsobj.verbose)

    return True


if __name__ == "__main__":
    MSG = "Cannot be called directly, please use the module:\nsapphyre OutlierCheck"
    raise RuntimeError(
        MSG,
    )
