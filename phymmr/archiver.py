import os
from multiprocessing.pool import Pool
import tarfile
from pathlib import Path
from shutil import rmtree
from .utils import printv
from .timekeeper import TimeKeeper, KeeperMode

def archive_worker(folder_to_archive: Path, verbosity) -> None:
    if folder_to_archive.is_dir():
        printv(f"Archiving {folder_to_archive}", verbosity, 2)
        with tarfile.open(str(folder_to_archive) + ".tar.gz", "w:gz") as tf:
            tf.add(str(folder_to_archive), arcname=folder_to_archive.parts[-1])
            
        rmtree(folder_to_archive)

def unarchive_worker(file_to_unarchive: Path, verbosity) -> None:
    if file_to_unarchive.exists():
        printv(f"Unarchiving {file_to_unarchive}", verbosity, 2)
        with tarfile.open(str(file_to_unarchive)) as tf:
            tf.extractall(str(file_to_unarchive.parent))

def process_folder(args, superfolder_path):
    tk = TimeKeeper(KeeperMode.DIRECT)
    if not args.specific_directories:
        directories_to_archive = ['*/aa', '*/nt', '*/mafft', '*/nt_aligned','*/outlier','*/trimmed','*/aa_merged','*/nt_merged','*/hmmsearch','*/rocksdb','*/blast']
    else:
        directories_to_archive = args.specific_directories

    if args.unarchive:
        directories_to_archive = [i+'.tar.gz' for i in directories_to_archive]
    
    arguments = []
    for gp in directories_to_archive:
        for folder in superfolder_path.glob(gp):
            arguments.append((folder, args.verbose))

    if args.unarchive:
        printv(f"Found {len(arguments)} directories to unarchive", args.verbose)
    else:
        printv(f"Found {len(arguments)} directories to archive", args.verbose)

    command = unarchive_worker if args.unarchive else archive_worker
    with Pool(args.processes) as pool:
        pool.starmap(command, arguments)
    
    printv(f"Done! Took {tk.lap():.2f}s overall.", args.verbose)
            
def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False
        
    for input_path in args.INPUT:
        process_folder(args, Path(input_path))
    if len(args.INPUT) > 1 or not args.verbose:
        printv(f"Done! Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True

if __name__ == "__main__":
    raise Exception(
        "Cannot be called directly, please use the module:\nphymmr toolbox"
    )
