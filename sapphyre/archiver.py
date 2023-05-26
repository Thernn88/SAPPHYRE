import os
import tarfile
from multiprocessing.pool import Pool
from shutil import rmtree

from .timekeeper import KeeperMode, TimeKeeper
from .utils import printv


def archive_worker(folder_to_archive, verbosity) -> None:
    if os.path.isdir(folder_to_archive) or os.path.isfile(folder_to_archive):
        printv(f"Archiving {folder_to_archive}", verbosity, 2)
        with tarfile.open(str(folder_to_archive) + ".tar.gz", "w:gz") as tf:
            tf.add(str(folder_to_archive), arcname=os.path.split(folder_to_archive)[-1])

        if os.path.isdir(folder_to_archive):
            rmtree(folder_to_archive)
        else:
            os.remove(folder_to_archive)


def unarchive_worker(file_to_unarchive, verbosity) -> None:
    if os.path.exists(file_to_unarchive):
        printv(f"Unarchiving {file_to_unarchive}", verbosity, 2)
        with tarfile.open(file_to_unarchive) as tf:
            tf.extractall(os.path.split(file_to_unarchive)[0])


def process_folder(args, superfolder_path):
    if not args.specific_directories:
        directories_to_archive = [
            "align",
            "nt_aligned",
            "outlier",
            "trimmed",
            "aa_merged",
            "nt_merged",
            "hmmsearch",
            "rocksdb",
            "blast",
            "aa",
            "nt",
        ]
    else:
        directories_to_archive = args.specific_directories

    if args.unarchive:
        directories_to_archive = [i + ".tar.gz" for i in directories_to_archive]

    arguments = []
    for root, _, files in os.walk(superfolder_path):
        #  unarchive logic
        if args.unarchive:
            for file in files:
                if file in directories_to_archive:
                    arguments.append(
                        (
                            os.path.join(root, file),
                            args.verbose,
                        ),
                    )
        #  archive logic
        elif os.path.basename(root) in directories_to_archive:
            if os.path.basename(os.path.split(root)[0]) in directories_to_archive:
                continue
            if (
                os.path.basename(os.path.split(root)[0]) == "sequences"
            ):  # rocksdb/sequences
                continue

            arguments.append(
                (
                    root,
                    args.verbose,
                ),
            )
        else:
            for file in files:
                if file in directories_to_archive:
                    arguments.append(
                        (
                            os.path.join(root, file),
                            args.verbose,
                        ),
                    )
    return arguments


def main(args):
    global_time = TimeKeeper(KeeperMode.DIRECT)
    if not all(os.path.exists(i) for i in args.INPUT):
        printv("ERROR: All folders passed as argument must exists.", args.verbose, 0)
        return False

    arguments = []
    for input_path in args.INPUT:
        arguments.extend(process_folder(args, input_path))
    command = unarchive_worker if args.unarchive else archive_worker
    if args.unarchive:
        printv(f"Found {len(arguments)} directories to unarchive", args.verbose)
    else:
        printv(f"Found {len(arguments)} directories to archive", args.verbose)
    with Pool(args.processes) as pool:
        pool.starmap(command, arguments)

    printv(f"Done! Took {global_time.differential():.2f}s overall.", args.verbose, 0)
    return True


if __name__ == "__main__":
    msg = "Cannot be called directly, please use the module:\nsapphyre archiver"
    raise Exception(
        msg,
    )
