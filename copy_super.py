import argparse
import os
import shutil
from sys import argv


def folder_check(folder: str) -> None:
    if not os.path.exists(folder):
        os.mkdir(folder)


def make_aa_nt(clone: str) -> None:
    folder_check(os.path.join(clone,'aa'))
    folder_check(os.path.join(clone,'nt'))


def get_filenames(folder: str) -> list:
    result = [item for item in os.listdir(folder) if item[0] not in ['.','$']]
    return result


def duplicate_files(input_dir: str, out_dir: str) -> None:
    fastas = get_filenames(input_dir)
    for fasta in fastas:
        in_path = os.path.join(input_dir, fasta)
        out_path = os.path.join(out_dir, fasta)
        shutil.copy(in_path, out_path)


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True, help="Super folder")
    parser.add_argument('-o','--output', required=True,
                        help="Location of cloned directory")
    args = parser.parse_args()

    folder_check(args.output)

    taxas = [taxa for taxa in os.listdir(args.input) if os.path.isdir(os.path.join(args.input, taxa))]
    for taxa in taxas:
        clone = os.path.join(args.output, taxa)
        folder_check(clone)
        make_aa_nt(clone)
        aa_in = os.path.join(os.path.join(args.input, taxa), 'aa')
        nt_in = os.path.join(os.path.join(args.input, taxa), 'nt')
        aa_out = os.path.join(clone, 'aa')
        nt_out = os.path.join(clone, 'nt')
        duplicate_files(aa_in, aa_out)
        duplicate_files(nt_in, nt_out)


if __name__ == '__main__':
    main(argv)
