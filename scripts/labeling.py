import os
import sys
import argparse
_cellbin_ = os.path.join(os.path.split(os.path.abspath(__file__))[0], '..')
sys.path.append(_cellbin_)
from cellbin.cell_labeling.GMMCorrectForv03 import CellCorrection


def args_parse():
    usage = """ Usage: %s Cell expression file (with background) path, multi-process """
    arg = argparse.ArgumentParser(usage=usage)
    arg.add_argument('-i', '--image_file', help='cell mask')
    arg.add_argument('-g', '--gene_exp_data', help='gem file')
    arg.add_argument('-o', '--output_path', help='output path', default='./')
    arg.add_argument('-p', '--process', help='n process', type=int, default=10)
    arg.add_argument('-t', '--threshold', help='threshold', type=int, default=20)

    return arg.parse_args()


def main():
    args = args_parse()
    correction = CellCorrection(args.image_file, args.gene_exp_data, args.output_path, args.threshold, args.process)
    correction.cell_correct()


if __name__ == '__main__':
    main()
