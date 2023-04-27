import sys
import os

import glog

_cellbin_ = os.path.join(os.path.split(os.path.abspath(__file__))[0], '..')
sys.path.append(_cellbin_)
import argparse
from cellbin.cell_segmentation.segment import cell_seg
from cellbin.tissue_segmentation import tissue_cut


def main(args, para):
    if args.type == 'tissue': tissue_cut(args.input, args.output)
    elif args.type == 'cell': cell_seg(args.input, args.output, flag=1)
    else: glog.info('unknown segment type')


'''
python .\segmentation.py  --type tissue --input D:\data\weights\135_D1.tif --output D:\data\weights\135_D1_mask.tif
'''
if __name__ == '__main__':
    usage = """Cell Segmentation (StereoCell)"""
    PROG_VERSION = 'v0.0.1'

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("--version", action="version", version=PROG_VERSION)
    parser.add_argument("-i", "--input", action="store", dest="input", type=str, required=True, help="Input image dir.")
    parser.add_argument("-o", "--output", action="store", dest="output", type=str, required=True,
                        help="Result output dir.")
    parser.add_argument("-t", "--type", action="store", dest="type", type=str, required=True,
                        help="Seg type.")
    parser.set_defaults(func=main)
    (para, args) = parser.parse_known_args()
    para.func(para, args)