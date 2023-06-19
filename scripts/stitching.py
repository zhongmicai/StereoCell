import os
import sys
_cellbin_ = os.path.join(os.path.split(os.path.abspath(__file__))[0], '..')
sys.path.append(_cellbin_)
import argparse
from cellbin.utils.util import search_files
from cellbin.stitching import Stitcher


def image_map(images):
    """
    Args:
        images: List of images
    Returns: dict -> Used to describe the number of row and column indexes that
            each file belongs to in the scanning space
    """
    dct = dict()
    for i in images:
        file_name = os.path.split(i)[-1]
        tags = os.path.splitext(file_name)[0].split('_')
        xy = list()
        for tag in tags:
            if (len(tag) == 4) and tag.isdigit(): xy.append(tag)
        x_str = xy[1]
        y_str = xy[0]

        dct['{}_{}'.format(y_str.zfill(4), x_str.zfill(4))] = i
    return dct


def main(args, para):
    tiles = search_files(args.input, ['.tif'])
    imap = image_map(tiles)
    s = Stitcher()
    s.stitch(imap, output_path=args.output)


if __name__ == '__main__':
    usage="""Stitching (StereoCell)"""
    PROG_VERSION='v0.0.1'

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("--version", action="version", version=PROG_VERSION)
    parser.add_argument("-i", "--input", action="store", dest="input", type=str, required=True, help="Input image dir.")
    parser.add_argument("-o", "--output", action="store", dest="output", type=str, required=True, help="Result output dir.")
    parser.set_defaults(func=main)
    (para, args) = parser.parse_known_args()
    para.func(para, args)

