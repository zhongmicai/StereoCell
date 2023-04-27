import argparse


def main(args, para):
    import tifffile
    a = tifffile.imread(args.src)
    b = tifffile.imread(args.dst)
    assert a.shape == b.shape
    tifffile.imwrite(args.output, a * b, compression="zlib", compressionargs={"level": 8})


if __name__ == '__main__':
    usage = """Mask merge """
    PROG_VERSION = 'v0.0.1'

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("--version", action="version", version=PROG_VERSION)
    parser.add_argument("-s", "--src", action="store", dest="src", type=str, required=True, help="Input mask1.")
    parser.add_argument("-d", "--dst", action="store", dest="dst", type=str, required=True, help="Input mask2.")
    parser.add_argument("-o", "--output", action="store", dest="output", type=str, required=True,
                        help="Result output file.")
    parser.set_defaults(func=main)
    (para, args) = parser.parse_known_args()
    para.func(para, args)
