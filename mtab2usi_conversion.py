import argparse

from converter.mtab2json import mtab2usi_conversion


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('idf',
                        help="name of MAGE-TAB IDF file")

    args = parser.parse_args()

    return args.idf


def main():
    idf_file = parse_args()
    mtab2usi_conversion(idf_file)


if __name__ == '__main__':
    main()
