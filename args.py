import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', help='Index in sequence to begin search')
    parser.add_argument('-e', help='Index in sequence to end search')
    parser.add_argument('-m', help='Starting motif')
    parser.add_argument('-c', help='Path to clustal formated alignment file')
    parser.add_argument('-o', default='.', help='Output path for png graph. Defauly = .')

    return parser.parse_args()