from sequence import Sequence
from freader import read_fq
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='append', dest='ifiles')
    parser.add_argument('-o', '--output', action='store', dest='ofile')
    args = parser.parse_args()
    res = open(args.ofile, 'w')
    fst = read_fq(open(args.ifiles[0], 'r'))
    sec = read_fq(open(args.ifiles[1], 'r'))
    for pair in zip(fst, sec):
        seq = Sequence(*pair)
        print(seq.result)
        res.write(seq.result)


def test():
    Sequence('AGCACACA', 'ACACACTA')


if __name__ == '__main__':
    main()
