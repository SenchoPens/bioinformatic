from sequence import Sequence
from freader import read_fa, read_fq
import argparse
import collections


def main():
    reads, sequences, result = make_args()
    sequences_mismatches = dict.fromkeys(sequences)
    for sequence in sequences_mismatches:
        sequences_mismatches[sequence] = collections.defaultdict(lambda: 0)
    print(len(sequences))
    for read in reads:
        scores = []
        for seq in sequences:
            product = Sequence(seq, read)
            product.make_lightweight_matrix()
            scores.append(product)
            print(product.max_score)
        max_scoring_seq = max(scores)
        product = max_scoring_seq
        product.make_matrix()
        product.build()
        product.mismatches()
        for key, value in product.mismatches.items():
            sequences_mismatches[product.a][key] += value
            table_row = '%s    %d    %d' % (sequences.index(product.a), key, value)
            print(table_row)
            result.write(table_row)

def make_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reads', action='store', dest='reads_file', type=argparse.FileType('r'))
    parser.add_argument('-s', '--sequence', action='store', dest='sequences_file', type=argparse.FileType('r'))
    parser.add_argument('-o', '--output', action='store', dest='output_file', type=argparse.FileType('w'))
    parser.add_argument('-d', '--debug', action='store_true', dest='debug')
    args = parser.parse_args()
    if args.debug:
        return [['ACCTG', 'TACCAAC'], ['ACCTATTGCAAC'], open('restest', 'w')]
    return read_fq(args.reads_file), read_fa(args.sequences_file), args.output_file


if __name__ == '__main__':
    main()
