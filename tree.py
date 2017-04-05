import argparse
import string
from collections import defaultdict
from freader import read_nex


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', dest='debug', action='store_true')
    parser.add_argument('-i', '--input', dest='input', action='store', type=argparse.FileType('r'))
    args = parser.parse_args()
    if args.debug:
        matrix = [[0, 1, 2, 3, 4], [1, 0, 3, 4, 5], [2, 3, 0, 5, 6], [3, 4, 5, 0, 7], [4, 5, 6, 7, 0]]
        tree = Tree(matrix=matrix)
        tree.seek_join()
        #print(tree.newick)
        print(tree.names)
        print(tree.newick_format())
    else:
        tree = Tree(reads=read_nex(args.input))
        tree.seek_join()
        print(tree.newick_format())

class Tree(object):
    def __init__(self, matrix=[], reads=[]):
        if matrix:
            self.matrix = matrix
        else:
            self.matrix = self.build_matrix(reads)
        self.newick = {}
        self.names = list(string.ascii_letters[:len(self.matrix)])
        self.extra_names = list(string.ascii_letters[len(self.matrix):])
        self.first_three = []

    def build_matrix(self, reads):
        size = len(reads)
        matrix = [[0 for i in range(size)] for i in range(size)]
        for i in range(size):
            for j in range(size):
                if i == j:
                    res = 0
                else:
                    res = self._find_s(reads[i], reads[j])
                matrix[i][j] = res
        return matrix

    def _find_s(self, read1, read2):
        s = 0
        for base1, base2 in zip(read1, read2):
            if base1 != base2:
                s += 1
        return s

    def _opt(self):
        """ returns the closest two peaks, which are the furthest from others """
        pairs = []
        edge = len(self.matrix)
        peaks = range(edge)
        for A in peaks:
            MA = sum(self.matrix[A])  # not to count this twice
            for B in range(A+1, edge):
                MB = sum(self.matrix[B])
                DAB = self.matrix[A][B]
                pairs.append((DAB / (edge - 2) - MA - MB, (A, B), (MA, MB)))
        return min(pairs)

    def seek_join(self):
        edge = len(self.matrix)
        peaks = range(edge)
        if edge == 3:  # terminate
            self._terminate_task()
            return None
        opt = self._opt()
        self.Ma, self.Mb = opt[2]  # sum of all s between A and others, B and others
        self.A, self.B = opt[1]  # A and B indices
        self.DA = self.matrix.pop(self.A)  # s between A and others
        self.DB = self.matrix.pop(self.B-1)
        self.DUA = 0.5 * (self.DA[self.B] + (self.Ma-self.Mb) / (edge-2))  # s between U and self.A
        self.DUB = 0.5 * (self.DA[self.B] + (self.Mb-self.Ma) / (edge-2))
        print(self.A, self.B, self.Ma, self.Mb, self.DUA, self.DUB, self.DA[self.B], self.DB[self.A])
        self._names_magick()
        self.matrix.append([0 for i in peaks])
        edge -= 1
        peaks = range(edge)
        for X in peaks:
            self._update(X)
        self.seek_join()

    def _update(self, X):
        """ updating rows """
        self.matrix[X].append(0)
        print(self.DA[X], self.DB[X])
        DUX = 0.5 * (self.DA[X] + self.DB[X] - self.DA[self.B])
        self.matrix[X][-1] = DUX
        self.matrix[-1][X] = DUX
        self.matrix[X].pop(self.A)
        self.matrix[X].pop(self.B - 1)

    def _names_magick(self):
        """ names for newick """
        new_name = self.extra_names.pop(0)
        self.newick[new_name] = {self.names[self.A]: self.DUA, self.names[self.B]: self.DUB}
        self.names.pop(self.A)
        self.names.pop(self.B - 1)  # because we poped A so all indices deceased by 1, B is always > A
        self.names.append(new_name)

    def _terminate_task(self):
        s = self.matrix[0][1], self.matrix[0][2], self.matrix[1][2]
        sum_of_s = sum(s)
        for i in range(len(self.matrix)):
            self.first_three.append(sum_of_s - s[i])

    def newick_format(self):
        res = '('
        for peak in range(3):
            res += '%s:%d,' % (self.names[peak], self.first_three[peak])
        res = res[:-1] + ');'
        while self.newick:
            for letter in res:
                if letter in self.newick:
                    letter_childs = self.newick[letter]
                    peaks = list(letter_childs.keys())
                    child_peak = '(%s:%d,%s:%d)' % (peaks[0], letter_childs[peaks[0]], peaks[1], letter_childs[peaks[1]])
                    res = res.replace(letter, child_peak)
                    self.newick.pop(letter)
        return res


if __name__ == '__main__':
    main()
