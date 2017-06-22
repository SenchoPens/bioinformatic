"""
Multiple Local Alignement
"""

import copy
import typing
from collections import defaultdict, Counter
from pprint import pprint
from typing import List, Dict, Tuple, Union
from debuging import perform, echo
from biocolor import color_nucleotids


CostMatrix = List[List[Union[int, float]]]
Pos = Tuple[int, int]
WayMatrix = List[List[Union[Pos, int]]]
Seq = List[str]
PNV = List[List[int]]


def test_multi():
    first = 'AC'
    a = MultiSequence(first, glob=True, debug=True, fancy=True)
    a.add('AGC')
    #a.add('ACTA')
    a.view()
    print(color_nucleotids(a.consensus()))


class MultiSequence(object):
    def __init__(self, start:str, glob:bool=False, debug:bool=False, fancy:bool=False, indel=-1.5, mismatch=-2, match=2):
        ### initializing variables
        self.matrix: CostMatrix = list()
        self.way: WayMatrix = list()
        self.max_pos: Pos = tuple()
        self.max_score: int = 0
        self.fancy = fancy

        self.match: int = match
        self.indel: int = indel
        self.mismatch: int = mismatch
        self.glob: bool = glob
        self.debug = debug

        self.debug_perform = perform(lambda: self.debug)
        self.debug_print = echo(lambda: self.debug)

        self.seq: str = str()
        self.cols: int = int()
        self.cons: List[typing.Counter[str, int]] = list()

        self.aligned: List[str] = list()

        ### actions and action-dependant variables
        self.make_cons(start)
        self.rows: int = len(self.cons) + 1

    def _fill_default(self):
        self.cons.append(Counter())
        for nucl in 'ACGT':
            self.cons[-1][nucl] = self.mismatch
        self.cons[-1]['-'] = self.indel

    def make_cons(self, start):
        self.aligned.append(start)
        for i, n in enumerate(start):
            self._fill_default()
            self.cons[-1][n] = self.match

    def _debug_view(self):
        for i in self.cons:
            print(dict(i), end=' ')
        print()

    def view(self):
        print_func = lambda x: print(color_nucleotids(x)) if self.fancy else print
        if self.debug:
            self._debug_view()
        for seq in self.aligned:
            print(color_nucleotids(seq))

    def add(self, seq: str):
        self.seq = seq
        self.cols = len(self.seq) + 1
        self.make_matrix()
        self.build()

    def make_matrix(self):
        """ construct score and way matrix to operate with it """
        self.matrix = [[0 for j in range(self.cols)] for i in range(self.rows)]
        self.way = copy.deepcopy(self.matrix)
        for i in range(1, self.rows):  # filling
            for j in range(1, self.cols):
                self._fill(i, j)


    def _fill_borders(self):
        for i in range(self.cols):
            self.matrix[0][i] = self.indel * i
        for i in range(self.rows):
            self.matrix[i][0] = self.indel * i

    def _fill(self, i, j):
        """ fill the M[i, j] cell with correct score and append its parent into way matrix """
        diag = self.matrix[i-1][j-1] + self.cons[i-1][self.seq[j-1]], (i-1, j-1)
        if len(self.cons) < j:
            self._fill_default()
        up = self.matrix[i-1][j] + self.cons[j-1]['-'], (i-1, j)
        left = self.matrix[i][j-1] + self.cons[i-1]['-'], (i, j-1)
        optimal = max(diag, up, left)

        if not self.glob and optimal[0] <= 0:
            optimal = (0, optimal[1])
        opt_score = optimal[0]
        self.matrix[i][j] = opt_score
        opt_pos = optimal[1]
        self.way[i][j] = opt_pos
        if opt_score > self.max_score:
            self.max_score, self.max_pos = optimal

    def build(self):
        """ give feedback to consensus and align sequence """
        seq_aligned = []
        if self.glob:
            move: Pos = (self.rows - 1, self.cols - 1)
        else:
            move: Pos = self.max_pos
        while move != 0:
            row_coord = move[0]
            col_coord = move[1]
            next_move = self.way[row_coord][col_coord]
            if next_move != 0:
                direction: Pos = (row_coord - next_move[0], col_coord - next_move[1])
                if direction == (1, 1):  # diagonal
                    self.cons[row_coord-1][self.seq[col_coord-1]] += self.match
                    seq_aligned.append(self.seq[col_coord-1])
                elif direction == (1, 0):  # up
                    self.cons[row_coord-1]['-'] += self.match
                    seq_aligned.append('-')
                elif direction == (0, 1):  # left
                    self.cons[col_coord-1]['-'] += self.match
                    seq_aligned.append(self.seq[col_coord-1])
            move = next_move
        seq_aligned.reverse()
        self.aligned.append(''.join(seq_aligned))

    def consensus(self):
        result = ''
        for i in self.cons:
            result += i.most_common(1)[0][0]
            self.debug_print(i.most_common())
        return result

    def make_pnv(self)->PNV:
        raise NotImplementedError
        pnv: PNV = list()
        nucl_numbers = {}
        for i in self.aligned:
            pnv
        return pnv

    @staticmethod
    def entropy(column: str):
        import math
        p: List[float] = list()
        for i in 'ACGT-':
            p.append(column.count(i) / len(column))
        res = 0
        for pi in p:
            if pi:
                res += pi * math.log(pi)
        return 1 - res / math.log(1 / 5)


if __name__ == '__main__':
    test_multi()
