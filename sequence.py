import copy


def test():
    a = Sequence('ACCTGACTCCAG', 'CCAGATTGCAA')
    a.make_matrix()
    a.build()
    a.make_cool_matrix()
    print(a.seq1, a.seq2)


class Sequence(object):
    def __init__(self, a, b, match=1, dismatch=-1, indel=-1):
        self.a = a
        self.b = b
        self.match = match
        self.dismatch = dismatch
        self.indel = indel
        self.rows = len(self.a) + 1
        self.cols = len(self.b) + 1

    def make_matrix(self):
        """ construct score and way matrix to operate with it """
        self.matrix = [[0 for i in range(self.cols)] for i in range(self.rows)]
        self.way = copy.deepcopy(self.matrix)
        self.max_score, self.max_pos = 0, None
        for i in range(1, self.rows):  # заполняем
            for j in range(1, self.cols):
                self._fill(i, j)

    def _fill(self, i, j):
        """ fill the M[i, j] cell with correct score and append its parent into way matrix """
        if self.a[i-1] == self.b[j-1]:
            k = self.match
        else:
            k = self.dismatch
        diag = self.matrix[i-1][j-1] + k, (i-1, j-1)
        up = self.matrix[i-1][j] + self.indel, (i-1, j)
        left = self.matrix[i][j-1] + self.indel, (i, j-1)
        optimal = max(diag, up, left)
        if optimal[0] <= 0:
            optimal = (0, 0)
        opt_score = optimal[0]
        self.matrix[i][j] = opt_score
        opt_pos = optimal[1]
        self.way[i][j] = opt_pos
        if optimal[0] > self.max_score:
            self.max_score, self.max_pos = optimal

    def make_lightweight_matrix(self):
        """ an alghoritm just to know the max_score"""
        self.matrix = [[0 for i in range(self.cols)] for i in range(2)]
        self.max_score = 0
        for i in range(1, self.rows):  # заполняем
            for j in range(1, self.cols):
                self._lightweight_fill(i, j)
            self.matrix[1], self.matrix[0] = self.matrix[0], self.matrix[1]

    def _lightweight_fill(self, i, j):
        if self.a[i-1] == self.b[j-1]:
            k = self.match
        else:
            k = self.dismatch
        diag = self.matrix[0][j-1] + k
        up = self.matrix[0][j] + self.indel
        left = self.matrix[1][j-1] + self.indel
        opt_score = max(0, diag, up, left)
        self.matrix[1][j] = opt_score
        if opt_score >= self.max_score:
            self.max_score = opt_score

    def visualize(self):
        """ a debug tool """
        for i in self.matrix:
            print(i)
        print()
        for i in self.way:
            print(i)
        print(self.max_pos)
        print(self.max_score)

    def make_cool_matrix(self):
        """ a way vizualization with cell scores """
        self.cool_matrix = [['✤' for i in range(self.cols)] for i in range(self.rows)]
        for i in range(1, self.rows):
            ci = i - 1
            # self.cool_matrix[ci][0] = '✤'
            for j in range(self.cols):
                next_move = self.way[i][j]
                if next_move == 0 or self.matrix[ci][j] == 0:
                    symbol = '✤'
                else:
                    direction = i - next_move[0], j - next_move[1]
                    if self.way[next_move[0]][next_move[1]] == 0:
                        symbol = '✤'
                    elif direction == (1, 1):  # диагонально
                        symbol = '↖'
                    elif direction == (1, 0):  # вверх
                        symbol = '↑'
                    elif direction == (0, 1):  # влево
                        symbol = '←'
                self.cool_matrix[ci][j - 1] = symbol
        for i in self.matrix[0]:
            print(i, end=' ')
        print()
        for points, directs in zip(self.matrix[1:], self.cool_matrix):
            print('    ', end='')
            for direct in directs[1:]:
                print(direct, end=' ')
            print()
            for point in points:
                print(point, end=' ')
            print()


    def build(self):
        """ find the most optimal local alingment part """
        self.seq1 = []
        self.seq2 = []
        move = self.max_pos
        self.seq2_end = move[1]
        while move != 0:
            row_coord = move[0]
            col_coord = move[1]
            next_move = self.way[row_coord][col_coord]
            if next_move != 0:
                direction = row_coord - next_move[0], col_coord - next_move[1]
                if direction == (1, 1):   # диагонально
                    self.seq1.append(self.a[row_coord - 1])
                    self.seq2.append(self.b[col_coord - 1])
                elif direction == (1, 0):  # вверх
                    self.seq1.append(self.a[row_coord - 1])
                    self.seq2.append('-')
                elif direction == (0, 1):  # влево
                    self.seq1.append('-')
                    self.seq2.append(self.b[col_coord - 1])
            else:
                self.seq1_start = row_coord
            move = next_move
        self.seq1.reverse()
        self.seq2.reverse()

    def concatenate(self):
        self.equal = ''.join(self.seq1)
        self.result = self.a[:self.seq1_start] + self.equal + self.b[self.seq2_end:]
    
    def mismatches(self):
        self.mismatches = dict()
        for i, (nucl1, nucl2) in enumerate(zip(self.seq1, self.seq2)):
            if nucl1 != nucl2 and nucl1 != '-' and nucl2 != '-':
                mismatch_index = i + self.seq1_start - self.seq1[:i].count('-')
                self.mismatches[mismatch_index] = 1

    def __gt__(self, other):
        """ overload to enable max() of sequences by their max_score"""
        if isinstance(other, Sequence):
            return self.max_score > other.max_score
        else:
            raise TypeError('unorderable types: Sequence > %s' % type(other))


def reverse_complemental(read):
    compl = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rc_read = list(map(lambda x: compl[x], read))
    rc_read.reverse()
    rc_read = ''.join(rc_read)
    return rc_read


if __name__ == '__main__':
    test()