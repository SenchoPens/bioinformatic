from freader import read_fq
from sequence import reverse_complemental
from collections import defaultdict
import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input_file', type=argparse.FileType('r'))
    parser.add_argument('-d', '--debug', action='store_true', dest='debug')
    parser.add_argument('-g', '--graph', action='store', dest='graph_viz')
    parser.add_argument('-o', '--output-graph-file', action='store', dest='output')
    parser.add_argument('-k', '--kmers', action='store', dest='k')
    return parser.parse_args()


def main():
    args = parse_args()
    if args.debug:
        k = 3
        reads2 = ['ACGGT', 'CGGTTT', 'GTTTTA', 'TTACGC', 'GCTTA', 'ACGCT']
        reads = ['ACTGA']
    else:
        reads = read_fq(args.input_file)
        k = int(args.k)
    graph = Graph(reads, k)
    graph.graph()
    graph.zip()
    if args.graph_viz:
        graph.graph_viz(args.graph_viz)
        if args.output:
            os.popen('dot %s -T %s -o %s' % (args.graph_viz, 'pdf', args.output))


class Graph(object):
    def __init__(self, reads, k):
        self.edges = defaultdict(lambda: 0)
        self.peaks = set()
        self.k = k
        self.k_edge = self.k + 1
        self.reads = reads

    def _process(self, read):
        for i in range(len(read) - self.k_edge + 1):
            self.edges[read[i:i+self.k_edge]] += 1
            self.peaks.add(read[i:i+self.k])
        self.peaks.add(read[-self.k:])

    def graph2(self):
        for read in self.reads:
            self._process(read)
            self._process(reverse_complemental(read))
        self.g = defaultdict(dict)
        self.g_reverse = defaultdict(list)
        for peak in self.peaks:
            for another_peak in self.peaks:
                if peak[1:] == another_peak[:-1]:
                    edge = peak + another_peak[-1]
                    if edge in self.edges:
                        self.g[peak][another_peak] = (edge, self.edges[edge])
                        self.g_reverse[another_peak].append(peak)

    def graph(self):
        for read in self.reads:
            self._process(read)
            self._process(reverse_complemental(read))
        self.g = defaultdict(dict)
        self.g_reverse = defaultdict(list)
        for peak in self.peaks:
            for another_peak in self.peaks:
                if peak[1:] == another_peak[:-1]:
                    edge = peak + another_peak[-1]
                    if edge in self.edges:
                        self.g[peak][another_peak] = (edge, self.edges[edge])
                        self.g_reverse[another_peak].append(peak)

    def graph_viz2(self, file_name):
        frame = ' %s -> %s [label="%s %d"];\n'
        with open(file_name, 'w') as file:
            file.write('digraph A{\n')
            for peak, outs in self.g.items():
                for out in outs:
                    file.write(frame % (peak, out, outs[out][0], outs[out][1]))
            file.write('}')

    def graph_viz(self, file_name):
        frame = ' %s[label=""]%s[label=""]\n %s -> %s[label=%s];\n'
        edge_chars_frame = '"L = %d Q = %d"'
        with open(file_name, 'w') as file:
            file.write('digraph A{\n')
            for peak, outs in self.g.items():
                for out in outs:
                    edge_chars = edge_chars_frame % (len(outs[out][0]), outs[out][1])
                    file.write(frame % (peak, out, peak, out, edge_chars))
            file.write('}')

    def zip(self):
        """ example of zipping a graph:
             / AGC     GCT   \         /  AGCT    \
            A----->B-------->C   =>   A---------->C
            \   3        2    \       \     2      \
        A and C are not allowed to connect
        B connects only to A and C
        """
        used = []
        for peak, outs in self.g.items():
            reversed_peak = self.g_reverse[peak]  # list of peaks, that are connected to past peak
            if (len(outs) == 1 and len(reversed_peak) == 1):  # if peak is connected only to two other peaks
                past = reversed_peak[0]  # because reversed peak is connected only to one past el, this el is past, str
                next = list(outs.keys())[0]  # the same as past, str
                next_ins = self.g_reverse[next]  # reversed connections of the next peak, list
                past_outs = list(self.g[past].keys())  # connections of past read, list
                if not ((past in next_ins) or (next in past_outs)):  # if not in a loop
                    out_edge = outs[next]  # tuple of edge code and its coverage
                    average_coverage = (out_edge[1] + list(self.g[past].values())[0][1]) // 2
                    new_edge = past[0] + out_edge[0]
                    used.append(peak)
                    self.g[past].pop(peak)
                    self.g[past][next] = (new_edge, average_coverage)
                    self.g_reverse[next].remove(peak)
                    self.g_reverse[next].append(past)
        for peak in used:
            self.g.pop(peak)
            self.g_reverse.pop(peak)


if __name__ == '__main__':
    main()
