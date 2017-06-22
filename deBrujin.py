# coding: utf-8

# работает только на python3.6 (typing)

from typing import Dict, Set, List, Iterable, DefaultDict, NamedTuple
from collections import defaultdict, Counter
from copy import deepcopy
import pickle
import os

from tqdm import tqdm

from freader import read_dna


# Тут тестирование (можно не смотреть)
def test_formater():
    formater = GraphVizFormater('read')
    formater.peak = 'ACGT'
    formater.outs = {'CGTA': Edge('ACGTA', 1)}
    formater.out = 'CGTA'
    
    res_read = formater.format_edge()
    res_read_peak = formater.format_peak()
    print(res_read_peak)
    print(res_read)
    assert res_read_peak == ' ACGT[label="ACGT"]\n'
    assert res_read == '  ACGT -> CGTA [label="ACGTA 1"];\n'
    formater.change_mode('info')
    res_info = formater.format_edge()
    res_info_peak = formater.format_peak()
    print(res_info_peak)
    print(res_info)
    assert res_info_peak == ' ACGT[label=""]\n'
    assert res_info == '  ACGT -> CGTA[label="L = 5 Q = 1"];\n'

def test_graph():
    reads = ['ABCDA', 'CDABD', 'DABCC']
    k = 3
    graph = Graph(reads, k)
    graph.cut_all()
    graph.make_graph()
    assert graph.graph == {'ABC': {'BCC': Edge(edge='ABCC', count=1), 'BCD': Edge(edge='ABCD', count=1)},
                           'DAB': {'ABD': Edge(edge='DABD', count=1), 'ABC': Edge(edge='DABC', count=1)},
                           'BCD': {'CDA': Edge(edge='BCDA', count=1)}, 'CDA': {'DAB': Edge(edge='CDAB', count=1)}}
    print(graph.graph)
    graph.graph_viz('deBrujin_test.dot', 'info')
    graph.save_to_pdf('deBrujin_test1.pdf', 'deBrujin_test.dot')

def test():
    test_formater()
    test_graph()


class GraphVizFormater:  # В принципе можно не смотреть - просто красиво сохраняет в .dot
    """Format strings into graphViz form."""

    def __init__(self, mode='read'):
        self.modes = {
            'read': (self._format_as_read, self._format_peak_as_read),
            'info': (self._format_as_info, self._format_peak_as_blank)
        }
        self.change_mode(mode)
        
        self.max_peak_len = 10
        self.peak: str = ''
        self.outs: Dict[str, Edge] = {}
        self.out: str = ''
    
    def change_mode(self, mode):
        if mode not in self.modes:
            raise ValueError('Unknown format mode. Known are: info, read')
        self.format_edge, self.format_peak = self.modes[mode]
    
    def _format_peak_as_read(self) -> str:
        """Format peak with its sequence displayed in max_peak_len."""
        return f' {self.peak}[label="{self.peak[:min(len(self.peak), self.max_peak_len)]}"]\n'
    
    def _format_peak_as_blank(self) -> str:
        """Format peak as blank."""
        return f' {self.peak}[label=""]\n'
    
    def _format_as_read(self) -> str:
        """Format an edge that connects two peaks.
        Peaks and the edge are displayed with their code.
        """
        edge = self.outs[self.out]
        frame = f'  {self.peak} -> {self.out} [label="{edge.edge} {edge.count}"];\n'
        return frame
    
    def _format_as_info(self) -> str:
        """Format an edge that connects two peaks.
        Peaks are blank, edge is displayed with its length and count.
        """
        edge = self.outs[self.out]
        edge_chars = f'"L = {len(edge.edge)} Q = {edge.count}"'
        frame = f'  {self.peak} -> {self.out}[label={edge_chars}];\n'
        return frame


class Edge(NamedTuple):
    edge: str
    count: int


class Graph:
    """De Brujin Graph implementation
    Minimal order to create a graph from reads:
      1) cut_all
      2) make_graph
    """
    
    def __init__(self, reads: Iterable[str], k: int):
        self.edges: Dict[str, int] = Counter()  # keys - edges, values - edge count  
        self.peaks: Set[str] = set()  # Множество вершин графа
        self.k: int = k  # kmer (peak) size
        self.k_edge: int = self.k + 1  # edge size
        self.reads: Iterable[str] = reads  # Это итератор, который возвращает риды (чтобы создать граф)
        self.graph: DefaultDict[Dict[str, Edge]] = defaultdict(dict)
        self.reversed_graph: DefaultDict[List[str]] = defaultdict(list)

    def _cut(self, read):
        """Cut read into kmers"""
        for i in range(len(read) - self.k_edge + 1):
            self.edges[read[i:i + self.k_edge]] += 1
            self.peaks.add(read[i:i + self.k])
        self.peaks.add(read[-self.k:])

    def cut_all(self):
        for read in self.reads:
            self._cut(read) 
    
    def make_graph(self):
        """Constructs the graph"""
        self.graph = defaultdict(dict)
        for peak in tqdm(self.peaks):
            for another_peak in self.peaks:
                if peak[1:] == another_peak[:-1]:
                    edge = peak + another_peak[-1]
                    if edge in self.edges:
                        self.graph[peak][another_peak] = Edge(edge, self.edges[edge])
                        self.reversed_graph[peak].append(peak)

    def graph_viz(self, filename, mode='read'):
        """Format graph into graphViz form with mode and save into file"""
        formater = GraphVizFormater(mode)
        with open(filename, 'w') as file:
            file.write('digraph A{\n')
            for peak in self.peaks.difference(set(self.graph.keys())):  # specify label for peaks without edges from them
                formater.peak = peak
                file.write(formater.format_peak())
            file.write('\n')
            for peak, outs in self.graph.items():
                formater.peak = peak
                file.write(formater.format_peak())
                formater.outs = outs
                for out in outs:
                    formater.out = out
                    file.write(formater.format_edge())
            file.write('}')
        return filename
    
    def save_to_pdf(self, pdf: str, dot: str = None, mode: str = 'read'):
        """Makes pdf file with graph from dot file."""
        print(f'Making pdf {pdf}...')
        if dot is None:
            dot = self.graph_viz('tmp.dot', mode)
        os.popen(f'dot {dot} -T pdf -o {pdf}')
        #os.popen('rm tmp.dot')

    def zip_(self):
        """Example of zipping a graph:
             / AGC     GCT   \         /  AGCT    \
            A----->B-------->C   =>   A---------->C
            \   3        2    \       \     2      \
        A and C are not allowed to be connected before zipping
        B connects only A and C.
        """
        used: List = []  # removed peaks
        for peak, outs in self.graph.items():
            ins = self.reversed_graph[peak]  # list of zpeaks connected to past peak
            if (len(outs) == 1 and len(ins) == 1):  # if peak is connected only to two other peaks
                print(peak)
                past_peak: str = ins[0]
                next_peak: str = list(outs.keys())[0]
                next_ins: List[str] = self.reversed_graph[next_peak]
                past_outs: List[str] = list(self.graph[past_peak].keys())

                if not ((past_peak in next_ins) or (next_peak in past_outs)):  # if not in a loop
                    out_edge: Edge = outs[next_peak]
                    average_coverage = (out_edge.count + list(self.graph[past].values())[0].count) // 2
                    new_edge = past_peak[0] + out_edge.read
                    used.append(peak)
                    self.graph[past_peak].pop(peak)
                    self.graph[past_peak][next_peak] = Edge(new_edge, average_coverage)
                    self.reversed_graph[next_peak].remove(peak)
                    self.reversed_graph[next_peak].append(past_peak)

        for peak in used:
            print(peak)
            self.graph.pop(peak)
            self.reversed_graph.pop(peak)

test_file = open('test.dna')
k = 7  # must be odd
graph = Graph(read_dna(test_file), k)
graph.cut_all()
print(len(graph.peaks))

graph.make_graph()
graph.graph_viz('test_alien.dot', 'info')
#graph.save_to_pdf('test_alien.pdf', 'test_alien.dot')

graph_ = deepcopy(graph.graph)
reversed_graph_ = deepcopy(graph.reversed_graph)

print(len(graph.graph))
outs_c = Counter()
for peak, outs in graph.graph.items():
    outs_c[len(outs)] += 1

print(outs_c.most_common()[::-1])

graph.zip()
graph.graph_viz('test_alien_zipped.dot', 'info')

count = 0
for peak, outs in graph.graph.items():
    if len(outs) == 1:
        count += 1
print(count)

pickle.dump((graph_, reversed_graph_), open('graph_k6.p', 'wb'))

file = open('test.dna')
k = 6
g = Graph(read_dna(file), k)
g.cut_all()
print(len(g.peaks))

g.graph = deepcopy(graph_)
g.reversed_graph = deepcopy(reversed_graph_)

g.zip_()

print(len(g.graph))

del g
