import debuging
from pprint import pprint
from collections import namedtuple, Counter
from typing import List, NamedTuple, TextIO


class Gene(NamedTuple):
    start: int
    end: int
    strand: str
    name: str

Genes = List[Gene]

class Read(NamedTuple):
    start: int
    end: int
    strand: str


class Quantifier:
    def __init__(self, gtf_filename:str, sam_filename:str):
        self.genes: Genes = self.read_gtf(open(gtf_filename))
        self.genes_counter = Counter()
        self.sam_filename = sam_filename

    def find_equal_genes(self):
        names_counter = Counter()
        for i in self.genes:
            names_counter[i.start] += 1
        print(names_counter.most_common(20))

    def count_reads(self):
        found_index = 0
        for read in self.read_sam(open(self.sam_filename)):
            for i in range(found_index, len(self.genes)):
                gene = self.genes[i]
                if read.strand != gene.strand:
                    continue

                if read.end < gene.start:
                    break
                if (gene.start <= read.start <= gene.end) or (gene.start <= read.end <= gene.end):
                    self.genes_counter[gene.name] += 1
                    found_index = i
                    break

    def print_reads(self, n):
        for i, read in enumerate(self.read_sam(open(self.sam_filename))):
            if i == n:
                break
            print(read)

    def print_genes(self, n):
        for gene in self.genes[:n+1]:
            print(gene)

    @staticmethod
    def read_sam(file):
        for line in file:
            line_list = line.split('\t')
            flag = int(line_list[1])
            if flag & 4:  # if unmapped
                continue

            if flag  & 16:  # if being reverse-complemental
                strand = '-'
            else:
                strand = '+'

            start = int(line_list[3])
            length = len(line_list[9])
            end = start + length

            read = Read(start=start, end=end, strand=strand)
            yield read

    @staticmethod
    def read_gtf(file)->Genes:
        genes: Genes = list()
        for line in file:
            line_list = line.strip().split('\t')
            name = line_list[8].split('; ')[2]
            name = name[name.find('"')+1:name.rfind('"')]
            start = int(line_list[3])
            end = int(line_list[4])
            strand = line_list[6]
            gene: Gene = Gene(start=start, end=end, strand=strand, name=name)
            genes.append(gene)
        return genes


def save_pi(file: TextIO, first: Counter, second: Counter):
    for gene in first.keys():
        file.write(f'{gene}\t{first[gene]}\t{second[gene]}\n')

def main():
    gtf_file = 'DE/genome_annotation.gtf'
    sam_file_thyp = 'DE/data/THYP2_22.sam'
    quant_thyp = Quantifier(gtf_file, sam_file_thyp)
    quant_thyp.count_reads()

    sam_file_tnor = 'DE/data/TNOR2_22.sam'
    quant_tnor = Quantifier(gtf_file, sam_file_tnor)
    quant_tnor.count_reads()
    #quant_tnor.print_reads(100)
    with open('result.pi', 'w') as pi_file:
        save_pi(pi_file, quant_thyp.genes_counter, quant_tnor.genes_counter)


if __name__ == '__main__':
    main()