from itertools import chain
from pathlib import Path


data_dir = Path('/home/sencho/Projects/bio/vdj_data/')
gene_dir = data_dir / 'Germline/IG/'
seq_dir = data_dir / 'rep_seq/'
dnr_dir = seq_dir / 'Dnr/'
twins_dir = seq_dir / 'Twins/'


def vdj_dir_iter():
    for vdj_dir in chain(dnr_dir.iterdir(), twins_dir.iterdir()):
        yield vdj_dir / 'igr/vjf/'

if __name__ == '__main__':
    print(list(next(vdj_dir_iter()).iterdir()))
