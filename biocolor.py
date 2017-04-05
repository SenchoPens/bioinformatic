from termcolor import colored

def color_nucleotids(seq: str)->str:
    seq_colored = ''
    colors = {'A': 'blue',
              'C': 'red',
              'G': 'green',
              'T': 'yellow',
              '-': 'white'}
    for i in seq:
        seq_colored += colored(i, color=colors[i])
    return seq_colored


if __name__ == '__main__':
    seq = 'ACGT-'
    print(color_nucleotids(seq))
