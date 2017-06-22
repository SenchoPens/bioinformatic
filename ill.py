# coding: utf-8

import pandas as pd
import numpy as np
from math import sqrt
from collections import Counter

n_case = 21
n_control = 23

df = pd.read_csv('DE/norm_camp.tsv', sep='\t')
genes = df['gene_symbol']
df.drop('gene_symbol', axis=1, inplace=True)

column_names = list(df)

control_names = column_names[:n_control]
case_names = column_names[n_control:]

ill = Counter()
healthy = Counter()


def t(ser):
    n_case = 21
    n_control = 23
    
    controls = ser[control_names]

    cases = ser[case_names]

    xA = cases.mean()
    xB = controls.mean()

    sA = cases.var()
    sB = controls.var()

    spooled = ((n_case - 1)*sA + (n_control - 1)*sB) / (n_case + n_control - 2)
    if spooled == 0:
        if xA == xB:
            return 0
        return float('inf')
    stat = abs(xA - xB) / (sqrt(spooled)*sqrt(1/n_case + 1/n_control))
    return stat


for i in range(299, 300):
    ser = df.loc[i]
    
    controls = ser[control_names]
    cases = ser[case_names]
    xA = cases.mean()
    xB = controls.mean()

    stat = t(ser)
    pcount = 0
    for _ in range(100):
        rand_ser = pd.Series(np.random.permutation(ser), index=column_names)
        rand_t = t(rand_ser)
        pcount += rand_t >= stat
    p = pcount / 100
    if p >= 0.05:
        r = xA - xB
        expr = abs(r)
        if r < 0:
            healthy[genes[i]] = expr
        else:
            ill[genes[i]] = expr


with open('ill.tsv', 'w') as file:
    for gene, i in ill.most_common():
        file.write(f'{gene}\t{i}\n')

with open('healthy.tsv', 'w') as file:
    for gene, i in healthy.most_common():
        file.write(f'{gene}\t{i}\n')
