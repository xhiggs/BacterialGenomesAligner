from src.align.approximate.global_settings import GlobalSettings
# from src.align.segmental_aligned_sequences import SegmentalAlignedSequences
from src.align.approximate.approximately_aligned_sequences import ApproximatelyAlignedSequences
from src.align.approximate.query_suffix_tree import QuerySuffixTree
from src.align.segmental_aligned_sequences import SegmentalAlignedSequences
from src.utils.fasta import FastaContent

import os
from matplotlib import pyplot as plt


GlobalSettings.init()

# sas = SegmentalAlignedSequences('data/large2/large_genome1.fasta', 'data/large2/large_genome2.fasta')
# sas.plot()


qsf = QuerySuffixTree()
qsf.insert(FastaContent('data/grouptest/LR595848.1.fasta')[0])
qsf.insert(FastaContent('data/grouptest/LR129840.1.fasta')[0])
print(qsf.get_descriptions())
aas = ApproximatelyAlignedSequences(qsf, FastaContent('data/grouptest/CP053210.1.fasta')[0])

x1, x2, y1, y2 = [], [], [], []

for k in aas.keys():
    for v in aas[k]:
        ([x1, x2][qsf.get_descriptions().index(v[1])]).append(k)
        ([y1, y2][qsf.get_descriptions().index(v[1])]).append(v[0])

for x, y in [[x1, y1], [x2, y2]]:
    plt.plot(x, y, '.')
    plt.grid()
    plt.show()
