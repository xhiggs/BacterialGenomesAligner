from src.align.approximate.align_global_settings import AlignGlobalSettings
from src.align.approximate.approximately_aligned_sequences import ApproximatelyAlignedSequences
from src.utils.fasta import FastaContent as Fasta
from matplotlib import pyplot as plt

AlignGlobalSettings.init(tree_depth=3)

x, y = [], []

query, target = Fasta('data/large7/large_genome1.fasta')[0], Fasta('data/large7/large_genome2.fasta')[0]
AAS = ApproximatelyAlignedSequences(query, target)
del query, target

print('Approximation is done!')

for k in AAS.keys():
    for v in AAS[k]:
        y.append(abs(k))
        x.append(abs(v))

del AAS

plt.plot(x, y, '.')
plt.grid()
plt.show()
