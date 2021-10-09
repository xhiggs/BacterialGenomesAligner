from src.align.approximate.align_global_settings import AlignGlobalSettings
from src.align.approximate.approximately_aligned_sequences import ApproximatelyAlignedSequences
from src.utils.fasta import FastaContent as FC

AlignGlobalSettings.CHUNK_LEN = 2
AlignGlobalSettings.TREE_DEPTH = 1

query, target = FC('data/debug/1.fasta')[0], FC('data/debug/2.fasta')[0]

print(str(ApproximatelyAlignedSequences(query, target)))
