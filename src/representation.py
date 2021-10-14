from src.align.approximate.global_settings import GlobalSettings
from src.align.segmental_aligned_sequences import SegmentalAlignedSequences


GlobalSettings.init(tree_depth=4, chunk_len=8, min_considering_segment_len=5)

sas = SegmentalAlignedSequences('data/large5/large_genome1.fasta', 'data/large5/large_genome2.fasta')
sas.plot()
