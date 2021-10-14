from src.align.approximate.align_global_settings import AlignGlobalSettings
from src.align.segmental_aligned_sequences import SegmentalAlignedSequences


AlignGlobalSettings.init()

sas = SegmentalAlignedSequences('data/large2/large_genome1.fasta', 'data/large2/large_genome2.fasta')
sas.plot()
