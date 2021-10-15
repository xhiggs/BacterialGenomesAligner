from src.utils.global_settings import GlobalSettings
from src.align.segmental_align import SegmentalAlign
from src.align.approximate.suffix_tree.query_suffix_tree import QuerySuffixTree
from src.utils.fasta import FastaContent

GlobalSettings.init()

f = FastaContent('large2/large_genome1.fasta')[0]

q = QuerySuffixTree()
q.supply(f)

x, y = [], []

aas = SegmentalAlign(FastaContent('large2/large_genome2.fasta')[0], q)
aas.plot()
