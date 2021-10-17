from src.utils.global_settings import GlobalSettings
from src.affinity_structure.bacterial import BacterialAffinityStructure


GlobalSettings.init(tree_depth=9, segment_min_size=int(2.5e3))

affinity_structure = BacterialAffinityStructure('grouptest')

for _ in range(2):
    affinity_structure.handle_next_genome()
