class AlignGlobalSettings:
    CHUNK_LEN = None
    TREE_DEPTH = None
    TREE_SEQUENCE_LEN = None
    MIN_SEGMENT_LEN = None

    @staticmethod
    def init(chunk_len=14, tree_depth=7, tree_sequence_len=1500000, segment_gap_factor=500):
        AlignGlobalSettings.CHUNK_LEN = chunk_len
        AlignGlobalSettings.TREE_DEPTH = tree_depth
        AlignGlobalSettings.TREE_SEQUENCE_LEN = tree_sequence_len
        AlignGlobalSettings.MIN_SEGMENT_LEN = int(1e3)
