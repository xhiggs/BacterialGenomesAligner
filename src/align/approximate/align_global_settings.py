class AlignGlobalSettings:
    CHUNK_LEN = None
    TREE_DEPTH = None
    TREE_SEQUENCE_LEN = None
    MIN_SEGMENT_LEN = None
    SEGMENTS_JOIN_SIZE = None
    SEGMENT_MIN_SIZE = None
    DOT_SKIP_RATE = None

    @staticmethod
    def init(chunk_len=14, tree_depth=7, tree_sequence_len=1500000):
        AlignGlobalSettings.CHUNK_LEN = chunk_len
        AlignGlobalSettings.TREE_DEPTH = tree_depth
        AlignGlobalSettings.TREE_SEQUENCE_LEN = tree_sequence_len
        AlignGlobalSettings.MIN_SEGMENT_LEN = int(1e3)
        AlignGlobalSettings.SEGMENT_MIN_SIZE = int(1e4)
        AlignGlobalSettings.SEGMENTS_JOIN_SIZE = AlignGlobalSettings.SEGMENT_MIN_SIZE + 3
        AlignGlobalSettings.DOT_SKIP_RATE = 10
