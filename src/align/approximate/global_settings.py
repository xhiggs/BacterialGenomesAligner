class GlobalSettings:
    CHUNK_LEN = None
    TREE_DEPTH = None
    TREE_SEQUENCE_LEN = None
    MIN_SEGMENT_LEN = None
    SEGMENTS_JOIN_SIZE = None
    SEGMENT_MIN_SIZE = None
    DOT_SKIP_RATE = None

    @staticmethod
    def init(chunk_len=14, tree_depth=7, tree_sequence_len=1500000, min_considering_segment_len=int(1e4)):
        GlobalSettings.CHUNK_LEN = chunk_len
        GlobalSettings.TREE_DEPTH = tree_depth
        GlobalSettings.TREE_SEQUENCE_LEN = tree_sequence_len
        GlobalSettings.MIN_CONSIDERING_SEGMENT_LEN = min_considering_segment_len
        GlobalSettings.SEGMENT_MIN_SIZE = int(1e4)
        GlobalSettings.SEGMENTS_JOIN_SIZE = GlobalSettings.SEGMENT_MIN_SIZE + 3
        GlobalSettings.DOT_SKIP_RATE = 10
