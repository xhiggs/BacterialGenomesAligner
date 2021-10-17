class GlobalSettings:
    CHUNK_LEN = None
    TREE_DEPTH = None
    MIN_CONSIDERING_SEGMENT_LEN = None
    SEGMENT_MIN_SIZE = None
    SEGMENTS_JOIN_SIZE = None
    DOT_SKIP_RATE = 1
    DATA_FOLDER = './data/'

    @staticmethod
    def init(chunk_len=14, tree_depth=7, min_considering_segment_amount=int(1e2), segment_min_size=int(1e4)):
        GlobalSettings.CHUNK_LEN = chunk_len
        GlobalSettings.TREE_DEPTH = tree_depth
        GlobalSettings.MIN_CONSIDERING_SEGMENT_LEN = min_considering_segment_amount

        GlobalSettings.SEGMENT_MIN_SIZE = segment_min_size * GlobalSettings.CHUNK_LEN
        GlobalSettings.SEGMENTS_JOIN_SIZE = (segment_min_size + 3) * GlobalSettings.CHUNK_LEN
