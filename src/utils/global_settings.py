class GlobalSettings:
    CHUNK_LEN = None
    TREE_DEPTH = None
    TREE_SEQUENCE_LEN = None
    MIN_CONSIDERING_SEGMENT_LEN = None
    SEGMENT_MIN_SIZE = int(1e4)
    SEGMENTS_JOIN_SIZE = int(1e4) + 3
    DOT_SKIP_RATE = 10
    DATA_FOLDER = './data/'

    @staticmethod
    def init(chunk_len=14, tree_depth=7, tree_sequence_len=1500000, min_considering_segment_amount=int(1e2)):
        GlobalSettings.CHUNK_LEN = chunk_len
        GlobalSettings.TREE_DEPTH = tree_depth
        GlobalSettings.TREE_SEQUENCE_LEN = tree_sequence_len
        GlobalSettings.MIN_CONSIDERING_SEGMENT_LEN = min_considering_segment_amount

        GlobalSettings.SEGMENT_MIN_SIZE *= GlobalSettings.CHUNK_LEN
        GlobalSettings.SEGMENTS_JOIN_SIZE *= GlobalSettings.CHUNK_LEN
