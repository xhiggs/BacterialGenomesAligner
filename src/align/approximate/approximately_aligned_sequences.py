from src.utils.fasta import FastaContent
from src.align.approximate.suffix_tree import SuffixTree
from src.utils.sliding_hash_frame import SlidingFrameHasher as Hasher
from src.align.approximate.align_global_settings import AlignGlobalSettings as Settings


class ApproximatelyAlignedSequences:
    def __init__(self, query: FastaContent.FastaSequence, target: FastaContent.FastaSequence):
        self.__chunk_matches = dict()
        for _i in range((len(query) - 1) // Settings.TREE_SEQUENCE_LEN + 1):
            query_tree = SuffixTree(query[_i * Settings.TREE_SEQUENCE_LEN:(_i + 1) * Settings.TREE_SEQUENCE_LEN])
            self.__compare_seq_with_tree(query_tree, target, _i * Settings.TREE_SEQUENCE_LEN)
            del query_tree

    def __compare_seq_with_tree(
            self, query_tree: SuffixTree, target: FastaContent.FastaSequence, additive_index: int) -> None:
        for _increasing in [True, False]:
            _hash_path = list()
            for _i in range(len(target) // Settings.CHUNK_LEN):
                _hash_path.append(Hasher(
                    target[_i * Settings.CHUNK_LEN: (_i + 1) * Settings.CHUNK_LEN] if _increasing else \
                        target[-_i * Settings.CHUNK_LEN - 1: -(_i + 1) * Settings.CHUNK_LEN - 1:-1]))
                if len(_hash_path) == Settings.TREE_DEPTH:
                    _entry_indexes = query_tree.get_leaf([i.get_value() for i in _hash_path]) + \
                                     query_tree.get_leaf([i.get_complimentary() for i in _hash_path])
                    _hash_path.pop(0)
                    if _entry_indexes:
                        self.__chunk_matches[
                            _i if _increasing else -(len(target) // Settings.CHUNK_LEN - 1 - _i)] = \
                            [_k + additive_index for _k in _entry_indexes]

    def __getitem__(self, item):
        return self.__chunk_matches[item]

    def keys(self):
        return list(self.__chunk_matches.keys())

    def __str__(self):
        return str(self.__chunk_matches)
