from src.align.approximate.query_suffix_tree import QuerySuffixTree
from src.utils.fasta import FastaContent
from src.utils.sliding_hash_frame import SlidingFrameHasher as SlidingHasher
from src.align.approximate.global_settings import GlobalSettings as Settings

from copy import deepcopy


class ApproximatelyAlignedSequences:
    def __init__(self, query_tree: QuerySuffixTree, target: FastaContent.FastaSequence):
        self.__chunk_matches = dict()
        self.__compare_seq_with_tree(target, query_tree)

    def __compare_seq_with_tree(self, target: FastaContent.FastaSequence, query_tree: QuerySuffixTree) -> None:
        _hash_path = [SlidingHasher(target[_i * Settings.CHUNK_LEN: (_i + 1) * Settings.CHUNK_LEN])
                      for _i in range(Settings.TREE_DEPTH)]
        self.__handle_path(_hash_path, query_tree, 0, target.description)
        for _i in range(Settings.CHUNK_LEN, len(target) - (Settings.TREE_DEPTH - 1) * Settings.CHUNK_LEN):
            if _i % 100000 == 0:
                print(_i // 100000, 'e+5...', sep='', end='')  # TODO Delete after debugging
            for _j in range(len(_hash_path)):
                _hash_path[_j].slide(target[_i + _j * Settings.CHUNK_LEN])
            self.__handle_path(_hash_path, query_tree, _i - Settings.CHUNK_LEN + 1, target.description)
        print()

        target.reverse()

        _hash_path = [SlidingHasher(target[_i * Settings.CHUNK_LEN: (_i + 1) * Settings.CHUNK_LEN])
                      for _i in range(Settings.TREE_DEPTH)]
        self.__handle_path(_hash_path, query_tree, 0, target.description)
        for _i in range(Settings.CHUNK_LEN, len(target) - (Settings.TREE_DEPTH - 1) * Settings.CHUNK_LEN):
            if _i % 100000 == 0:
                print(_i // 100000, 'e+5...', sep='', end='')  # TODO Delete after debugging
            for _j in range(len(_hash_path)):
                _hash_path[_j].slide(target[_i + _j * Settings.CHUNK_LEN])
            self.__handle_path(_hash_path, query_tree, len(target) - (_i - Settings.CHUNK_LEN + 1), target.description)
        print()

    def __handle_path(self, target_path: list, query_tree: QuerySuffixTree, target_entry_index: int, description: str):
        _leaf = query_tree.get_leaf([_p.get_value() for _p in target_path]) + \
                query_tree.get_leaf([_p.get_complimentary() for _p in target_path])
        if _leaf:
            if target_entry_index not in self.__chunk_matches.keys():
                self.__chunk_matches[target_entry_index] = list()
            for _query_entry_index_and_description in _leaf:
                self.__chunk_matches[target_entry_index].append(deepcopy(_query_entry_index_and_description))

    def keys(self):
        return list(self.__chunk_matches.keys())

    def __getitem__(self, item):
        return self.__chunk_matches[item]

    def __str__(self):
        return str(self.__chunk_matches)
