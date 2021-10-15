from src.utils.fasta import FastaContent
from src.utils.global_settings import GlobalSettings as Settings
from src.align.approximate.suffix_tree.hash.sliding_hash_framer import SlidingHashFramer as HashFramer

from copy import copy, deepcopy


class QuerySuffixTree:
    def __init__(self):
        self.__tree = dict()
        self.__seqs_descriptions = list()

    def supply(self, sequence: FastaContent.FastaSequence) -> None:
        self.__seqs_descriptions.append(sequence.description)

        _hash_framer = HashFramer(sequence[:Settings.CHUNK_LEN * Settings.TREE_DEPTH])
        self.__update_leaf(_hash_framer.values, 0, sequence.description)

        for _i in range(Settings.CHUNK_LEN * Settings.TREE_DEPTH, len(sequence) - Settings.CHUNK_LEN,
                        Settings.CHUNK_LEN):
            if _i % 100000 == 0:
                print(f'{_i // 100000}e+5..', end='')
            _hash_framer.slide_by_frame(sequence[_i:_i + Settings.CHUNK_LEN])
            self.__update_leaf(
                _hash_framer.values, _i - Settings.CHUNK_LEN * (Settings.TREE_DEPTH - 1), sequence.description)
        print()
        del _hash_framer
        print('Suffix tree is built!')  # TODO Delete after debugging

    def __update_leaf(self, _hash_path: list, entry_index: int, seq_description: str) -> None:
        _current_node = self.__tree
        for _h in _hash_path:
            if _h not in _current_node.keys():
                _current_node[_h] = dict()
            _current_node = _current_node[_h]
        if seq_description not in _current_node.keys():
            _current_node[seq_description] = list()
        if entry_index not in _current_node[seq_description]:
            _current_node[seq_description].append(entry_index)

    def get_leaf(self, _hash_path: list) -> dict:
        _current_node = self.__tree
        for _h in _hash_path:
            if _h not in _current_node.keys():
                return dict()
            else:
                _current_node = _current_node[_h]
        return deepcopy(_current_node)

    @property
    def descriptions(self):
        return copy(self.__seqs_descriptions)

    def __repr__(self):
        return str(self.__tree)
