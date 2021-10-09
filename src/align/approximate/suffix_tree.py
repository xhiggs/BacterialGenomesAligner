from src.utils.fasta import FastaContent
from src.align.approximate.align_global_settings import AlignGlobalSettings as Settings
from src.utils.sliding_hash_frame import SlidingFrameHasher


class SuffixTree:
    def __init__(self, sequence: FastaContent.FastaSequence):
        self.__tree = {}
        self.__build(sequence)
        print('Suffix tree is built!')  # TODO Delete after debugging

    def __build(self, sequence: FastaContent.FastaSequence):
        _descent_stack = [SlidingFrameHasher(sequence[_i * Settings.CHUNK_LEN:(_i + 1) * Settings.CHUNK_LEN])
                          for _i in range(Settings.TREE_DEPTH)]
        self.__add_descent(_descent_stack, self.__tree, 0)
        for _i in range(Settings.CHUNK_LEN, len(sequence) - (Settings.TREE_DEPTH - 1) * Settings.CHUNK_LEN):
            if _i % 100000 == 0:
                print(_i)  # TODO Delete after debugging
            for _j in range(len(_descent_stack)):
                _descent_stack[_j].slide(sequence[_i + _j * Settings.CHUNK_LEN])
            self.__add_descent(_descent_stack, self.__tree, _i - Settings.CHUNK_LEN + 1)

    def __add_descent(self, descent_stack: list, tree_node: dict, entry_index: int, index=0) -> None:
        if index + 1 == len(descent_stack):
            if descent_stack[index].get_value() not in tree_node.keys():
                tree_node[descent_stack[index].get_value()] = []
            tree_node[descent_stack[index].get_value()].append(entry_index)
        elif index + 1 < len(descent_stack):
            if descent_stack[index].get_value() not in tree_node.keys():
                tree_node[descent_stack[index].get_value()] = {}
            tree_node = tree_node[descent_stack[index].get_value()]
            self.__add_descent(descent_stack, tree_node, entry_index, index + 1)

    def get_leaf(self, path: list) -> list:
        _current_node = self.__tree
        for _i in range(len(path)):
            if path[_i] not in _current_node.keys():
                return []
            else:
                _current_node = _current_node[path[_i]]
        return _current_node
