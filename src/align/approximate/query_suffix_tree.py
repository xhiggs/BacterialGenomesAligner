from src.utils.fasta import FastaContent
from src.align.approximate.global_settings import GlobalSettings as Settings
from src.utils.sliding_hash_frame import SlidingFrameHasher as Hasher


class QuerySuffixTree:
    def __init__(self):
        self.__tree = {}
        self.__descriptions = []

    def insert(self, sequence: FastaContent.FastaSequence):
        self.__descriptions.append(sequence.description)
        _descent_stack = list()
        for _i in range(len(sequence) // Settings.CHUNK_LEN):
            if _i * Settings.CHUNK_LEN % 100000 == 0:
                print(_i * Settings.CHUNK_LEN // 100000, 'e+5...', sep='', end='')  # TODO Delete after debugging
            _descent_stack.append(Hasher(sequence[_i * Settings.CHUNK_LEN: (_i + 1) * Settings.CHUNK_LEN]))
            if len(_descent_stack) == Settings.TREE_DEPTH:
                self.__create_empty_leaf(_descent_stack, self.__tree) \
                    .append(tuple([(_i - Settings.TREE_DEPTH + 1) * Settings.CHUNK_LEN, sequence.description]))
                _descent_stack.pop(0)
        print('\nSuffix tree is built!')  # TODO Delete after debugging

    def __create_empty_leaf(self, descent_stack: list, tree_node: dict, _stack_index=0) -> list:
        descent_stack_value = descent_stack[_stack_index].get_value()
        if _stack_index + 1 == len(descent_stack):
            if descent_stack_value not in tree_node.keys():
                tree_node[descent_stack_value] = []
            return tree_node[descent_stack_value]
        elif _stack_index + 1 < len(descent_stack):
            if descent_stack_value not in tree_node.keys():
                tree_node[descent_stack_value] = {}
            tree_node = tree_node[descent_stack_value]
            return self.__create_empty_leaf(descent_stack, tree_node, _stack_index + 1)
        else:
            raise Exception()

    def get_leaf(self, path: list) -> list:
        _current_node = self.__tree
        for _i in range(len(path)):
            if path[_i] not in _current_node.keys():
                return []
            else:
                _current_node = _current_node[path[_i]]
        return _current_node

    def create_leaf(self, path: list, entry_index: int, description: str):
        self.__create_empty_leaf(path, self.__tree).append(tuple([entry_index, description]))

    def __repr__(self):
        return str(self.__tree)

    def get_descriptions(self):
        return self.__descriptions[:]
