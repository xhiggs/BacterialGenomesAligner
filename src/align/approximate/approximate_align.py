from src.align.approximate.suffix_tree.query_suffix_tree import QuerySuffixTree
from src.utils.fasta import FastaContent
from src.align.approximate.suffix_tree.hash.sliding_hash_framer import SlidingHashFramer as HashFramer
from src.utils.global_settings import GlobalSettings as Settings


class ApproximateAlign:
    def __init__(self, target: FastaContent.FastaSequence, query_tree: QuerySuffixTree):
        self.__chunk_matches = {_description: dict() for _description in query_tree.descriptions}
        self.__compare_seq_with_tree(target, query_tree)

    def __compare_seq_with_tree(self, target: FastaContent.FastaSequence, query_tree: QuerySuffixTree) -> None:
        for _reversed in [False, True]:
            _hash_framer = HashFramer(target[:Settings.CHUNK_LEN * Settings.TREE_DEPTH])
            self.__handle_hash_framer(_hash_framer, 0, query_tree)
            for _i in range(Settings.CHUNK_LEN * Settings.TREE_DEPTH, len(target)):
                if _i % 100000 == 0:
                    print(f'{_i // 100000}e+5..', end='')
                _hash_framer.slide_by_nucleo(target[_i])
                _entry_index = _i - Settings.CHUNK_LEN * Settings.TREE_DEPTH + 1
                if _reversed:
                    _entry_index = len(target) - 1 - _entry_index
                    _entry_index *= -1
                self.__handle_hash_framer(_hash_framer, _entry_index, query_tree)
            del _hash_framer
            target.reverse()
            print()

    def __handle_hash_framer(self, hash_framer: HashFramer, entry_index: int, query_tree: QuerySuffixTree):
        for _hash_path in [hash_framer.values, hash_framer.complimentary]:
            _query_leaf = query_tree.get_leaf(_hash_path)
            if _query_leaf:
                for _leaf_description in _query_leaf.keys():
                    if entry_index not in self.__chunk_matches[_leaf_description].keys():
                        self.__chunk_matches[_leaf_description][entry_index] = list()
                    self.__chunk_matches[_leaf_description][entry_index] += _query_leaf[_leaf_description]

    @property
    def keys(self) -> list:
        return list(self.__chunk_matches.keys())

    def __getitem__(self, item) -> dict:
        return self.__chunk_matches[item]
