from src.utils.fasta import FastaContent
from src.align.approximate.suffix_tree import SuffixTree


class ApproximatelyAlignedSequences:
    def __init__(self, query: FastaContent.FastaSequence, target: FastaContent.FastaSequence):
        self.__chunk_matches = dict()
        query_tree, target_tree = SuffixTree(query), SuffixTree(target)
        self.__compare_suffix_trees(query_tree.as_dict(), target_tree.as_dict())

    def __compare_suffix_trees(self, query_tree, target_tree) -> None:
        if isinstance(query_tree, dict) and isinstance(target_tree, dict):
            for _q_key in query_tree.keys():
                if _q_key in target_tree.keys():
                    self.__compare_suffix_trees(query_tree[_q_key], target_tree[_q_key])
        elif isinstance(query_tree, list) and isinstance(target_tree, list):
            for _q in query_tree:
                if _q not in self.__chunk_matches.keys():
                    self.__chunk_matches[_q] = list()
                self.__chunk_matches[_q] += target_tree[:]
        else:
            raise Exception("Cannot compare suffix trees with types '{}' and '{}'".format(
                query_tree.__class__, target_tree.__class__))

    def __str__(self):
        return str(self.__chunk_matches)
