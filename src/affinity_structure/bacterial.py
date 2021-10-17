from src.utils.global_settings import GlobalSettings as Settings
from src.utils.fasta import FastaContent
from src.align.segmental import SegmentalAlign
from src.align.suffix_tree.query import QuerySuffixTree

import os


class BacterialAffinityStructure:
    def __init__(self, folder_name: str):
        self.__filenames = [folder_name + '/' + _file for _file in os.listdir(Settings.DATA_FOLDER + folder_name)]
        self.__query_tree = QuerySuffixTree()

    def handle_next_genome(self) -> None:
        if self.__query_tree.is_empty():
            self.__build_in_tree_next(FastaContent(self.__filenames.pop(0)))
        _fasta_compared = self.__compare_with_tree_next()
        print('Comparing : {}'.format(_fasta_compared.first_seq.description))  # TODO delete after debugging
        self.__build_in_tree_next(_fasta_compared)
        del _fasta_compared

    def __build_in_tree_next(self, fasta: FastaContent) -> None:
        self.__query_tree.supply(fasta.first_seq)

    def __compare_with_tree_next(self) -> FastaContent:
        _fasta = FastaContent(self.__filenames.pop(0))
        _segmental = SegmentalAlign(_fasta.first_seq, self.__query_tree)
        _segmental.plot()
        del _segmental
        return _fasta
