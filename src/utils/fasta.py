from src.utils.global_settings import GlobalSettings

from copy import copy


class FastaContent:
    class FastaSequence:
        def __init__(self, _description: str, _sequence: str):
            self.__description = _description
            self.__sequence = _sequence

        def reverse(self):
            self.__sequence = self.__sequence[::-1]

        @property
        def description(self):
            return copy(self.__description)

        def __getitem__(self, item): return self.__sequence[item]

        def __len__(self): return len(self.__sequence)

    def __init__(self, _filepath: str):
        self.__sequences = list()
        with open(GlobalSettings.DATA_FOLDER + _filepath, 'r') as _f:
            for _raw_seq in ''.join(_f.readlines()).split('>')[1:]:
                _splitted = _raw_seq.split('\n')
                self.__sequences.append(FastaContent.FastaSequence(_splitted[0], ''.join(_splitted[1:])))

    def __getitem__(self, item): return self.__sequences[item]

    def __len__(self): return len(self.__sequences)
