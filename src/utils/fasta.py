class FastaContent:
    class FastaSequence:
        def __init__(self, _description: str, _sequence: str):
            self.__description = _description
            self.__sequence = _sequence

        def __getitem__(self, item): return self.__sequence[item]

        def __len__(self): return len(self.__sequence)

        def __str__(self): return self.__description

    def __init__(self, _filepath: str):
        self.__sequences = []
        with open(_filepath, 'r') as _f:
            for _raw_seq in ''.join(_f.readlines()).split('>')[1:]:
                _splitted = _raw_seq.split('\n')
                self.__sequences.append(FastaContent.FastaSequence(
                    _splitted[0], ''.join(_splitted[1:])))

    def __getitem__(self, item): return self.__sequences[item]

    def __len__(self): return len(self.__sequences)

    def __str__(self): return str(list(map(str, self.__sequences)))

    @staticmethod
    def get_instance(filepath: str):
        return FastaContent(filepath)
