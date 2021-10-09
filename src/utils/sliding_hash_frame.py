from src.align.approximate.align_global_settings import AlignGlobalSettings


class SlidingFrameHasher:
    def __init__(self, init_subsequence: str):
        if len(init_subsequence) == AlignGlobalSettings.CHUNK_LEN:
            self.__hash_value = sum([self.nucleo_to_num(init_subsequence[i]) * (4 ** i)
                                     for i in range(len(init_subsequence))])
        else:
            raise Exception('Cannot hash init subsequence which len is not as in AlignGlobalSettings')

    def slide(self, adding_nucleo: str) -> None:
        self.__hash_value //= 4
        self.__hash_value += self.nucleo_to_num(adding_nucleo) * (4 ** (AlignGlobalSettings.CHUNK_LEN - 1))

    def get_value(self) -> int:
        return self.__hash_value

    def get_complimentary(self) -> int:
        return 4 ** AlignGlobalSettings.CHUNK_LEN - 1 - self.get_value()

    @staticmethod
    def nucleo_to_num(nucleo: str) -> int:
        if nucleo == 'A':
            return 0
        elif nucleo == 'C':
            return 1
        elif nucleo == 'G':
            return 2
        elif nucleo == 'T':
            return 3
        else:
            raise Exception('Unknown nucleo \'{}\''.format(nucleo))
