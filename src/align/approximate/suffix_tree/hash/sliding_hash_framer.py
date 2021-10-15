from src.utils.global_settings import GlobalSettings as Settings


class SlidingHashFramer:
    class HashFrame:
        def __init__(self, initial_subsequence: str):
            if len(initial_subsequence) == Settings.CHUNK_LEN:
                self.value = sum(
                    [self.__nucleo_to_num(initial_subsequence[i]) * (4 ** i) for i in range(len(initial_subsequence))])
            else:
                raise Exception('Cannot hash init subsequence which len is not as in AlignGlobalSettings')

        def slide_by_value(self, add_value: int) -> int:
            _remainder = self.value % 4
            self.value //= 4
            self.value += add_value * (4 ** (Settings.CHUNK_LEN - 1))
            return _remainder

        def slide_by_nucleo(self, add_nucleo: str) -> int:
            return self.slide_by_value(self.__nucleo_to_num(add_nucleo))

        @property
        def complimentary(self) -> int:
            return 4 ** Settings.CHUNK_LEN - 1 - self.value

        @staticmethod
        def __nucleo_to_num(nucleo: str) -> int:
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

    def __init__(self, initial_subsequence: str):
        if len(initial_subsequence) == Settings.CHUNK_LEN * Settings.TREE_DEPTH:
            self.__hash_path = [
                self.HashFrame(initial_subsequence[_i * Settings.CHUNK_LEN:(_i + 1) * Settings.CHUNK_LEN]) for _i in
                range(Settings.TREE_DEPTH)]
        else:
            raise Exception('Cannot hash init subsequence which len is not as in AlignGlobalSettings')

    def slide_by_nucleo(self, add_nucleo: str):
        _value_additive = self.__hash_path[-1].slide_by_nucleo(add_nucleo)
        for _i in range(len(self.__hash_path) - 2, -1, -1):
            _value_additive = self.__hash_path[_i].slide_by_value(_value_additive)

    def slide_by_frame(self, initial_subsequence: str):
        self.__hash_path.pop(0)
        self.__hash_path.append(self.HashFrame(initial_subsequence))

    @property
    def values(self):
        return [_h.value for _h in self.__hash_path]

    @property
    def complimentary(self):
        return [_h.complimentary for _h in self.__hash_path]

    def __getitem__(self, item):
        return self.__hash_path[item].value

    def __len__(self):
        return len(self.__hash_path)
