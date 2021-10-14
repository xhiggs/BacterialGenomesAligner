from src.utils.approximately_aligned_sequences import ApproximatelyAlignedSequences
from src.align.approximate.global_settings import GlobalSettings as Settings
from src.utils.fasta import FastaContent as Fasta
from src.utils.segment import Segment

from matplotlib import pyplot as plt


class SegmentalAlignedSequences:
    def __init__(self, query_filepath: str, target_filepath: str):
        self.__segments = list()
        _query_seq, _target_seq = Fasta(query_filepath)[0], Fasta(target_filepath)[0]
        _approx = ApproximatelyAlignedSequences(_query_seq, _target_seq)
        del _query_seq, _target_seq

        _graph = [list() for _ in range(0, max(_approx.keys()) + 1)]
        for _k in _approx.keys():
            for _v in [(1 if _k >= 0 else -1) * _i for _i in _approx[_k]]:
                _graph[abs(_k)].append(_v)
        self.__find_segments(_graph)
        del _approx, _graph

    def __find_segments(self, graph: list) -> None:
        print('Counting ... ', end='')

        _segments_join_size = Settings.SEGMENTS_JOIN_SIZE ** 2
        __segment_min_size = Settings.SEGMENT_MIN_SIZE ** 2

        for _x in range(0, len(graph), Settings.DOT_SKIP_RATE):
            for _y in graph[_x]:
                for _segment in self.__segments:
                    if Segment.distance2(_x, _y, *_segment.dots[-1]) <= _segments_join_size and \
                            (len(_segment.dots) == 1 or Segment.distance2(_x, _y, *_segment.dots[-2]) <=
                             _segments_join_size):
                        _segment.dots.append([_x, _y])
                        break
                else:
                    self.__segments.append(Segment(dots=[[_x, _y]]))

        for _segment in self.__segments:
            _segment.dots.sort()

            _segment.start_x, _segment.start_y = _segment.dots[0]
            _segment.end_x, _segment.end_y = _segment.dots[-1]

            if len(_segment.dots) >= 2:
                k, b = Segment.linear_approx_dots(_segment.dots)  # \
                _segment.start_y = int(k * _segment.start_x + b)  # |--> Approximation  TODO: int
                _segment.end_y = int(k * _segment.end_x + b)  # /
            # _segment[4] = _segment[4][::settings["dot_skip_rate"]]  # Optional compress

        self.__segments = [_segment for _segment in self.__segments if Segment.distance2(
            _segment.start_x, _segment.start_y, _segment.end_x, _segment.end_y) >= __segment_min_size]
        self.__segments.sort(key=lambda _segment: (_segment.start_x, _segment.start_y))

        for _segment in self.__segments:
            print(_segment, len(_segment.dots))
            if len(_segment.dots) < Settings.MIN_CONSIDERING_SEGMENT_LEN:
                self.__segments.remove(_segment)

        print(" {} segments :".format(len(self.__segments)))
        print(*self.__segments, sep='\n')

    def __len__(self) -> int:
        return len(self.__segments)

    def __getitem__(self, item) -> list:
        return self.__segments[item]

    def plot(self):  # TODO delete after debugging
        for _segment in self.__segments:
            plt.plot([_segment.start_x, _segment.end_x], [_segment.start_y, _segment.end_y])
        plt.grid()
        plt.show()
