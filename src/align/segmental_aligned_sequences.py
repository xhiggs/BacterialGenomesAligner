from src.align.approximate.approximately_aligned_sequences import ApproximatelyAlignedSequences
from src.utils.fasta import FastaContent as Fasta
from src.align.approximate.align_global_settings import AlignGlobalSettings as Settings

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
            if len(_segment.dots) < Settings.MIN_SEGMENT_LEN:
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


class Segment:
    def __init__(self, start_x=None, start_y=None, end_x=None, end_y=None, dots=None):
        self.start_x = start_x
        self.start_y = start_y
        self.end_x = end_x
        self.end_y = end_y
        self.dots = dots if dots is not None else list()

    def __repr__(self):
        return "Segment(start_x={}, start_y={}, end_x={}, end_y={}, dots=[{}])".format(
            self.start_x, self.start_y, self.end_x, self.end_y, len(self.dots))

    @property
    def coords(self):
        return self.start_x, self.start_y, self.end_x, self.end_y

    @property
    def center_x(self):
        return (self.start_x + self.end_x) // 2

    @property
    def center_y(self):
        return (self.start_y + self.end_y) // 2

    @property
    def size_x(self):
        return abs(self.start_x - self.end_x)

    @property
    def size_y(self):
        return abs(self.start_y - self.end_y)

    def is_tilted_correctly(self):
        return self.start_y <= self.end_y

    @property
    def k(self):
        return (self.end_y - self.start_y) / (self.end_x - self.start_x)

    @property
    def b(self):
        return self.end_y - self.end_x * self.k

    def cope_coords(self):
        return Segment(self.start_x, self.start_y, self.end_x, self.end_y, dots=[])

    def shift(self, dx=0, dy=0):
        self.start_x += dx
        self.start_y += dy
        self.end_x += dx
        self.end_y += dy
        for _i in range(len(self.dots)):
            self.dots[_i][0] += dx
            self.dots[_i][1] += dy
        return self

    def rotate_y(self, rotation_center, segment=True, dots=None):
        if segment:
            self.start_y -= (self.start_y - rotation_center) * 2
            self.end_y -= (self.end_y - rotation_center) * 2
        if dots is not None:
            for i in range(len(self.dots)):
                self.dots[i][1] -= (self.dots[i][1] - rotation_center) * 2
        return self

    @staticmethod
    def linear_approx_dots(dots):
        _n, _sum_x, _sum_y, _sum_x2, _sum_xy = len(dots), 0, 0, 0, 0
        for _x, _y in dots:
            _sum_x += _x
            _sum_y += _y
            _sum_x2 += _x ** 2
            _sum_xy += _x * _y

        _k = (_n * _sum_xy - _sum_x * _sum_y) / (_n * _sum_x2 - _sum_x * _sum_x)
        return _k, (_sum_y - _k * _sum_x) / _n

    @staticmethod
    def distance2(x1, y1, x2, y2):
        return (x1 - x2) ** 2 + (y1 - y2) ** 2
