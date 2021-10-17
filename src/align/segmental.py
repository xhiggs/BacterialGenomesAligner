from src.align.suffix_tree.query import QuerySuffixTree
from src.align.approximate import ApproximateAlign
from src.utils.global_settings import GlobalSettings as Settings
from src.utils.fasta import FastaContent

from matplotlib import pyplot as plt


class SegmentalAlign:
    class Segment:
        def __init__(self, start_x=None, start_y=None, end_x=None, end_y=None, dots=None):
            self.start_x = start_x
            self.start_y = start_y
            self.end_x = end_x
            self.end_y = end_y
            self.dots = dots if dots is not None else list()

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
            return SegmentalAlign.Segment(self.start_x, self.start_y, self.end_x, self.end_y, dots=[])

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

        def __repr__(self):
            return "Segment(start_x={}, start_y={}, end_x={}, end_y={}, dots=[{}])".format(
                self.start_x, self.start_y, self.end_x, self.end_y, len(self.dots))

    def __init__(self, target_sequence: FastaContent.FastaSequence, query_tree: QuerySuffixTree):
        self.__seqs_segments = dict()
        _approx = ApproximateAlign(target_sequence, query_tree)
        self.__plot_approx(_approx)
        for _seq_d in _approx.keys:
            _graph = [list() for _ in range(len(target_sequence) + 1)]
            for _k in _approx[_seq_d].keys():
                for _v in _approx[_seq_d][_k]:
                    _graph[abs(_k)].append(abs(_v))
            self.__seqs_segments[_seq_d] = self.__find_segments(_graph)
            del _graph
        del _approx

    @staticmethod
    def __find_segments(graph: list) -> list:
        print('Searching for segments ...')

        _all_segments = list()

        _segments_join_size = Settings.SEGMENTS_JOIN_SIZE ** 2
        _segment_min_size = Settings.SEGMENT_MIN_SIZE ** 2

        for _x in range(0, len(graph), Settings.DOT_SKIP_RATE):
            for _y in graph[_x]:
                for _segment in _all_segments:
                    if SegmentalAlign.Segment.distance2(_x, _y, *_segment.dots[-1]) <= _segments_join_size and \
                            (len(_segment.dots) == 1 or SegmentalAlign.Segment.distance2(_x, _y, *_segment.dots[-2]) <=
                             _segments_join_size):
                        _segment.dots.append([_x, _y])
                        break
                else:
                    _all_segments.append(SegmentalAlign.Segment(dots=[[_x, _y]]))

        for _segment in _all_segments:
            _segment.dots.sort()

            _segment.start_x, _segment.start_y = _segment.dots[0]
            _segment.end_x, _segment.end_y = _segment.dots[-1]

            if len(_segment.dots) >= 2:
                k, b = SegmentalAlign.Segment.linear_approx_dots(_segment.dots)  # \
                _segment.start_y = int(k * _segment.start_x + b)  # |--> Approximation  TODO: int
                _segment.end_y = int(k * _segment.end_x + b)  # /
            # _segment[4] = _segment[4][::settings["dot_skip_rate"]]  # Optional compress

        _all_segments = [_segment for _segment in _all_segments if SegmentalAlign.Segment.distance2(
            _segment.start_x, _segment.start_y, _segment.end_x, _segment.end_y) >= _segment_min_size]
        _all_segments.sort(key=lambda _segment: (_segment.start_x, _segment.start_y))

        for _segment in _all_segments:
            # print(_segment, len(_segment.dots))
            if len(_segment.dots) < Settings.MIN_CONSIDERING_SEGMENT_LEN:
                _all_segments.remove(_segment)

        # print(" {} segments :".format(len(_all_segments)))
        # print(*_all_segments, sep='\n')

        return _all_segments

    def __len__(self) -> int:
        return len(self.__seqs_segments.keys())

    def __getitem__(self, item) -> list:
        return self.__seqs_segments[item]

    def plot(self):  # TODO delete after debugging
        for _s_q in self.__seqs_segments.keys():
            plt.figure(figsize=(8, 6))
            for _segment in self.__seqs_segments[_s_q]:
                plt.plot([_segment.start_x, _segment.end_x], [abs(_segment.start_y), abs(_segment.end_y)])
            plt.grid()
            plt.title(_s_q)
            plt.show()

    @staticmethod
    def __plot_approx(approx: ApproximateAlign):
        for _seq_q in approx.keys:
            plt.figure(figsize=(8, 6))
            _x, _y = list(), list()
            for _k in approx[_seq_q]:
                for _v in approx[_seq_q][_k]:
                    _x.append(abs(_k))
                    _y.append(abs(_v))
            plt.plot(_x, _y, '.')
            plt.title(_seq_q)
            plt.grid()
            plt.show()
