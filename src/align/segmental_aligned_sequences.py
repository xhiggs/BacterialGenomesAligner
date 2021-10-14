from src.align.approximate.approximately_aligned_sequences import ApproximatelyAlignedSequences
from src.utils.fasta import FastaContent as Fasta
from src.align.approximate.align_global_settings import AlignGlobalSettings as Settings

from matplotlib import pyplot as plt
from copy import deepcopy


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
        print(graph)
        print("Counting lines...", end="")

        lines_join_size2 = settings["lines_join_size"] ** 2
        line_min_size2 = settings["line_min_size"] ** 2

        lines = []

        for x in range(0, len(graph), settings["dot_skip_rate"]):
            for y in graph[x]:
                for line in lines:
                    if distance2(x, y, *line.dots[-1]) <= lines_join_size2 and \
                            (len(line.dots) == 1 or distance2(x, y, *line.dots[-2]) <= lines_join_size2):
                        line.dots.append([x, y])
                        break
                else:
                    lines.append(Line(dots=[[x, y]]))

        for line in lines:
            line.dots.sort()

            line.start_x, line.start_y = line.dots[0]
            line.end_x, line.end_y = line.dots[-1]

            if len(line.dots) >= 2:
                k, b = linearApproxDots(line.dots)  # \
                line.start_y = int(k * line.start_x + b)  # |--> Approximation  TODO: int
                line.end_y = int(k * line.end_x + b)  # /

            # line[4] = line[4][::settings["dot_skip_rate"]]  # Optional compress

        lines = [line for line in lines if
                 distance2(line.start_x, line.start_y, line.end_x, line.end_y) >= line_min_size2]

        lines.sort(key=lambda line: (line.start_x, line.start_y))

        for line in lines:
            if len(line.dots) < Settings.MIN_SEGMENT_LEN:
                lines.remove(line)

        print(" {} lines".format(len(lines)))
        print("Lines:", *lines, sep='\n')

        self.__segments = deepcopy(lines)

    def __len__(self) -> int:
        return len(self.__segments)

    def __getitem__(self, item) -> list:
        return self.__segments[item]

    def plot(self):  # TODO delete after debugging
        for line in self.__segments:
            plt.plot([line.start_x, line.end_x], [abs(line.start_y), abs(line.end_y)])
        plt.grid()
        plt.show()


from typing import List
from collections import deque


class Line:
    def __init__(self, start_x=None, start_y=None, end_x=None, end_y=None, dots=[]):
        self.start_x = start_x
        self.start_y = start_y
        self.end_x = end_x
        self.end_y = end_y
        self.dots = dots

    def __repr__(self):
        return "Line(start_x={}, start_y={}, end_x={}, end_y={}, dots=[{}])".format(
            self.start_x, self.start_y, self.end_x, self.end_y, len(self.dots)
        )

    @property
    def coords(self):
        return self.start_x, self.start_y, self.end_x, self.end_y

    # @property
    # def x1(self):
    #     return self.start_x

    # @property
    # def y1(self):
    #     return self.start_y

    # @property
    # def x2(self):
    #     return self.end_x

    # @property
    # def y2(self):
    #     return self.end_y

    @property
    def center_x(self):
        return (self.start_x + self.end_x) // 2

    @property
    def center_y(self):
        return (self.start_y + self.end_y) // 2

    @property
    def sizeX(self):
        return abs(self.start_x - self.end_x)

    @property
    def sizeY(self):
        return abs(self.start_y - self.end_y)

    def isTiltedCorrectly(self):
        return self.start_y <= self.end_y

    @property
    def k(self):
        return (self.end_y - self.start_y) / (self.end_x - self.start_x)

    @property
    def b(self):
        return self.end_y - self.end_x * self.k

    def copyCoords(self):
        return Line(self.start_x, self.start_y, self.end_x, self.end_y, dots=[])

    def shift(self, dx=0, dy=0):
        self.start_x += dx
        self.start_y += dy
        self.end_x += dx
        self.end_y += dy
        for i in range(len(self.dots)):
            self.dots[i][0] += dx
            self.dots[i][1] += dy
        return self

    def rotateY(self, rotation_center, line=True, dots=False):
        if line:
            self.start_y -= (self.start_y - rotation_center) * 2
            self.end_y -= (self.end_y - rotation_center) * 2

        if dots:
            for i in range(len(self.dots)):
                self.dots[i][1] -= (self.dots[i][1] - rotation_center) * 2
        return self


def shiftLines(lines, count) -> List[Line]:
    result = deque(lines)
    for _ in range(count):
        result.append(result.popleft())
    return list(result)


settings = {
    "grid_size": int(1e5),
    "min_block_size": int(1e3),
    "dot_skip_rate": 10,
    "dotsize": 0.1,
    "fontsize": 8,
    "figsize": (10, 7),

    "min_event_size": int(1e4),
    "lines_join_size": int(1e4) + 3,
    "line_min_size": int(1e4)
}


def linearApproxDots(dots):
    n, sumx, sumy, sumx2, sumxy = len(dots), 0, 0, 0, 0
    for x, y in dots:
        sumx += x
        sumy += y
        sumx2 += x ** 2
        sumxy += x * y

    k = (n * sumxy - (sumx * sumy)) / (n * sumx2 - sumx * sumx)
    b = (sumy - k * sumx) / n
    return k, b


def distance2(x1, y1, x2, y2):
    return (x1 - x2) ** 2 + (y1 - y2) ** 2
