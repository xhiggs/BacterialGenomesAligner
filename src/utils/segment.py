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
