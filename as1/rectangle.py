class Rectangle(object):

    def __init__(self, width, height, bottom_left_x, bottom_left_y):
        self._width = float(width)
        self._height = float(height)

        self._center_x = float(bottom_left_x) + width / 2.0
        self._center_y = float(bottom_left_y) + height / 2.0

        self._top = bottom_left_y + self._height
        self._bottom = bottom_left_y
        self._right = bottom_left_x + self._width
        self._left = bottom_left_x

    def __repr__(self):
        return "Rectangle{s.shape} @ {s.center}".format(s=self)

    def __str__(self):
        return self.__repr__()

    def __contains__(self, coord):
        if isinstance(coord, tuple):
            x, y = coord
            return (self._left <= x <= self._right and
                    self._bottom <= y <= self._top)
        raise NotImplementedError()

    def on_edge(self, coord):
        x, y = coord
        return (x in (self._left, self._right) or
                y in (self._bottom, self._top))

    @property
    def width(self):
        return self._width

    @property
    def height(self):
        return self._height

    @property
    def shape(self):
        return self._width, self._height

    @property
    def center(self):
        return self._center_x, self._center_y

    @property
    def top(self):
        return self._top

    @property
    def bottom(self):
        return self._bottom

    @property
    def right(self):
        return self._right

    @property
    def left(self):
        return self._left

    @property
    def top_right(self):
        return self._top, self._right

    @property
    def top_left(self):
        return self._top, self._left

    @property
    def bottom_right(self):
        return self._bottom, self._right

    @property
    def bottom_left(self):
        return self._bottom, self._left
