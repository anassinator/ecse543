import copy
import operator
from collections import defaultdict


class Matrix2D(object):

    def __init__(self, arr, rows=None, cols=None, dtype=None):
        if not isinstance(arr, list):
            raise ValueError("incompatible type: {}".format(type(arr)))

        if len(arr) == 0:
            raise ValueError("incompatible size: 0")

        if not isinstance(arr[0], list):
            arr = [arr]

        self._arr = arr

        if rows is not None:
            self._rows = rows
        else:
            self._rows = len(arr)

        if cols is not None:
            self._cols = cols
        else:
            self._cols = len(arr[0])

        if dtype is not None:
            self._dtype = dtype
            if type(arr[0][0]) != dtype:
                self._arr = self.map(dtype)._arr
        else:
            self._dtype = type(arr[0][0])

    @property
    def rows(self):
        return self._rows

    @property
    def cols(self):
        return self._cols

    @property
    def shape(self):
        return self._rows, self._cols

    @property
    def dtype(self):
        return self._dtype

    @property
    def T(self):
        return self.transpose()

    def issquare(self):
        return self._rows == self._cols

    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self._arr[idx]

        if isinstance(idx, slice):
            return self.__class__(self._arr[idx])

        if isinstance(idx, tuple):
            x, y = idx
            if isinstance(x, int):
                return self._arr[x][y]
            else:
                m = []
                rows = self._arr[x]
                for a in rows:
                    if isinstance(y, int):
                        m.append([a[y]])
                    else:
                        m.append(a[y])
                return self.__class__(m)

        raise NotImplementedError()

    def __setitem__(self, idx, value):
        if isinstance(idx, int):
            self._arr[idx] = list(map(lambda _: value, self._arr[idx]))
            return

        if isinstance(idx, slice):
            self._arr[idx] = list(
                map(lambda row: list(map(lambda _: value, row)),
                    self._arr[idx]))
            return

        if isinstance(idx, tuple):
            x, y = idx
            if isinstance(x, int):
                if isinstance(y, int):
                    self._arr[x][y] = value
                    return
                self._arr[x][y] = list(map(lambda _: value, self._arr[x][y]))
                return

            def apply(row):
                row[y] = list(map(lambda _: value, row[y]))
                return row

            self._arr[x] = list(map(apply, self._arr[x]))

    def __lt__(self, other):
        return self._arr < other._arr

    def __le__(self, other):
        return self._arr <= other._arr

    def __eq__(self, other):
        return self._arr == other._arr

    def __ne__(self, other):
        return self._arr != other._arr

    def __gt__(self, other):
        return self._arr > other._arr

    def __ge__(self, other):
        return self._arr >= other._arr

    def __add__(self, other):
        return self.__operate(other, operator.add)

    def __sub__(self, other):
        return self.__operate(other, operator.sub)

    def __mul__(self, other):
        if isinstance(other, Matrix2D):
            return self.dot(other)
        return self.__operate(other, operator.mul)

    def __truediv__(self, other):
        return self.__operate(other, operator.truediv)

    def __mod__(self, other):
        return self.__operate(other, operator.mod)

    def __neg__(self):
        return self.map(lambda x: -x)

    def __len__(self):
        return len(self._arr)

    def __repr__(self):
        return "Matrix2D({})".format(self._arr.__repr__())

    def __str__(self):
        return self.__repr__()

    def __operate(self, other, op):
        if isinstance(other, Matrix2D):
            if self.shape != other.shape:
                raise ValueError("incompatible shapes: {} != {}".format(
                    self.shape, other.shape))

            m = self.clone()
            for x in range(m.rows):
                for y in range(m.cols):
                    m[x, y] = op(self[x, y], other[x, y])
            return m

        if isinstance(other, list):
            raise NotImplementedError()

        return self.map(lambda x: op(x, other))

    def tosparse(self):
        return SparseMatrix2D(self._arr, self._rows, self._cols, self._dtype)

    def tomatrix(self):
        return self.clone()

    def tolist(self):
        return copy.deepcopy(self._arr)

    def flatten(self):
        arr = [
            self._arr[x][y]
            for y in range(self.cols)
            for x in range(self.rows)
        ]
        return arr

    def dot(self, other):
        if isinstance(other, Matrix2D):
            if self.cols != other.rows:
                raise ValueError("incompatible shape: {} and {}".format(
                    self.shape, other.shape))

            m = self.__class__.zeros(self.rows, other.cols)
            for x in range(m.rows):
                for y in range(m.cols):
                    m[x, y] = sum(
                        map(lambda x: x[0] * x[1],
                            zip(self[x], other[:, y].flatten())))
            return m

        if isinstance(other, list):
            m = Matrix2D([other]).T
            return self.dot(m)

        raise NotImplementedError()

    def map(self, f):
        m = self.clone()
        m = list(map(lambda row: list(map(f, row)), m))
        return self.__class__(m)

    def transpose(self):
        m = self.__class__.zeros(self.cols, self.rows)
        for x in range(m.rows):
            for y in range(m.cols):
                m[x, y] = self[y, x]
        return m

    def clone(self):
        return self.__class__(copy.deepcopy(self._arr), self._rows, self._cols)

    @classmethod
    def zeros(cls, rows, cols, dtype=float):
        arr = [[dtype(0.) for _ in range(cols)] for _ in range(rows)]
        return cls(arr, rows, cols)

    @classmethod
    def ones(cls, rows, cols, dtype=float):
        arr = [[dtype(1.) for _ in range(cols)] for _ in range(rows)]
        return cls(arr, rows, cols)

    @classmethod
    def eye(cls, n, dtype=float):
        arr = [
            [dtype(1.) if i == j else dtype(0.) for i in range(n)]
            for j in range(n)
        ]
        return cls(arr, n, n)


class SparseMatrix2D(Matrix2D):

    def __init__(self, arr, rows=None, cols=None, dtype=None):
        if not isinstance(arr, list):
            raise ValueError("incompatible type: {}".format(type(arr)))

        if len(arr) == 0:
            raise ValueError("incompatible size: 0")

        if not isinstance(arr[0], list):
            arr = [arr]

        self._arr = arr

        if rows is not None:
            self._rows = rows
        else:
            self._rows = len(arr)

        if cols is not None:
            self._cols = cols
        else:
            self._cols = len(arr[0])

        if dtype is not None:
            self._dtype = dtype
            if type(arr[0][0]) != dtype:
                self._arr = self.map(dtype)._arr
        else:
            self._dtype = type(arr[0][0])

        self._elements = defaultdict(lambda: defaultdict(self._dtype))
        for x in range(self._rows):
            for y in range(self._cols):
                if arr[x][y]:
                    self._elements[x][y] = arr[x][y]

    def __setitem__(self, idx, value):
        super().__setitem__(idx, value)

        if isinstance(idx, int):
            if value == 0:
                if idx in self._elements:
                    del self._elements[idx]
            else:
                self._elements[idx] = {
                    i: self._arr[idx][i]
                    for i in self._cols
                    if self._arr[idx][i] != 0
                }

        if isinstance(idx, slice):
            for x in _slice_to_range(idx):
                self._elements[x] = {
                    y: self._arr[x][y]
                    for y in self._cols
                    if self._arr[x][y] != 0
                }

        if isinstance(idx, tuple):
            x, y = idx
            if isinstance(x, int):
                if isinstance(y, int):
                    self._elements[x][y] = value
                    if value == 0:
                        del self._elements[x][y]
                    return

                self._elements[x] = {
                    y: self._arr[x][y]
                    for y in self._cols
                    if self._arr[x][y] != 0
                }
                return

            for x in _slice_to_range(idx):
                self._elements[x] = {
                    y: self._arr[x][y]
                    for y in self._cols
                    if self._arr[x][y] != 0
                }

    def __repr__(self):
        return "SparseMatrix2D({})".format(self._arr.__repr__())

    def __operate(self, other, op):
        if isinstance(other, SparseMatrix2D):
            seen = set()
            for x in self._elements:
                for y in self._elements[x]:
                    self._arr[x][y] = op(self._arr[x][y], other[x][y])
                    self._elements[x][y] = self._arr[x][y]
                    seen.add((x, y))

            for x in other._elements:
                for y in other._elements[x]:
                    if (x, y) in seen:
                        continue
                    self._arr[x][y] = op(self._arr[x][y], other[x][y])
                    self._elements[x][y] = self._arr[x][y]

        return super().__operate(other, op)

    def dot(self, other):
        if isinstance(other, Matrix2D):
            if self.cols != other.rows:
                raise ValueError("incompatible shape: {} and {}".format(
                    self.shape, other.shape))

            if False and isinstance(other, SparseMatrix2D):
                m = SparseMatrix2D.zeros(self.rows, other.cols)
                for x in other._elements:
                    row = self._elements[x]
                    for y in other._elements[x]:
                        for _x in row:
                            m[x, y] += row[_x] * other._arr[_x][y]
            else:
                m = Matrix2D.zeros(self.rows, other.cols)
                for x in range(m.rows):
                    row = self._elements[x]
                    for y in range(m.cols):
                        for _x in row:
                            m[x, y] += row[_x] * other[_x, y]
            return m

        if isinstance(other, list):
            m = Matrix2D([other]).T
            return self.dot(m)

        raise NotImplementedError()

    def map(self, f):
        zero_value = f(0.)
        if zero_value == 0:
            m = self.clone()
            for x in m._elements:
                for y in m._elements[x]:
                    value = f(m._elements[x][y])
                    m._elements[x][y] = value
                    m._arr[x][y] = value
                    if value == 0:
                        del m._elements[x][y]
            return m
        else:
            return super().map(f)

    def tosparse(self):
        return self.clone()

    def tomatrix(self):
        return Matrix2D(self._arr, self._rows, self._cols, self._dtype)


def _slice_to_range(s):
    start = s.start if s.start is not None else 0
    stop = s.stop if s.stop is not None else -1
    step = s.step if s.start is not None else 1
    return range(start, stop, step)
