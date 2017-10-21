import copy
import operator


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
            return Matrix2D(self._arr[idx])

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
                return Matrix2D(m)

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

            m = Matrix2D.zeros(self.rows, other.cols)
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
        return Matrix2D(m)

    def transpose(self):
        m = Matrix2D.zeros(self.cols, self.rows)
        for x in range(m.rows):
            for y in range(m.cols):
                m[x, y] = self[y, x]
        return m

    def clone(self):
        return Matrix2D(copy.deepcopy(self._arr), self._rows, self._cols)

    @staticmethod
    def zeros(rows, cols, dtype=float):
        arr = [[dtype(0.) for _ in range(cols)] for _ in range(rows)]
        return Matrix2D(arr, rows, cols)

    @staticmethod
    def ones(rows, cols, dtype=float):
        arr = [[dtype(1.) for _ in range(cols)] for _ in range(rows)]
        return Matrix2D(arr, rows, cols)

    @staticmethod
    def eye(n, dtype=float):
        arr = [
            [dtype(1.) if i == j else dtype(0.) for i in range(n)]
            for j in range(n)
        ]
        return Matrix2D(arr, n, n)
