from matrix import Matrix2D

data = {
    0.0: 0.0,
    0.2: 14.7,
    0.4: 36.5,
    0.6: 71.7,
    0.8: 121.4,
    1.0: 197.4,
    1.1: 256.2,
    1.2: 348.7,
    1.3: 540.6,
    1.4: 1062.8,
    1.5: 2318.0,
    1.6: 4781.9,
    1.7: 8687.4,
    1.8: 13924.3,
    1.9: 22650.2,
}

ALL_B = sorted(data.keys())


def get_data(B=ALL_B):
    X = Matrix2D(B)
    Y = Matrix2D.zeros(X.rows, X.cols)
    for i, b in enumerate(B):
        Y[0, i] = data[b]
    return X.T, Y.T
