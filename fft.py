def hadamard(x, y):
    return [a * b for a, b in zip(x, y)]


def vec_add(x, y):
    return [a + b for a, b in zip(x, y)]


def fft(x):
    assert len(x) & (len(x) - 1) == 0, "length of x must be power of 2"
    EXP = 2.7182818284590452353602874713526624977
    PI = 3.1415926535897932384626433832795028841
    n = len(x)
    if n == 1:
        return x
    X_even = fft(x[::2])
    X_odd = fft(x[1::2])
    factor = [EXP ** (-2j * PI * k / n) for k in range(n)]
    return vec_add(X_even, hadamard(X_odd, factor[:n // 2])) + \
        vec_add(X_even, hadamard(X_odd, factor[n // 2:]))

