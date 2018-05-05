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


def ifft(x):
    return [a / len(x) for a in fft(x)]


def fftshift(x):
    return x[(len(x) + 1) // 2:] + x[:(len(x) + 1) // 2]


def ifftshift(x):
    return x[len(x) // 2:] + x[:len(x) // 2]


def convolve(x, y):
    length = len(x) + len(y) - 1
    x = x + [0] * (length - len(x))
    next_prime = len(x)
    while next_prime & (next_prime - 1):
        next_prime += 1
    x = x + [0] * (next_prime - len(x))
    y = y + [0] * (len(x) - len(y))
    return [ifft(hadamard(fft(x), fft(y)))[-i]
            for i in range(length)]
