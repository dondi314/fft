def hadamard(x, y):
    return [a * b for a, b in zip(x, y)]


def vec_add(x, y):
    return [a + b for a, b in zip(x, y)]


def fft(x, inverse=False):
    coef = 2j if inverse else -2j
    assert len(x) & (len(x) - 1) == 0, "length of x must be power of 2"
    EXP = 2.7182818284590452353602874713526624977
    PI = 3.1415926535897932384626433832795028841
    n = len(x)
    if n == 1:
        return x
    X_even = fft(x[::2])
    X_odd = fft(x[1::2])
    factor = [EXP ** (coef * PI * k / n) for k in range(n)]
    return vec_add(X_even, hadamard(X_odd, factor[:n // 2])) + \
        vec_add(X_even, hadamard(X_odd, factor[n // 2:]))


def ifft(x):
    return [a / len(x) for a in fft(x, True)]


def fftshift(x):
    return x[(len(x) + 1) // 2:] + x[:(len(x) + 1) // 2]


def ifftshift(x):
    return x[len(x) // 2:] + x[:len(x) // 2]


def transform(vector, inverse=False):
    n = len(vector)
    if (n & (n - 1)) == 0:
        return fft(vector, inverse)
    else:
        return bluestein(vector, inverse)

def convolve(x, y, realoutput=True):
    pad = len(x)
    while pad & (pad - 1):
        pad += 1
    x += [0] * (pad - len(x))
    pad = len(y)
    while pad & (pad - 1):
        pad += 1
    y += [0] * (pad - len(y))
    x = transform(x)
    y = transform(y)
    z = [a * b for a, b in zip(x, y)]
    result = [transform(z, inverse=True)[-i] 
              for i in range(len(z))]
    if realoutput:
        return [x.real for x in result]
    else:
        return result


def bluestein(vector, inverse=False):
    EXP = 2.7182818284590452353602874713526624977
    PI = 3.1415926535897932384626433832795028841
    n = len(vector)
    m = 2**((n * 2).bit_length())
    coef = (1j if inverse else -1j) * (PI / n)
    exptable = [EXP ** (((i * i) % (n * 2)) * coef) for i in range(n)]
    a = [(x * y) for x, y in zip(vector, exptable)] + [0] * (m - n)
    b = exptable[ : n] + [0] * (m - (n * 2 - 1)) + exptable[ : 0: -1]
    b = [x.conjugate() for x in b]
    c = convolve(a, b, False)[ : n]
    return [(x * y) for x, y in zip(c, exptable)]
