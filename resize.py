from matrix import transpose
from fft import *

def convert_int(image):
    for i, row in enumerate(image):
        for j, col in enumerate(row):
            image[i][j] = round(row[j].real)
    return image

def zero_seperate(image):
    def helper(row):
        new_row = []
        for i in row:
            new_row += [0, i]
        return row
    if isinstance(image[0], list):
        new_image = []
        row = [0] * (2 * len(image[0]))
        for i in range(len(image)):
            new_row = []
            for j in image[i]:
                new_row.append(j)
                new_row.append(0)
            new_image += [new_row]
            new_image += [row]
        return new_image
        
    else:
        return helper(image)

def resize_image(image, factor):
    image = zero_seperate(image)
    new_image = []
    for row in image:
        vec = [1] * factor + [0] * (len(row) - factor)
        new_row = ifft(hadamard(fft(row), fft(vec)))
        new_image.append(new_row)
    image = transpose(new_image)
    new_image = []
    for col in image:
        vec = [1] * factor + [0] * (len(col) - factor)
        new_col = ifft(hadamard(fft(col), fft(vec)))
        new_image.append(new_col)
    image = transpose(new_image)
    return convert_int(image)
