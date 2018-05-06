from PIL import Image
import numpy as np


def save_image(data, filename):
    data = np.matrix(data)
    rescaled = (255.0 / data.max() * (data - data.min())).astype(np.uint8)
    im = Image.fromarray(rescaled)
    im.save(filename + '.png')
    return None
