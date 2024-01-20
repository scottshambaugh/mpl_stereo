import numpy as np
import matplotlib.image as mpimg
from pathlib import Path

MPL_STEREO_DATA_DIR = Path(__file__).parent.resolve() / 'data'

def trefoil(i: int = 0, n_points: int = 100, n_steps: int = 10):
    dt = 2*np.pi*i/n_points/n_steps
    t = np.linspace(0, 2*np.pi, n_points) + dt
    x = np.cos(2*t) * (3 + np.cos(3*t))
    y = np.sin(2*t) * (3 + np.cos(3*t))
    z = np.sin(3*t)
    return x, y, z

def sun_left_right():
    # Get only the red channel for a grayscale image
    sun_left_data = mpimg.imread(MPL_STEREO_DATA_DIR / 'sun_left.png')[:, :, 0]
    sun_right_data = mpimg.imread(MPL_STEREO_DATA_DIR / 'sun_right.png')[:, :, 0]
    return sun_left_data, sun_right_data

def church_left_right():
    church_left_data = mpimg.imread(MPL_STEREO_DATA_DIR / 'church_left.jpg')
    church_right_data = mpimg.imread(MPL_STEREO_DATA_DIR / 'church_right.jpg')
    return church_left_data, church_right_data

def church_left_cropped():
    church_left_data = mpimg.imread(MPL_STEREO_DATA_DIR / 'church_left_cropped.jpg')
    return church_left_data
