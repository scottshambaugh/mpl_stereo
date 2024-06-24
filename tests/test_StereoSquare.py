import matplotlib.pyplot as plt

from pathlib import Path
from mpl_toolkits.mplot3d.axes3d import get_test_data
from mpl_stereo import StereoSquare2D, StereoSquare3D
from mpl_stereo.example_data import trefoil, sun_left_right, church_left_right, church_left_cropped


def _testdata():
    data = dict()
    # Parametric equations for a (3,2) trefoil knot
    x, y, z = trefoil()
    data['trefoil'] = (x, y, z)

    # 3D test data
    X, Y, Z = get_test_data(0.05)
    data['3d_default'] = (X, Y, Z)

    # Sun data
    sun_left_data, sun_right_data = sun_left_right()
    data['sun'] = (sun_left_data, sun_right_data)

    # Church data
    church_left_data, church_right_data = church_left_right()
    data['church'] = (church_left_data, church_right_data)

    # Church data cropped
    church_left_data = church_left_cropped()
    data['church_cropped'] = (church_left_data, church_right_data)
    return data


# Smoke test StereoSquare2D
def test_StereoSquare2D():
    wiggle_filepath = Path('test.gif')
    x, y, z = _testdata()['trefoil']
    stereosquare = StereoSquare2D()
    stereosquare.plot(x, y, z)
    stereosquare.wiggle(wiggle_filepath)
    assert wiggle_filepath.exists()
    wiggle_filepath.unlink()  # remove the file
    assert True

    church_left_data, church_right_data = _testdata()['church_cropped']
    stereosquare = StereoSquare2D()
    stereosquare.imshow_stereo(church_left_data, church_right_data, crop=True)
    stereosquare.wiggle(wiggle_filepath)
    assert wiggle_filepath.exists()
    wiggle_filepath.unlink()  # remove the file
    assert True


# Smoke test StereoSquare3D
def test_StereoSquare3D():
    wiggle_filepath = Path('test.gif')
    x, y, z = _testdata()['trefoil']
    stereosquare = StereoSquare3D()
    stereosquare.plot(x, y, z)
    stereosquare.wiggle(wiggle_filepath)
    assert wiggle_filepath.exists()
    wiggle_filepath.unlink()  # remove the file
    assert True


## The following tests are for visual inspection only
def plotting_tests_stereo_square2d(filepath):
    x, y, z = _testdata()['trefoil']
    stereosquare = StereoSquare2D()
    stereosquare.plot(x, y, z, c='k', alpha=0.5)
    stereosquare.wiggle(filepath)


def plotting_tests_stereo_square3d(filepath):
    x, y, z = _testdata()['trefoil']
    stereosquare = StereoSquare3D()
    stereosquare.plot(x, y, z, c='k', alpha=0.5)
    stereosquare.wiggle(filepath)


if __name__ == '__main__':
    filepath = Path('test.gif')

    plotting_tests_stereo_square2d(filepath)
    plotting_tests_stereo_square3d(filepath)

    plt.show()

    filepath.unlink()  # remove the file
