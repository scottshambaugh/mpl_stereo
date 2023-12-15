import matplotlib.pyplot as plt
import numpy as np
import inspect

from mpl_stereo import AxesStereo2D, AxesStereo3D


def _testdata():
    data = dict()
    # Parametric equations for a (3,2) trefoil knot
    t = np.linspace(0, 2*np.pi, 100)
    x = np.cos(2*t) * (3 + np.cos(3*t))
    y = np.sin(2*t) * (3 + np.cos(3*t))
    z = np.sin(3*t)
    data['trefoil'] = (x, y, z)
    return data


def test_AxesStereo2D():
    # Smoke test
    x = np.linspace(0, 2*np.pi, 100)
    y = np.sin(x)
    z = np.cos(x)
    axstereo = AxesStereo2D()
    for method in axstereo.known_methods:
        getattr(axstereo, method)(x, y, z)
        getattr(axstereo, method)(x=x, y=y, z=z)
    assert True


def test_AxesStereo3D():
    # Smoke test
    x = np.linspace(0, 2*np.pi, 100)
    y = np.sin(x)
    z = np.cos(x)
    axstereo = AxesStereo3D()
    for method in axstereo.known_methods:
        parameters = inspect.signature(getattr(axstereo.ax_left, method, None)).parameters
        getattr(axstereo, method)(x, y, z)
        if all([keyword in parameters for keyword in ['x', 'y', 'z']]):
            getattr(axstereo, method)(x=x, y=y, z=z)
        elif all([keyword in parameters for keyword in ['xs', 'ys', 'zs']]):
            getattr(axstereo, method)(xs=x, ys=y, zs=z)
        elif all([keyword in parameters for keyword in ['X', 'Y', 'Z']]):
            getattr(axstereo, method)(X=x, Y=y, Z=z)
    assert True


def plotting_tests_2d_pairwise():
    # test plot and scatter
    x, y, z = _testdata()['trefoil']
    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)

    # test bar and stem
    x = y = z = np.arange(10)
    axstereo = AxesStereo2D()
    axstereo.bar(x, y, z, width=1, edgecolor="white", linewidth=0.7)
    axstereo.stem(x, y, z, 'k')


def plotting_tests_3d_pairwise():

    x, y, z = _testdata()['trefoil']
    axstereo = AxesStereo3D()
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)


if __name__ == '__main__':
    plotting_tests_2d_pairwise()
    plotting_tests_3d_pairwise()
    plt.show()
