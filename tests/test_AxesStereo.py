import matplotlib.pyplot as plt
import numpy as np
import inspect
from mpl_toolkits.mplot3d.axes3d import get_test_data

from mpl_stereo import AxesStereo2D, AxesStereo3D, AxesAnaglyph
from mpl_stereo.example_data import trefoil, sun_left_right, church_left_right


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
    return data


def test_existing_fig():
    # Smoke tests
    fig = plt.figure()
    AxesStereo2D(fig=fig)
    fig = plt.figure()
    AxesStereo3D(fig=fig)
    fig = plt.figure()
    AxesAnaglyph(fig=fig)
    assert True


def test_existing_ax():
    # Smoke tests
    _, axs = plt.subplots(2, 1)
    AxesStereo2D(axs=axs)
    _, axs = plt.subplots(2, 1, subplot_kw=dict(projection='3d'))
    AxesStereo3D(axs=axs)
    _, ax = plt.subplots()
    AxesAnaglyph(ax=ax)
    assert True


def test_AxesStereo2D():
    # Smoke test plotting
    x, y, z = _testdata()['trefoil']
    axstereo = AxesStereo2D()
    for method in axstereo.known_methods:
        if method == 'text':
            getattr(axstereo, method)(0, 0, 0, 'text')
            getattr(axstereo, method)(x=0, y=0, z=0, s='text')
        else:
            getattr(axstereo, method)(x, y, z)
            getattr(axstereo, method)(x=x, y=y, z=z)
    assert True

    # Smoke test eye balance
    axstereo = AxesStereo2D(eye_balance=1.0)
    axstereo.plot(x, y, z)
    assert True

    # Smoke test zzero
    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z, zzero=min(z))
    assert True

    # Smoke test scatter sorting
    axstereo = AxesStereo2D()
    axstereo.scatter(x, y, z, c=z, cmap='viridis')
    assert True

    # Smoke test non-plotting methods
    axstereo.set_title('title')
    assert True

def test_AxesStereo2D_zlim():
    x = y = z = np.arange(10)
    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z/2 + 1)
    assert axstereo.zlim == (1, 5.5)

    axstereo.plot(x, y, z)
    assert axstereo.zlim == (0, 9)

    axstereo.plot(x, y, z/2)
    assert axstereo.zlim == (0, 9)

    axstereo.set_zlim((0, 5), zautoscale=True)
    assert axstereo.zlim == (0, 5)

    axstereo.redraw()  # we set autoscale=True, so redraw() should reset zlim
    assert axstereo.zlim == (0, 9)

    axstereo.set_zlim((0, 5), zautoscale=False)
    axstereo.plot(x, y, z)
    assert axstereo.zlim == (0, 5)

    axstereo.plot(x, y, z, zlim=(1, 2))
    assert axstereo.zlim == (1, 2)

    assert axstereo.get_zlim() == (1, 2)

    n_artists = 5  # from above plot()'s
    assert len(axstereo.artists_left) == n_artists
    assert len(axstereo.artists_right) == n_artists
    assert len(axstereo.artist_args) == n_artists

def test_AxesAnaglyph_zlim():
    x = y = z = np.arange(10)
    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z/2 + 1)
    assert axstereo.zlim == (1, 5.5)

    axstereo.plot(x, y, z)
    assert axstereo.zlim == (0, 9)

    axstereo.plot(x, y, z/2)
    assert axstereo.zlim == (0, 9)

    axstereo.set_zlim((0, 5), zautoscale=True)
    assert axstereo.zlim == (0, 5)

    axstereo.redraw()  # we set autoscale=True, so redraw() should reset zlim
    assert axstereo.zlim == (0, 9)

    axstereo.set_zlim((0, 5), zautoscale=False)
    axstereo.plot(x, y, z)
    assert axstereo.zlim == (0, 5)

    axstereo.plot(x, y, z, zlim=(1, 2))
    assert axstereo.zlim == (1, 2)

    assert axstereo.get_zlim() == (1, 2)

    n_artists = 5  # from above plot()'s
    assert len(axstereo.artists_left) == n_artists
    assert len(axstereo.artists_right) == n_artists
    assert len(axstereo.artist_args) == n_artists

def test_AxesStereo3D():
    # Smoke test plotting
    x, y, z = _testdata()['trefoil']
    X, Y, Z = _testdata()['3d_default']
    axstereo = AxesStereo3D()
    for method in axstereo.known_methods:
        parameters = inspect.signature(getattr(axstereo.ax_left, method, None)).parameters
        if all([keyword in parameters for keyword in ['x', 'y', 'z']]):
            getattr(axstereo, method)(x, y, z)
            getattr(axstereo, method)(x=x, y=y, z=z)
        elif all([keyword in parameters for keyword in ['xs', 'ys', 'zs']]):
            getattr(axstereo, method)(x, y, z)
            getattr(axstereo, method)(xs=x, ys=y, zs=z)
        elif all([keyword in parameters for keyword in ['X', 'Y', 'Z']]):
            getattr(axstereo, method)(X, Y, Z)
            getattr(axstereo, method)(X=X, Y=Y, Z=Z)
    assert True

    # Smoke test eye balance
    axstereo = AxesStereo3D(eye_balance=1.0)
    axstereo.plot(x, y, z)
    assert True

    # Smoke test non-plotting methods
    axstereo.set_title('title')
    assert True


def test_AxesAnaglyph():
    # Smoke test plotting
    x, y, z = _testdata()['trefoil']
    axstereo = AxesAnaglyph()
    for method in axstereo.known_methods:
        if method == 'text':
            getattr(axstereo, method)(0, 0, 0, 'text')
            getattr(axstereo, method)(x=0, y=0, z=0, s='text')
        else:
            getattr(axstereo, method)(x, y, z)
            getattr(axstereo, method)(x=x, y=y, z=z)
    assert True

    # Smoke test eye balance
    axstereo = AxesAnaglyph(eye_balance=1.0)
    axstereo.plot(x, y, z)
    assert True

    # Smoke test zzero
    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z, zzero=min(z))
    axstereo = AxesAnaglyph(zzero=min(z))
    axstereo.plot(x, y, z)
    assert True

    # Smoke test zscale
    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z, zscale=np.ptp(z))
    axstereo = AxesAnaglyph(zscale=np.ptp(z))
    axstereo.plot(x, y, z)
    assert True

    # Smoke test non-plotting methods
    axstereo.set_title('title')
    assert True


def test_AxesAnaglyph_imshow_steroe():
    # Smoke test imshow_stereo
    axstereo = AxesAnaglyph()
    axstereo.imshow_stereo(np.array([[1]]), np.array([[1]]))
    assert True


def plotting_tests_2d_pairwise():
    # test plot and scatter
    x, y, z = _testdata()['trefoil']
    axstereo = AxesStereo2D(eye_balance=-1.0)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.grid(True)
    axstereo.set_title('eye_balance=-1')

    axstereo = AxesStereo2D(eye_balance=1.0)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.grid(True)
    axstereo.set_title('eye_balance=1')

    # test bar and stem
    x = y = z = np.arange(10)
    axstereo = AxesStereo2D()
    axstereo.bar(x, y, z, width=1, edgecolor="white", linewidth=0.7)
    axstereo.stem(x, y, z, 'k')


def plotting_tests_2d_pairwise_zlim():
    # test plot and scatter
    x, y, z = _testdata()['trefoil']
    axstereo = AxesStereo2D(eye_balance=-1.0)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.grid(True)
    axstereo.set_zlim((-0.5, 0.5))
    axstereo.set_title('zlim=(-0.5, 0.5)')

    axstereo = AxesStereo2D(eye_balance=1.0)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.grid(True)
    axstereo.set_zlim((-2, 2))
    axstereo.set_title('zlim=(-2, 2)')


def plotting_tests_anaglyph_pairwise_zzero():
    # test plot and scatter
    x, y, z = _testdata()['trefoil']
    axstereo = AxesAnaglyph(zzero=min(z))
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.grid(True)
    axstereo.set_title('zzero=min(zzero)')

    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z, c='k', alpha=0.2, zzero=None)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10, zzero=None)
    axstereo.grid(True)
    axstereo.set_title('zzero=None')

    axstereo = AxesAnaglyph(zscale=np.ptp(z))
    axstereo.plot(x, y, z, c='k', alpha=0.2, zzero=max(z))
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10, zzero=max(z))
    axstereo.grid(True)
    axstereo.set_title('zzero=max(zzero)')


def plotting_tests_3d():
    x, y, z = _testdata()['trefoil']
    axstereo = AxesStereo3D()
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)


def plotting_tests_anaglyph_pairwise():
    x, y, z = _testdata()['trefoil']
    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)

    axstereo = AxesAnaglyph(eye_balance=1.0)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)


def plotting_tests_anaglyph_imshow_stereo():
    sun_left_data, sun_right_data = _testdata()['sun']
    axstereo = AxesAnaglyph()
    axstereo.imshow_stereo(sun_left_data, sun_right_data)

    church_left_data, church_right_data = _testdata()['church']
    axstereo = AxesAnaglyph()
    axstereo.imshow_stereo(church_left_data, church_right_data)


if __name__ == '__main__':
    plotting_tests_2d_pairwise()
    plotting_tests_2d_pairwise_zlim()
    plotting_tests_3d()
    plotting_tests_anaglyph_pairwise()
    plotting_tests_anaglyph_pairwise_zzero()
    plotting_tests_anaglyph_imshow_stereo()
    plt.show()
