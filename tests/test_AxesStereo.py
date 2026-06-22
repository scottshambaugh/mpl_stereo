import matplotlib.pyplot as plt
import numpy as np
import inspect
import pytest

from pathlib import Path
from mpl_toolkits.mplot3d.axes3d import get_test_data
from mpl_stereo import AxesStereo2D, AxesStereo3D, AxesAnaglyph, AxesAnaglyph3D
from mpl_stereo.example_data import trefoil, sun_left_right, church_left_right, church_left_cropped


def _testdata():
    data = dict()
    # Parametric equations for a (3,2) trefoil knot
    x, y, z = trefoil()
    data["trefoil"] = (x, y, z)

    # 3D test data
    X, Y, Z = get_test_data(0.05)
    data["3d_default"] = (X, Y, Z)

    # Sun data
    sun_left_data, sun_right_data = sun_left_right()
    data["sun"] = (sun_left_data, sun_right_data)

    # Church data
    church_left_data, church_right_data = church_left_right()
    data["church"] = (church_left_data, church_right_data)

    # Church data cropped
    church_left_data = church_left_cropped()
    data["church_cropped"] = (church_left_data, church_right_data)
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
    _, axs = plt.subplots(2, 1, subplot_kw=dict(projection="3d"))
    AxesStereo3D(axs=axs)
    _, ax = plt.subplots()
    AxesAnaglyph(ax=ax)
    assert True


def test_AxesStereo2D():
    # Smoke test plotting
    x, y, z = _testdata()["trefoil"]
    axstereo = AxesStereo2D()
    for method in axstereo.known_methods:
        if method == "text":
            getattr(axstereo, method)(0, 0, 0, "text")
            getattr(axstereo, method)(x=0, y=0, z=0, s="text")
        elif method == "fill_between":
            # fill_between takes (x, y1, y2) with the depth as a `z` keyword
            getattr(axstereo, method)(x, y, z=z)
            getattr(axstereo, method)(x, y, y * 0.5, z=z)
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
    axstereo.scatter(x, y, z, c=z, cmap="viridis")
    assert True

    # Smoke test non-plotting methods
    axstereo.set_title("title")
    assert True


def test_AxesStereo2D_zlim():
    x = y = z = np.arange(10)
    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z / 2 + 1)
    assert axstereo.zlim == (1, 5.5)

    axstereo.plot(x, y, z)
    assert axstereo.zlim == (0, 9)

    axstereo.plot(x, y, z / 2)
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


def test_fill_between():
    x = np.linspace(0, 1, 20)
    y1 = np.zeros_like(x)
    y2 = np.ones_like(x)
    z = x  # depth increases along x

    # The depth (passed as a `z` keyword) offsets the band between the two eyes,
    # while y1 and y2 are preserved as the fill boundaries.
    axstereo = AxesStereo2D(ipd=65)
    res_left, res_right = axstereo.fill_between(x, y1, y2, z=z)
    assert len(axstereo.artists_left) == 1
    assert len(axstereo.artists_right) == 1
    vleft = res_left.get_paths()[0].vertices
    vright = res_right.get_paths()[0].vertices
    assert not np.allclose(vleft[:, 0], vright[:, 0])  # horizontally shifted
    assert axstereo.zlim == (0, 1)  # z drove the z limits

    # Without a z keyword there is no stereo offset, both eyes are identical,
    # and a warning is raised.
    axstereo = AxesStereo2D(ipd=65)
    with pytest.warns(UserWarning, match="z` keyword"):
        res_left, res_right = axstereo.fill_between(x, y1, y2)
    assert np.allclose(res_left.get_paths()[0].vertices, res_right.get_paths()[0].vertices)

    # A redraw (triggered by a later plot) reproduces the band correctly.
    axstereo = AxesStereo2D(ipd=65)
    axstereo.fill_between(x, y1, y2, z=z)
    axstereo.plot(x, y2, z)  # forces a redraw of earlier artists
    assert len(axstereo.artists_left) == 2


def test_transformed_xaxis():
    x = np.geomspace(1, 1000, 25)
    y = np.log10(x)
    z = np.linspace(0, 1, 25)
    zlim, zscale = (0.0, 1.0), 2.0

    # Linear x-axis: the transform is the identity, so the shift is in data space.
    axlin = AxesStereo2D(ipd=65)
    axlin.set_zlim(zlim, zscale, zautoscale=False)
    _, rr_lin = axlin.plot(x, y, z)
    lin_shift = x - rr_lin[0].get_xdata()
    assert np.any(lin_shift != 0)  # there is a parallax offset

    # Log x-axis with the same zscale: the offset is applied in the transformed
    # (log) space, so the screen-uniform parallax matches the linear shift.
    axlog = AxesStereo2D(ipd=65)
    axlog.set_xscale("log")
    axlog.set_zlim(zlim, zscale, zautoscale=False)
    rl_log, rr_log = axlog.plot(x, y, z)
    tf = axlog.ax_left.xaxis.get_transform()
    log_shift = tf.transform(x) - tf.transform(rr_log[0].get_xdata())
    assert np.allclose(log_shift, lin_shift)
    assert np.allclose(rl_log[0].get_xdata(), x)  # left eye unshifted (eye_balance -1)
    # the data-space shift is non-uniform on a log axis, unlike the linear case
    assert not np.allclose(x - rr_log[0].get_xdata(), lin_shift)

    # Setting the scale after plotting redraws and re-places the artists.
    ax = AxesStereo2D(ipd=65)
    ax.set_zlim(zlim, zscale, zautoscale=False)
    ax.plot(x, y, z)
    before = ax.artists_right[0].get_xdata().copy()
    ax.set_xscale("log")
    assert not np.allclose(before, ax.artists_right[0].get_xdata())


def test_semilogx_loglog():
    x = np.geomspace(1, 1000, 25)
    y = np.log10(x)
    z = np.linspace(0, 1, 25)
    zlim, zscale = (0.0, 1.0), 2.0

    # semilogx is equivalent to set_xscale("log") followed by plot
    a = AxesStereo2D(ipd=65)
    a.set_zlim(zlim, zscale, zautoscale=False)
    ra = a.semilogx(x, y, z)
    b = AxesStereo2D(ipd=65)
    b.set_zlim(zlim, zscale, zautoscale=False)
    b.set_xscale("log")
    rb = b.plot(x, y, z)
    assert a.ax_left.get_xscale() == "log"
    assert np.allclose(ra[1][0].get_xdata(), rb[1][0].get_xdata())

    # loglog log-scales both axes on both eyes
    c = AxesStereo2D(ipd=65)
    c.loglog(x, y, z)
    assert c.ax_left.get_xscale() == "log" and c.ax_left.get_yscale() == "log"
    assert c.ax_right.get_xscale() == "log" and c.ax_right.get_yscale() == "log"

    # the base keyword is forwarded through to the scale
    d = AxesStereo2D(ipd=65)
    d.semilogx(x, y, z, base=2)
    assert d.ax_left.xaxis.get_transform().base == 2

    # the anaglyph supports them too
    f = AxesAnaglyph(ipd=80)
    f.loglog(x, y, z)
    assert f.ax.get_xscale() == "log" and f.ax.get_yscale() == "log"


def test_AxesAnaglyph_zlim():
    x = y = z = np.arange(10)
    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z / 2 + 1)
    assert axstereo.zlim == (1, 5.5)

    axstereo.plot(x, y, z)
    assert axstereo.zlim == (0, 9)

    axstereo.plot(x, y, z / 2)
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

    axstereo.autoscale_z()  # will redraw
    assert axstereo.zlim == (0, 9)

    n_artists = 5  # from above plot()'s
    assert len(axstereo.artists_left) == n_artists
    assert len(axstereo.artists_right) == n_artists
    assert len(axstereo.artist_args) == n_artists


def test_AxesStereo3D():
    # Smoke test plotting
    x, y, z = _testdata()["trefoil"]
    X, Y, Z = _testdata()["3d_default"]
    axstereo = AxesStereo3D()
    for method in axstereo.known_methods:
        parameters = inspect.signature(getattr(axstereo.ax_left, method, None)).parameters
        if all([keyword in parameters for keyword in ["x", "y", "z"]]):
            getattr(axstereo, method)(x, y, z)
            getattr(axstereo, method)(x=x, y=y, z=z)
        elif all([keyword in parameters for keyword in ["xs", "ys", "zs"]]):
            getattr(axstereo, method)(x, y, z)
            getattr(axstereo, method)(xs=x, ys=y, zs=z)
        elif all([keyword in parameters for keyword in ["X", "Y", "Z"]]):
            getattr(axstereo, method)(X, Y, Z)
            getattr(axstereo, method)(X=X, Y=Y, Z=Z)
    assert True

    # Smoke test eye balance
    axstereo = AxesStereo3D(eye_balance=1.0)
    axstereo.plot(x, y, z)
    assert True

    # Smoke test non-plotting methods
    axstereo.set_title("title")
    assert True


def test_anaglyph_stem():
    x = np.linspace(0, 6, 20)
    y = np.sin(x) + 1.5
    z = np.cos(x)
    axstereo = AxesAnaglyph()
    assert "stem" in axstereo.known_methods
    rl, rr = axstereo.stem(x, y, z)
    # stem takes no color/alpha keyword, so its container is colored per eye
    assert rl.markerline.get_color() == axstereo.colors[1]
    assert rr.markerline.get_color() == axstereo.colors[0]
    assert rl.markerline.get_alpha() == axstereo.alpha
    assert rr.stemlines.get_alpha() == axstereo.alpha


def test_AxesAnaglyph3D(tmp_path):
    x, y, z = _testdata()["trefoil"]

    axstereo = AxesAnaglyph3D(ipd=65)
    axstereo.plot(x, y, z)

    # The two eye axes are stacked on top of each other (overlapping)
    assert axstereo.ax_left.get_position().bounds == axstereo.ax_right.get_position().bounds

    # Each eye's data is drawn in its glasses-lens color with blending alpha
    line_left = axstereo.ax_left.get_lines()[0]
    line_right = axstereo.ax_right.get_lines()[0]
    assert line_left.get_color() == axstereo.colors[1]
    assert line_right.get_color() == axstereo.colors[0]
    assert line_left.get_alpha() == axstereo.alpha

    # Negative ipd is forced positive (anaglyphs are not cross-view)
    assert AxesAnaglyph3D(ipd=-65).ipd == 65

    # Smoke test saving to file
    filepath = tmp_path / "anaglyph3d.png"
    axstereo.save(filepath)
    assert filepath.exists()


def test_AxesAnaglyph():
    # Smoke test plotting
    x, y, z = _testdata()["trefoil"]
    axstereo = AxesAnaglyph()
    for method in axstereo.known_methods:
        if method == "text":
            getattr(axstereo, method)(0, 0, 0, "text")
            getattr(axstereo, method)(x=0, y=0, z=0, s="text")
        elif method == "fill_between":
            # fill_between takes (x, y1, y2) with the depth as a `z` keyword
            getattr(axstereo, method)(x, y, z=z)
            getattr(axstereo, method)(x, y, y * 0.5, z=z)
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
    axstereo.set_title("title")
    assert True


def test_AxesAnaglyph_imshow_stereo():
    # Smoke test imshow_stereo for scalar data
    sun_left_data, sun_right_data = _testdata()["sun"]
    for cmap in [None, "Oranges_r", "viridis"]:
        axstereo = AxesAnaglyph()
        axstereo.imshow_stereo([sun_left_data, sun_right_data], cmap=cmap)
    assert True

    # Smoke test imshow_stereo for RGB data
    church_left_data, church_right_data = _testdata()["church"]
    for method in ["dubois", "photoshop", "photoshop2"]:
        axstereo = AxesAnaglyph()
        axstereo.imshow_stereo([church_left_data, church_right_data], method=method)
    assert True

    # Smoke test cropping
    church_left_data, church_right_data = _testdata()["church_cropped"]
    axstereo = AxesAnaglyph()
    axstereo.imshow_stereo([church_left_data, church_right_data], crop=True)
    assert True


def test_wiggle():
    # Smoke test wiggle
    wiggle_filepath = Path("test.gif")
    x, y, z = _testdata()["trefoil"]

    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z)
    axstereo.wiggle(wiggle_filepath)
    assert wiggle_filepath.exists()
    wiggle_filepath.unlink()  # remove the file
    assert True

    axstereo = AxesStereo3D()
    axstereo.plot(x, y, z)
    axstereo.wiggle(wiggle_filepath)
    assert wiggle_filepath.exists()
    wiggle_filepath.unlink()  # remove the file
    assert True

    church_left_data, church_right_data = _testdata()["church"]
    axstereo = AxesStereo2D()
    axstereo.ax_left.imshow(church_left_data)
    axstereo.ax_right.imshow(church_right_data)
    axstereo.wiggle(wiggle_filepath)
    assert wiggle_filepath.exists()
    wiggle_filepath.unlink()  # remove the file
    assert True


def test_wiggle_frames():
    from PIL import Image
    from mpl_stereo.AxesStereo import boomerang_sequence

    assert boomerang_sequence(2) == [0, 1]
    assert boomerang_sequence(4) == [0, 1, 2, 3, 2, 1]

    wiggle_filepath = Path("test.gif")
    x, y, z = _testdata()["trefoil"]

    # Both 2D and 3D plotted data support an arbitrary number of frames, played
    # as a boomerang loop.
    for cls in (AxesStereo2D, AxesStereo3D):
        for frames in (2, 5):
            axstereo = cls()
            axstereo.plot(x, y, z)
            axstereo.wiggle(wiggle_filepath, frames=frames)
            assert wiggle_filepath.exists()
            with Image.open(wiggle_filepath) as im:
                assert im.n_frames == len(boomerang_sequence(frames))
            wiggle_filepath.unlink()

    # Image-based wiggles have no intermediate viewpoints, so frames > 2 errors
    church_left_data, church_right_data = _testdata()["church"]
    axstereo = AxesStereo2D()
    axstereo.ax_left.imshow(church_left_data)
    axstereo.ax_right.imshow(church_right_data)
    with pytest.raises(NotImplementedError):
        axstereo.wiggle(wiggle_filepath, frames=4)

    # frames must be an integer >= 2
    with pytest.raises(ValueError):
        AxesStereo3D().wiggle(wiggle_filepath, frames=1)

    assert not wiggle_filepath.exists()


def test_wiggle_images():
    from PIL import Image
    from mpl_stereo.AxesStereo import boomerang_sequence

    wiggle_filepath = Path("test.gif")

    # Build several synthetic viewpoints by shifting a gradient.
    base = np.tile(np.linspace(0, 1, 60), (40, 1))
    views = [np.roll(base, k, axis=1) for k in range(4)]

    axstereo = AxesStereo2D()
    axstereo.imshow_stereo(views, cmap="gray")
    assert len(axstereo.wiggle_images) == 4

    # frames=None wiggles through every supplied image.
    axstereo.wiggle(wiggle_filepath, frames=None)
    with Image.open(wiggle_filepath) as im:
        assert im.n_frames == len(boomerang_sequence(4))
    wiggle_filepath.unlink()

    # An integer frames count samples a subset of the images.
    axstereo = AxesStereo2D()
    axstereo.imshow_stereo(views, cmap="gray")
    axstereo.wiggle(wiggle_filepath, frames=3)
    with Image.open(wiggle_filepath) as im:
        assert im.n_frames == len(boomerang_sequence(3))
    wiggle_filepath.unlink()

    # imshow_stereo requires at least two images.
    with pytest.raises(ValueError):
        AxesStereo2D().imshow_stereo([base])


def test_save_plot_area(tmp_path):
    from PIL import Image

    x, y, z = _testdata()["trefoil"]
    dpi = 100

    # Anaglyph: a normal save keeps decorations, plot_area=True crops to the
    # plot area and leaves the original figure untouched.
    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z)
    axstereo.fig.set_size_inches(3, 3)

    axstereo.save(tmp_path / "full.png", dpi=dpi)
    with Image.open(tmp_path / "full.png") as im:
        assert im.size == (300, 300)  # exact: 3in * 100dpi

    # The cropped plot-area sizes come from the axes' rendered extent, which can
    # vary by a pixel across platforms (fonts/freetype), so allow a 1px margin.
    axstereo.save(tmp_path / "area.png", plot_area=True, dpi=dpi)
    with Image.open(tmp_path / "area.png") as im:
        assert tuple(im.size) == pytest.approx((232, 231), abs=1)
    assert axstereo.ax.axison  # original figure untouched

    # Side-by-side 2D: cropped to the two-panel plot area, with no gap between.
    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z)
    axstereo.fig.set_size_inches(6, 3)
    axstereo.save(tmp_path / "sbs.png", plot_area=True, dpi=dpi)
    with Image.open(tmp_path / "sbs.png") as im:
        assert tuple(im.size) == pytest.approx((462, 231), abs=1)

    # plot_area is not supported for 3D static plots.
    axstereo = AxesStereo3D()
    axstereo.plot(x, y, z)
    with pytest.raises(NotImplementedError):
        axstereo.save(tmp_path / "no.png", plot_area=True)


def test_save_animate_dispatch(tmp_path):
    from PIL import Image

    x, y, z = _testdata()["trefoil"]

    # A .png extension saves a static image (one frame).
    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z)
    axstereo.save(tmp_path / "static.png", dpi=100)
    with Image.open(tmp_path / "static.png") as im:
        assert getattr(im, "n_frames", 1) == 1

    # A .gif extension animates a wiggle.
    axstereo.save(tmp_path / "anim.gif", dpi=100)
    with Image.open(tmp_path / "anim.gif") as im:
        assert im.n_frames == 2


def test_save_wiggle_plot_area(tmp_path):
    from PIL import Image

    x, y, z = _testdata()["trefoil"]

    # 2D plotted wiggle with plot_area keeps sensible (landscape) proportions
    # rather than the narrow side-by-side panel.
    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z)
    axstereo.save(tmp_path / "plotted.gif", plot_area=True, dpi=100)
    with Image.open(tmp_path / "plotted.gif") as im:
        assert im.n_frames == 2
        assert im.size == (640, 480)

    # Image wiggle of a square image yields a square frame (no letterboxing).
    square = np.tile(np.linspace(0, 1, 50), (50, 1))
    axstereo = AxesStereo2D()
    axstereo.imshow_stereo([square, np.roll(square, 5, axis=1)], cmap="gray")
    axstereo.save(tmp_path / "image.gif", plot_area=True, dpi=100)
    with Image.open(tmp_path / "image.gif") as im:
        assert im.size == (400, 400)

    # plot_area is not supported for 3D wiggles.
    axstereo = AxesStereo3D()
    axstereo.plot(x, y, z)
    with pytest.raises(NotImplementedError):
        axstereo.save(tmp_path / "no.gif", plot_area=True)


## The following tests are for visual inspection only
def plotting_tests_2d_pairwise():
    # test plot and scatter
    x, y, z = _testdata()["trefoil"]
    axstereo = AxesStereo2D(eye_balance=-1.0)
    axstereo.plot(x, y, z, c="k", alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10)
    axstereo.grid(True)
    axstereo.set_title("eye_balance=-1")

    axstereo = AxesStereo2D(eye_balance=1.0)
    axstereo.plot(x, y, z, c="k", alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10)
    axstereo.grid(True)
    axstereo.set_title("eye_balance=1")

    # test bar and stem
    x = y = z = np.arange(10)
    axstereo = AxesStereo2D()
    axstereo.bar(x, y, z, width=1, edgecolor="white", linewidth=0.7)
    axstereo.stem(x, y, z, "k")


def plotting_tests_2d_pairwise_zlim():
    # test plot and scatter
    x, y, z = _testdata()["trefoil"]
    axstereo = AxesStereo2D(eye_balance=-1.0)
    axstereo.plot(x, y, z, c="k", alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10)
    axstereo.grid(True)
    axstereo.set_zlim((-0.5, 0.5))
    axstereo.set_title("zlim=(-0.5, 0.5)")

    axstereo = AxesStereo2D(eye_balance=1.0)
    axstereo.plot(x, y, z, c="k", alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10)
    axstereo.grid(True)
    axstereo.set_zlim((-2, 2))
    axstereo.set_title("zlim=(-2, 2)")


def plotting_tests_anaglyph_pairwise_zzero():
    # test plot and scatter
    x, y, z = _testdata()["trefoil"]
    axstereo = AxesAnaglyph(zzero=min(z))
    axstereo.plot(x, y, z, c="k", alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10)
    axstereo.grid(True)
    axstereo.set_title("zzero=min(zzero)")

    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z, c="k", alpha=0.2, zzero=None)
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10, zzero=None)
    axstereo.grid(True)
    axstereo.set_title("zzero=None")

    axstereo = AxesAnaglyph(zscale=np.ptp(z))
    axstereo.plot(x, y, z, c="k", alpha=0.2, zzero=max(z))
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10, zzero=max(z))
    axstereo.grid(True)
    axstereo.set_title("zzero=max(zzero)")


def plotting_tests_3d():
    x, y, z = _testdata()["trefoil"]
    axstereo = AxesStereo3D()
    axstereo.plot(x, y, z, c="k", alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10)


def plotting_tests_anaglyph_pairwise():
    x, y, z = _testdata()["trefoil"]
    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z, c="k", alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10)

    axstereo = AxesAnaglyph(eye_balance=1.0)
    axstereo.plot(x, y, z, c="k", alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap="viridis", s=10)


def plotting_tests_anaglyph_imshow_stereo():
    sun_left_data, sun_right_data = _testdata()["sun"]
    for cmap in [None, "Oranges_r", "viridis"]:
        axstereo = AxesAnaglyph()
        axstereo.imshow_stereo([sun_left_data, sun_right_data], cmap=cmap)
        axstereo.set_title(cmap)

    church_left_data, church_right_data = _testdata()["church"]
    for method in ["dubois", "photoshop", "photoshop2"]:
        axstereo = AxesAnaglyph()
        axstereo.imshow_stereo([church_left_data, church_right_data], method=method)
        axstereo.set_title(method)
        axstereo.fig.set_size_inches(4, 3)

    church_left_data, church_right_data = _testdata()["church_cropped"]
    axstereo = AxesAnaglyph()
    axstereo.imshow_stereo([church_left_data, church_right_data], crop=True)
    axstereo.set_title("cropped")


def plotting_tests_wiggle(filepath):
    church_left_data, church_right_data = _testdata()["church"]
    axstereo = AxesStereo2D()
    axstereo.ax_left.imshow(church_left_data)
    axstereo.ax_right.imshow(church_right_data)
    axstereo.wiggle(filepath, ax=None)


if __name__ == "__main__":
    filepath = Path("test.gif")

    plotting_tests_2d_pairwise()
    plotting_tests_2d_pairwise_zlim()
    plotting_tests_3d()
    plotting_tests_anaglyph_pairwise()
    plotting_tests_anaglyph_pairwise_zzero()
    plotting_tests_anaglyph_imshow_stereo()
    plotting_tests_wiggle(filepath)

    plt.show()

    filepath.unlink()  # remove the file
