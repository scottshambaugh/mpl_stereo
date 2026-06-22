import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import inspect
import copy
import warnings
import pickle
import io

from abc import ABC
from typing import Optional, Union, Any
from types import MethodType
from pathlib import Path
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d.axes3d import Axes3D, _Quaternion


## Functions
def sort_by_z(x: np.ndarray, y: np.ndarray, z: np.ndarray, kwargs: dict[str, Any]):
    """
    Sorts the provided data arrays based on the z values, to avoid improper
    occlusion in visualizations.

    Parameters
    ----------
    x : np.ndarray
        An array of x-coordinates for the data points.
    y : np.ndarray
        An array of y-coordinates for the data points.
    z : np.ndarray
        An array of z-coordinates for the data points.
        Data will be sorted based on these values.
    kwargs : dict[str, Any], optional
        A dictionary of additional keyword arguments.
        If it contains a 'c' key with an array of the same shape as z,
        this array will also be sorted alongside x, y, and z.
    """
    sort_idx = np.argsort(z)
    x = x[sort_idx]
    y = y[sort_idx]
    z = z[sort_idx]
    if "c" in kwargs and np.array(kwargs["c"]).shape == np.array(z).shape:
        c = kwargs.pop("c")
        c = c[sort_idx]
        kwargs["c"] = c
    return x, y, z, kwargs


def process_fill_between_args(args: Any, kwargs: dict[str, Any]):
    """
    Process the arguments to a fill_between call. fill_between has an
    (x, y1, y2) signature, so the depth must come from a separate `z` keyword.
    y1 is the value the offset is applied around, and y2 (plus anything else)
    stays in args as the second fill boundary.

    Parameters
    ----------
    args : Any
        The positional arguments passed to fill_between.
    kwargs : dict[str, Any]
        The keyword arguments passed to fill_between.
    """
    if "x" in kwargs:
        x = kwargs.pop("x")
    else:
        x, *args = args
    if "y1" in kwargs:
        y = kwargs.pop("y1")
    elif "y" in kwargs:  # reconstructed call during a redraw
        y = kwargs.pop("y")
    else:
        y, *args = args
    z = kwargs.pop("z", None)
    if z is None:
        warnings.warn(
            "fill_between needs the depth data passed as a `z` keyword argument "
            "to produce a stereoscopic effect, e.g. fill_between(x, y1, y2, z=z)."
        )
    return x, y, z, args, kwargs


def process_args(ax_method: Any, known_methods: list[str], args: Any, kwargs: dict[str, Any]):
    """
    Process the arguments to a method call to determine if the method is
    plotting x-y data and if there is a z argument or keyword argument.

    Parameters
    ----------
    ax_method : Any
        The matplotlib axes method for which the arguments are being processed.
    known_methods : list[str]
        A list of method names that are known to plot x-y data.
    args : Any
        The positional arguments passed to the ax_method.
    kwargs : dict[str, Any]
        The keyword arguments passed to the ax_method.
    """
    x = y = z = None

    if ax_method.__name__ == "fill_between":
        return process_fill_between_args(args, kwargs)

    parameters = inspect.signature(ax_method).parameters
    if "x" in kwargs:
        x = kwargs.pop("x")
    elif "x" in parameters or ax_method.__name__ in known_methods:
        x, *args = args
    if "y" in kwargs:
        y = kwargs.pop("y")
    elif "y" in parameters or ax_method.__name__ in known_methods:
        y, *args = args

    # Check if 'z' is in the keyword arguments or if there is a third
    # argument of the same shape as x
    if "z" in kwargs:
        z = kwargs.pop("z")
    elif len(args) > 0 and (np.array(args[0]).shape == np.array(x).shape):
        z, *args = args
    return x, y, z, args, kwargs


def calc_2d_offsets(
    eye_balance: float,
    z: np.ndarray,
    d: float,
    ipd: float,
    zautoscale: bool = True,
    zscale: Optional[float] = None,
    zlim: Optional[tuple[float, float]] = None,
    zzero: Optional[float] = None,
    xlim: Optional[tuple[float, float]] = None,
):
    """
    Calculates the x-offsets for a 2D plot to create a stereoscopic effect
    based on the z-coordinates of the data points.

    Parameters
    ----------
    eye_balance : float
        The eye balance, a value between -1 and 1.
    z : np.ndarray
        An array of z-coordinates for the data points.
    d : float
        The distance from the focal plane to the viewer, in millimeters.
    ipd : float
        The interpupillary distance, in millimeters.
    zscale : Optional[float]
        Scaling factor for the visual depth of the z-data (in x-axis units).
        If None (default), then will be set to 1/4 the x-axis range.
    zlim : Optional[tuple[float, float]]
        Z-axis limits. If None (default), then z range will be autoscaled to
        (min(z), max(z)).
    zzero : Optional[float]
        The z-coordinate of the focal plane. Set to min(z) to have all the data
        float above the page, or set to max(z) to have all the data sink into
        the page. If None (default), will be set to the midpoint of the z range.
    xlim : Optional[tuple[float, float]]
        The x-axis limits, for calculating the zscale. If None, then an
        undefined zscale will be set to 1.
    """
    if zautoscale or zscale is None:
        if zautoscale and xlim is not None:
            zscale = (max(xlim) - min(xlim)) / 4
        else:
            zscale = 1

    if zautoscale or zlim is None:
        zlim_new = (np.min(z), np.max(z))
        if zlim is None:
            zlim = zlim_new
        else:
            zlim = (min(zlim[0], zlim_new[0]), max(zlim[1], zlim_new[1]))

    zrange = zlim[1] - zlim[0]
    if zrange == 0:  # If all the z values are the same
        zrange = np.max(abs(z))
        if zrange == 0:
            zrange = 1
    z_midpoint = zlim[0] + zrange / 2
    if zzero is None:
        zzero = z_midpoint

    zscaled = (z + z_midpoint - zzero) / zrange * zscale
    offset = ipd * zscaled / (d + zscaled)
    offset_left = (eye_balance + 1) / 2 * offset
    offset_right = (1 - eye_balance) / 2 * offset
    return offset_left, offset_right, zlim, zscale


def apply_x_offset(ax: Axes, x: np.ndarray, offset: np.ndarray) -> np.ndarray:
    """
    Shift ``x`` by ``offset`` within the x-axis's scaled (transformed) space, so
    the parallax is screen-uniform for log and other non-linear x scales. For a
    linear axis the transform is the identity, giving the original x + offset.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes whose x-axis scale (transform) the offset is applied within.
    x : np.ndarray
        The x-coordinates of the data, in data units.
    offset : np.ndarray
        The horizontal offset to apply, in the x-axis's transformed units.
    """
    tf = ax.xaxis.get_transform()
    return tf.inverted().transform(tf.transform(np.asarray(x, dtype=float)) + offset)


def color_stem_container(container: Any, color: str, alpha: float):
    """
    Color the marker line, stem lines, and baseline of a stem plot's container.
    Used to apply an anaglyph eye color to a stem plot, which takes no ``color``
    keyword of its own (it uses ``linefmt``/``markerfmt`` format strings).

    Parameters
    ----------
    container : matplotlib.container.StemContainer
        The container returned by a ``stem`` call.
    color : str
        The color to apply.
    alpha : float
        The transparency to apply.
    """
    for artist in (container.markerline, container.stemlines, container.baseline):
        artist.set_color(color)
        artist.set_alpha(alpha)


def plot_anaglyph_eye(
    ax: Axes,
    name: str,
    x: np.ndarray,
    offset: np.ndarray,
    y: np.ndarray,
    color: str,
    alpha: float,
    args: Any,
    kwargs: dict[str, Any],
) -> Any:
    """
    Draw one eye's copy of the data for an anaglyph: offset x by ``offset`` in
    the x-axis's scaled space and draw it in ``color``. Most methods take a
    ``color`` keyword, but ``stem`` does not (it uses format strings), so it is
    drawn plainly and its returned container colored afterwards.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to draw on.
    name : str
        The name of the plotting method.
    x, y : np.ndarray
        The data coordinates.
    offset : np.ndarray
        The horizontal offset to apply, in the x-axis's transformed units.
    color : str
        The eye's color.
    alpha : float
        The transparency to apply.
    args, kwargs : Any
        Additional arguments forwarded to the plotting method.
    """
    x_offset = apply_x_offset(ax, x, offset)
    if name == "stem":
        # stem() takes no color/alpha keywords, so color the container after.
        res = getattr(ax, name)(x_offset, y, *args, **kwargs)
        color_stem_container(res, color, alpha)
    else:
        res = getattr(ax, name)(x_offset, y, color=color, alpha=alpha, *args, **kwargs)
    return res


def log_scale_kwargs(kwargs: dict[str, Any], axis: str) -> dict[str, Any]:
    """
    Pick the log-scale keyword arguments for one axis out of ``kwargs``, to
    forward to ``set_xscale``/``set_yscale``. Mirrors how matplotlib's
    semilogx/semilogy/loglog split their keywords: the shared ``base``,
    ``subs``, and ``nonpositive`` apply to either axis, plus the axis-suffixed
    aliases (e.g. ``basex``/``basey``).

    Parameters
    ----------
    kwargs : dict[str, Any]
        The keyword arguments passed to the plotting method.
    axis : str
        Which axis to collect keywords for, "x" or "y".
    """
    keys = ("base", "subs", "nonpositive", f"base{axis}", f"subs{axis}", f"nonpos{axis}")
    return {k: v for k, v in kwargs.items() if k in keys}


def view_init(self, elev=None, azim=None, roll=None, vertical_axis="z", share=False):
    """
    Overrides the Axes3D.view_init method to link the views of the left and
    right 3D axes while maintaining the correct offset. See that method's
    documentation for more information.
    """

    self._dist = 10  # The camera distance from origin. Behaves like zoom

    if elev is None:
        elev = self.initial_elev
    if azim is None:
        azim = self.initial_azim
    if roll is None:
        roll = self.initial_roll
    _vertical_axes = {"x": 0, "y": 1, "z": 2}
    try:
        vertical_axis = _vertical_axes[vertical_axis]
    except KeyError:
        raise ValueError(
            f"{vertical_axis!r} is not a valid value for vertical_axis. "
            f"Supported values are {', '.join(map(repr, _vertical_axes))}"
        ) from None

    if share:
        axes = {sibling for sibling in self._shared_axes["view"].get_siblings(self)}
        # Ensure self is at the beginning
        axes.remove(self)
        axes = [self] + list(axes)
    else:
        axes = [self]

    for ax in axes:
        if hasattr(ax, "stereo_offset") and hasattr(self, "stereo_offset") and ax is not self:
            q = _Quaternion.from_cardan_angles(*np.deg2rad((elev, azim, roll)))
            th = np.deg2rad(self.stereo_offset - ax.stereo_offset)
            k = np.array([0, 0, 1])
            dq = _Quaternion(np.cos(th), k * np.sin(th))
            q = dq * q
            elev, azim, roll = np.rad2deg(q.as_cardan_angles())
        ax.elev = elev
        ax.azim = azim
        ax.roll = roll
        ax._vertical_axis = vertical_axis


def boomerang_sequence(n: int) -> list[int]:
    """
    Build a "boomerang" index sequence over ``n`` frames that rocks back and
    forth without repeating the endpoints, so the animation loops seamlessly.

    For example, ``n=2`` gives ``[0, 1]`` (the classic two-frame wiggle) and
    ``n=4`` gives ``[0, 1, 2, 3, 2, 1]``.

    Parameters
    ----------
    n : int
        The number of distinct frames.
    """
    return list(range(n)) + list(range(n - 2, 0, -1))


# File suffixes that `save` treats as animated wiggle stereograms.
ANIMATION_SUFFIXES = {".gif", ".mp4", ".mov", ".avi", ".mkv", ".webm", ".apng", ".html"}


def crop_fig_to_plot_area(fig: Figure, ax_indices: list[int]) -> Figure:
    """
    Modify ``fig`` in place so it contains only the plot area of the given axes:
    strip those axes' spines, ticks, and labels, drop every other axis, and pack
    the kept axes left-to-right with no gap or surrounding padding. Each axis
    keeps its original size (nothing is squeezed).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to crop. Usually a throwaway copy of a live figure.
    ax_indices : list[int]
        Indices into ``fig.axes`` of the axes whose plot area to keep, in
        left-to-right order.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The same figure, cropped to the plot area.
    """
    fig.canvas.draw()
    keep = [fig.axes[i] for i in ax_indices]
    for ax in list(fig.axes):
        if ax not in keep:
            fig.delaxes(ax)

    # Pack the kept axes left-to-right with no gaps, each keeping its own width
    # and height, aligned along a common row.
    ordered = sorted(zip(keep, [a.get_position() for a in keep]), key=lambda kp: kp[1].x0)
    total_w = sum(p.width for _, p in ordered)
    max_h = max(p.height for _, p in ordered)

    x = 0.0
    for ax, p in ordered:
        w = p.width / total_w
        ax.set_axis_off()
        ax.set_anchor("C")  # keep the (aspect-constrained) image centered in its slot
        ax.set_position((x, 0.0, w, p.height / max_h))
        x += w

    fw, fh = fig.get_size_inches()
    fig.set_size_inches(total_w * fw, max_h * fh)
    return fig


def crop_image_center(data: np.ndarray, shape: tuple[int, int]):
    """
    Crop an image array to a specified shape, trimming equally on each side
    such that the center of the array is kept.

    Parameters
    ----------
    data : np.ndarray
        The array to crop.
    shape : tuple[int, int]
        The shape to crop to.
    """
    ndim = len(data.shape)
    data = np.atleast_3d(data)
    x, y, z = data.shape
    if len(shape) == 2:
        shape = (shape[0], shape[1], 1)
    cropx, cropy, _ = shape
    if cropx > x or cropy > y:
        raise ValueError(f"Crop shape ({cropx}, {cropy}) is larger than data shape ({x}, {y})")
    startx = x // 2 - (cropx // 2)
    starty = y // 2 - (cropy // 2)
    cropped_data = data[startx : startx + cropx, starty : starty + cropy, :]
    if ndim == 2:
        cropped_data = cropped_data[:, :, 0]
    return cropped_data


def sanitize_data_left_right(
    data_left: np.ndarray, data_right: np.ndarray, crop: bool
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sanitize the data for the left and right images, ensuring that they have
    the same shape and dtype, and that they are in the range 0.0 - 1.0.

    Parameters
    ----------
    data_left : np.ndarray
        The data from the left image.
    data_right : np.ndarray
        The data for the right image.
    crop : bool
        Whether to crop the image to the minimum size of the two images,
        keeping the images centered.

    Returns
    -------
    data_left : np.ndarray
        The sanitized data for the left axes.
    data_right : np.ndarray
        The sanitized data for the right axes.
    """
    # Check that the data is valid
    data_left = np.array(data_left)
    data_right = np.array(data_right)
    if data_left.shape != data_right.shape:
        if crop:
            min_shape = tuple(np.min([data_left.shape, data_right.shape], axis=0).tolist())
            data_left = crop_image_center(data_left, min_shape)
            data_right = crop_image_center(data_right, min_shape)
        else:
            raise ValueError("data_left and data_right must have the same shape")
    if data_left.dtype != data_right.dtype:
        raise ValueError("data_left and data_right must have the same dtype")

    # Accept both 0.0 - 1.0 and 0 - 255 data, map to 0.0 - 1.0
    if np.issubdtype(data_left.dtype, np.integer):
        data_left = data_left.astype(float) / 255
        data_right = data_right.astype(float) / 255
    return data_left, data_right


def sanitize_images(images: list[np.ndarray], crop: bool) -> list[np.ndarray]:
    """
    Sanitize a sequence of images, ensuring that they all share the same shape
    and dtype, and that they are in the range 0.0 - 1.0.

    Parameters
    ----------
    images : list[np.ndarray]
        The images to sanitize.
    crop : bool
        Whether to crop the images to the minimum common size, keeping them
        centered. If False, mismatched shapes raise a ValueError.

    Returns
    -------
    images : list[np.ndarray]
        The sanitized images.
    """
    arrays = [np.array(image) for image in images]
    if len({a.shape for a in arrays}) > 1:
        if crop:
            min_shape = tuple(np.min([a.shape for a in arrays], axis=0).tolist())
            arrays = [crop_image_center(a, min_shape) for a in arrays]
        else:
            raise ValueError("all images must have the same shape")
    if len({a.dtype for a in arrays}) > 1:
        raise ValueError("all images must have the same dtype")

    # Accept both 0.0 - 1.0 and 0 - 255 data, map to 0.0 - 1.0
    if np.issubdtype(arrays[0].dtype, np.integer):
        arrays = [a.astype(float) / 255 for a in arrays]
    return arrays


## Classes
class AxesStereoBase(ABC):
    def __init__(
        self,
        eye_balance: float = -1,
        d: float = 350,
        ipd: float = 65,
        zscale: Optional[float] = None,
        zlim: Optional[tuple[float, float]] = None,
        zzero: Optional[float] = None,
        is_3d: bool = False,
    ):
        """
        Parameters
        ----------
        eye_balance : float
            The eye balance parameter, from -1 to 1. A value of -1 means the
            left plot will have accurate x-axis labels, and a value of 1 means
            the right plot will. For any other value, both plots will have
            inaccurate x-axis labels.
        d : float
            Distance from the focal plane to the viewer (in millimeters).
        ipd : float
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        zscale : Optional[float]
            Scaling factor for the visual depth of the z-data (in x-axis units).
            If None (default), then will be set to 1/4 the x-axis range.
        zlim : Optional[tuple[float, float]]
            Z-axis limits. If None (default), then z range will be autoscaled to
            (min(z), max(z)).
        zzero : Optional[float]
            The z-coordinate of the focal plane. Set to min(z) to have all the
            data float above the page, or set to max(z) to have all the data
            sink into the page. If None (default), will be set to the midpoint
            of the z range.
        is_3d : bool
            Whether the axes are 3D. Default is False.
        """
        self.eye_balance = eye_balance
        self.d = d
        self.ipd = ipd
        self.zscale = zscale
        self.zlim = zlim
        self.zzero = zzero
        self.is_3d = is_3d
        self.known_methods: list[str] = []

        self.zautoscale = True
        if self.zlim is not None or zscale is not None:
            self.zautoscale = False
        self.is_redrawing = False

        self.artists_left = []
        self.artists_right = []
        self.artist_args = []

        # All image frames supplied to imshow_stereo, used for wiggle animations,
        # along with the imshow arguments they were displayed with.
        self.wiggle_images: list[np.ndarray] = []
        self.wiggle_imshow_args: tuple = ()
        self.wiggle_imshow_kwargs: dict[str, Any] = {}

    def _copy_fig(self) -> Figure:
        """
        Return a deep copy of ``self.fig`` (via pickle), so it can be modified
        without affecting the live figure.
        """
        buf = io.BytesIO()
        pickle.dump(self.fig, buf)
        buf.seek(0)
        return pickle.load(buf)

    def save(
        self,
        filepath: Union[str, Path],
        plot_area: bool = False,
        animate: Optional[bool] = None,
        *args: Any,
        **kwargs: Any,
    ):
        """
        Save the stereogram to a file. This is the single entry point for both
        static images and animated wiggle stereograms.

        Whether to write an animation is chosen from the file extension (e.g.
        ``.gif`` or ``.mp4`` animate; ``.png`` does not), or forced with the
        ``animate`` argument.

        Parameters
        ----------
        filepath : str | pathlib.Path
            The filepath to save to. The extension selects the format.
        plot_area : bool
            If True, save only the plot area: strip all spines, ticks, and
            labels and remove the surrounding padding. Default is False.
        animate : Optional[bool]
            Force (True) or forbid (False) writing a wiggle animation. If None
            (default), inferred from the file extension.
        *args : Any
            Additional positional arguments. For animations, passed to
            `animation.save`; for static images, passed to `Figure.savefig`.
        **kwargs : Any
            For animations, the wiggle options (``interval``, ``frames``, ...,
            see `wiggle`) plus any `animation.save` arguments. For static images,
            arguments passed to `Figure.savefig`.
        """
        filepath = Path(filepath)
        if animate is None:
            animate = filepath.suffix.lower() in ANIMATION_SUFFIXES

        if animate:
            self._run_wiggle(filepath, *args, plot_area=plot_area, **kwargs)
        elif plot_area:
            self._plot_area_fig().savefig(filepath, *args, **kwargs)
        else:
            self.fig.savefig(filepath, *args, **kwargs)

    def _run_wiggle(self, filepath: Path, *args: Any, **kwargs: Any):
        """Animate a wiggle stereogram. Only supported for side-by-side plots."""
        raise NotImplementedError(f"{type(self).__name__} cannot be animated")

    def _plot_area_fig(self) -> Figure:
        """Return a copy of the figure cropped to its plot area(s)."""
        raise NotImplementedError(f"plot_area is not supported for {type(self).__name__}")

    def log_artists(
        self, res_left: Any, res_right: Any, name: str, args: Any, kwargs: dict[str, Any]
    ):
        """
        Log artists in each of the self.artists_left and self.artists_right
        lists, and log the arguments in self.artist_args.

        Arguments
        ---------
        res_left : Any
            The return value of the left axes method call.
        res_right : Any
            The return value of the right axes method call.
        name : str
            Name of the plotting method.
        args: Any
            Arguments to the plotting method.
        kwargs: dict[str, Any]
            Keyword arguments to the plotting method.
        """
        if isinstance(res_left, list):
            self.artists_left.extend(res_left)
        else:
            self.artists_left.append(res_left)
        if isinstance(res_right, list):
            self.artists_right.extend(res_right)
        else:
            self.artists_right.append(res_right)

        self.artist_args.append((name, args, kwargs))


class AxesStereoSideBySide(AxesStereoBase):
    def __init__(
        self,
        fig: Optional[Figure] = None,
        axs: Optional[Union[tuple[Axes, Axes], tuple[Axes3D, Axes3D]]] = None,
        eye_balance: float = -1,
        d: float = 350,
        ipd: float = 65,
        zscale: Optional[float] = None,
        zlim: Optional[tuple[float, float]] = None,
        zzero: Optional[float] = None,
        is_3d: bool = False,
    ):
        """
        Parameters
        ----------
        fig : matplotlib.figure.Figure, optional
            The figure object to plot on.
        axs : tuple[matplotlib.axes.Axes, matplotlib.axes.Axes], optional
            The axes objects to plot on (ax_left, ax_right).
        eye_balance : float
            The eye balance parameter, from -1 to 1. A value of -1 means the
            left plot will have accurate x-axis labels, and a value of 1 means
            the right plot will. For any other value, both plots will have
            inaccurate x-axis labels.
        d : float
            Distance from the focal plane to the viewer (in millimeters).
        ipd : float
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        zscale : Optional[float]
            Scaling factor for the visual depth of the z-data (in x-axis units).
            If None (default), then will be set to 1/4 the x-axis range.
        zlim : Optional[tuple[float, float]]
            Z-axis limits. If None (default), then z range will be autoscaled to
            (min(z), max(z)).
        zzero : Optional[float]
            The z-coordinate of the focal plane. Set to min(z) to have all the
            data float above the page, or set to max(z) to have all the data
            sink into the page. If None (default), will be set to the midpoint
            of the z range.
        is_3d : bool
            Whether the axes are 3D. Default is False.
        """
        super().__init__(
            eye_balance=eye_balance,
            d=d,
            ipd=ipd,
            zscale=zscale,
            zlim=zlim,
            zzero=zzero,
            is_3d=is_3d,
        )

        # Generate two side-by-side subplots
        if fig is None and axs is None:
            if not is_3d:
                fig, axs = plt.subplots(1, 2)
            else:
                fig, axs = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
            self.ax_left = axs[0]
            self.ax_right = axs[1]
        elif axs is None:
            if not is_3d:
                self.ax_left = fig.add_subplot(121)
                self.ax_right = fig.add_subplot(122)
            else:
                self.ax_left = fig.add_subplot(121, projection="3d")
                self.ax_right = fig.add_subplot(122, projection="3d")
        else:
            fig = axs[0].figure
            self.ax_left = axs[0]
            self.ax_right = axs[1]

        self.ax_left.sharex(self.ax_right)
        self.ax_left.sharey(self.ax_right)

        self.fig = fig
        self.axs = (self.ax_left, self.ax_right)

    def wiggle(
        self,
        filepath: Union[str, Path],
        interval: float = 125,
        frames: Optional[int] = 2,
        ax: Optional[Axes] = None,
        yaxis_off: bool = False,
        plot_area: bool = False,
        *args: Any,
        **kwargs: Any,
    ):
        """
        Save the figure as a wiggle stereogram. This is a convenience wrapper,
        equivalent to ``save`` with ``animate=True``. See `_run_wiggle` for the
        meaning of each argument.
        """
        self._run_wiggle(
            Path(filepath), interval, frames, ax, yaxis_off, plot_area, *args, **kwargs
        )

    def _run_wiggle(
        self,
        filepath: Path,
        interval: float = 125,
        frames: Optional[int] = 2,
        ax: Optional[Axes] = None,
        yaxis_off: bool = False,
        plot_area: bool = False,
        *args: Any,
        **kwargs: dict[str, Any],
    ):
        """
        Write a wiggle-stereogram animation to ``filepath``.

        Parameters
        ----------
        filepath : pathlib.Path
            The filepath to save the animation to.
        interval : float
            The interval between frames in milliseconds, default 125.
        frames : Optional[int]
            The number of distinct viewpoints to sample across the stereo
            baseline, default 2. The animation rocks back and forth between the
            outermost views (a boomerang loop), so with ``frames > 2`` the
            intermediate viewpoints produce a smooth rocking motion. Any number
            of frames is supported for plotted 2D and 3D data. For image-based
            wiggles created with `imshow_stereo`, set ``frames=None`` to use
            all the supplied images.
        ax : matplotlib.axes.Axes, optional
            The target axes to plot the wiggle stereogram on. If None (default),
            then will plot on the axes of a new Figure.
        yaxis_off : bool
            Whether to hide the y-axis. Default is False.
        plot_area : bool
            If True, save only the plot area: strip all spines, ticks, and
            labels. When ``ax`` is None, the output figure is also sized to the
            content so it fills the frame without being squished. Only supported
            for 2D wiggles. Default is False.
        **kwargs : dict[str, Any]
            Additional keyword arguments passed to `animation.save`.
        """
        # frames=None means "use every available image" for image wiggles, and
        # falls back to the default of 2 for plotted data.
        if frames is None:
            frames = len(self.wiggle_images) if self.wiggle_images else 2
        if not isinstance(frames, int) or frames < 2:
            raise ValueError(f"frames must be an integer >= 2 or None, got {frames!r}")

        if plot_area and self.is_3d:
            raise NotImplementedError("plot_area=True is only supported for 2D wiggles")

        filepath = Path(filepath)
        if ax is None and plot_area:
            # Build a single full-bleed axis. Match the figure's aspect to the
            # content so the data fills the frame without being squished: for
            # images use the image's aspect, otherwise use the default figure
            # proportions (the same proportions a normal wiggle is drawn at).
            if self.wiggle_images:
                h, w = self.wiggle_images[0].shape[:2]
                base = 4.0
                fig = plt.figure(figsize=(base * w / max(w, h), base * h / max(w, h)))
            else:
                fig = plt.figure()
            ax_target = fig.add_axes((0, 0, 1, 1))
        elif ax is None:
            fig, ax_target = plt.subplots()
        else:
            ax_target = ax
            fig = ax_target.figure

        pos = ax_target.get_position()
        fig.delaxes(ax_target)

        if self.is_3d:
            ani = self._wiggle_3d(fig, pos, frames, interval, yaxis_off, kwargs)
        elif self.wiggle_images:
            ani = self._wiggle_images(fig, pos, frames, interval, yaxis_off, plot_area, kwargs)
        else:
            ani = self._wiggle_2d(fig, pos, frames, interval, yaxis_off, plot_area, kwargs)

        ani.save(filepath, *args, **kwargs)

    def _duplicate_fig(self, fig: Figure, kwargs: dict[str, Any]) -> Figure:
        """
        Return a deep copy of ``self.fig``, sized and dpi-matched to the target
        figure ``fig``.
        """
        fig_buffer = self._copy_fig()
        fig_buffer.set_size_inches(fig.get_size_inches())
        fig_buffer.set_dpi(kwargs.get("dpi", fig.get_dpi()))
        return fig_buffer

    def _attach_axis(
        self, ax: Axes, fig: Figure, pos: Any, yaxis_off: bool, plot_area: bool = False
    ):
        """
        Move a duplicated axis from a pickled buffer figure onto the target
        figure ``fig``, sizing it to fill the position ``pos``. When
        ``plot_area`` is True, strip all of the axis decorations so only the
        plot area remains.
        """
        ax._parent_figure = None
        ax.figure = fig
        fig.axes.append(ax)
        fig.add_axes(ax)
        ax.set_position(pos.bounds)
        if plot_area:
            ax.set_axis_off()
        elif yaxis_off:
            ax.yaxis.set_visible(False)

    def _wiggle_2d(
        self,
        fig: Figure,
        pos: Any,
        frames: int,
        interval: float,
        yaxis_off: bool,
        plot_area: bool,
        kwargs: dict[str, Any],
    ) -> FuncAnimation:
        """
        Build a 2D wiggle animation. For ``frames == 2`` the duplicated left and
        right axes are toggled on a shared position. For ``frames > 2``, the
        parallax is re-rendered at intermediate viewpoints by sweeping the
        ``eye_balance`` from -1 to +1 (which slides the left-eye axis across the
        full parallax range) and capturing the left axis at each step.
        """
        if frames == 2:
            # Fast path, byte-compatible with the original two-frame behavior.
            fig_buffer = self._duplicate_fig(fig, kwargs)
            axs = fig_buffer.axes[0:2]
            for ax in axs:
                self._attach_axis(ax, fig, pos, yaxis_off, plot_area)
                ax.set_visible(False)

            def update(frame):
                axs[frame].set_visible(True)
                return (axs[frame],)

            return FuncAnimation(fig, update, frames=2, interval=interval)

        # frames > 2: re-render the plotted data at intermediate viewpoints.
        if not self.artist_args:
            raise NotImplementedError(
                "frames > 2 is only supported for plotted 2D data. For image "
                "wiggles, pass all frames to `imshow_stereo` so they can be "
                "animated (intermediate viewpoints cannot be synthesized from a "
                "fixed pair of images)."
            )

        frame_axes = []
        original_eye_balance = self.eye_balance
        try:
            for eye_balance in np.linspace(-1, 1, frames):
                self.eye_balance = eye_balance
                self.redraw()
                fig_buffer = self._duplicate_fig(fig, kwargs)
                ax_wiggle = fig_buffer.axes[0]  # the left-eye axis
                self._attach_axis(ax_wiggle, fig, pos, yaxis_off, plot_area)
                ax_wiggle.set_visible(False)
                frame_axes.append(ax_wiggle)
        finally:
            # Leave the live figure as we found it.
            self.eye_balance = original_eye_balance
            self.redraw()

        sequence = boomerang_sequence(frames)

        def update(i):
            for k, ax in enumerate(frame_axes):
                ax.set_visible(k == sequence[i])
            return (frame_axes[sequence[i]],)

        return FuncAnimation(fig, update, frames=len(sequence), interval=interval)

    def _wiggle_3d(
        self,
        fig: Figure,
        pos: Any,
        frames: int,
        interval: float,
        yaxis_off: bool,
        kwargs: dict[str, Any],
    ) -> FuncAnimation:
        """
        Build an ``frames``-frame wiggle animation for a 3D plot by sweeping the
        azimuth of a single duplicated axis between the left- and right-eye
        views and rocking back and forth (a boomerang loop).
        """
        # A single faithful copy of the left-eye axis is enough: every frame is
        # just a different projection of the same 3D artists.
        fig_buffer = self._duplicate_fig(fig, kwargs)
        ax_wiggle = fig_buffer.axes[0]
        self._attach_axis(ax_wiggle, fig, pos, yaxis_off)

        # Sweep the azimuth between the two eye views.
        azims = np.linspace(self.ax_left.azim, self.ax_right.azim, frames)
        elev, roll = ax_wiggle.elev, ax_wiggle.roll
        sequence = boomerang_sequence(frames)

        def update(i):
            ax_wiggle.view_init(elev=elev, azim=azims[sequence[i]], roll=roll)
            return (ax_wiggle,)

        return FuncAnimation(fig, update, frames=len(sequence), interval=interval)

    def imshow_stereo(
        self,
        images: list[np.ndarray],
        cmap: Optional[str] = None,
        crop: bool = False,
        *args: Any,
        **kwargs: Any,
    ) -> tuple[Any, Any]:
        """
        Display existing stereo images side-by-side. A list of two or more images
        may be passed; the first is shown on the left axes and the last on the
        right axes. All of the images are kept for use by `wiggle` (see its
        ``frames`` argument). Any further args or kwargs are passed to `imshow`.

        Parameters
        ----------
        images : list[numpy.ndarray]
            The stereo images, ordered from the left-eye view to the right-eye
            view. Each may be of shape (M, N), (M, N, 3), or (M, N, 4).
        cmap : str, optional
            The colormap to use, default None. Recommended for (M, N) scalar
            data rather than color images.
        crop : bool
            Whether to crop the images to the minimum common size, keeping them
            centered. Default is False.
        """
        if len(images) < 2:
            raise ValueError("imshow_stereo requires a list of at least two images")

        self.wiggle_images = sanitize_images(list(images), crop)
        if cmap is not None:
            kwargs["cmap"] = cmap
        self.wiggle_imshow_args = args
        self.wiggle_imshow_kwargs = dict(kwargs)

        res_left = self.ax_left.imshow(self.wiggle_images[0], *args, **kwargs)
        res_right = self.ax_right.imshow(self.wiggle_images[-1], *args, **kwargs)
        return (res_left, res_right)

    def _wiggle_images(
        self,
        fig: Figure,
        pos: Any,
        frames: int,
        interval: float,
        yaxis_off: bool,
        plot_area: bool,
        kwargs: dict[str, Any],
    ) -> FuncAnimation:
        """
        Build a wiggle animation from the images supplied to `imshow_stereo`. The
        images are sampled evenly on a single axis as a boomerang loop.
        """
        n = len(self.wiggle_images)
        if frames >= n:
            selected = self.wiggle_images
        else:
            idx = np.unique(np.linspace(0, n - 1, frames).round().astype(int))
            selected = [self.wiggle_images[i] for i in idx]

        fig_buffer = self._duplicate_fig(fig, kwargs)
        ax_wiggle = fig_buffer.axes[0]  # left-eye axis
        self._attach_axis(ax_wiggle, fig, pos, yaxis_off, plot_area)

        # Replace the displayed image with one artist per sampled frame.
        for image in list(ax_wiggle.images):
            image.set_visible(False)
        artists = [
            ax_wiggle.imshow(image, *self.wiggle_imshow_args, **self.wiggle_imshow_kwargs)
            for image in selected
        ]
        for artist in artists:
            artist.set_visible(False)

        sequence = boomerang_sequence(len(artists))

        def update(i):
            for k, artist in enumerate(artists):
                artist.set_visible(k == sequence[i])
            return (artists[sequence[i]],)

        return FuncAnimation(fig, update, frames=len(sequence), interval=interval)


class AxesStereo2DBase(ABC):
    def set_zlim(
        self,
        zlim: tuple[float, float],
        zscale: Optional[float] = None,
        zautoscale: bool = False,
        redraw: Optional[bool] = None,
    ):
        """
        Set the z limits of both axes to the same value.

        Parameters
        ----------
        zlim : tuple[float, float]
            The new z limits.
        zscale : Optional[float]
            Scaling factor for the visual depth of the z-data (in x-axis units).
        zautoscale : bool
            Whether to later automatically scale the z limits based on the data.
            Default is False.
        redraw : bool
            Whether to redraw the plot. If None, then will redraw.
        """
        if redraw is None and (self.zlim != zlim or (self.zscale != zscale and zscale is not None)):
            redraw = True
        self.zlim = zlim
        if zscale is not None:
            self.zscale = zscale
        self.zautoscale = False  # Do not autoscale for this redraw
        if redraw:
            self.redraw()
        self.zautoscale = zautoscale

    def get_zlim(self) -> tuple[float, float]:
        """
        Return the z limit of the axes.
        """
        return self.zlim

    def autoscale_z(self):
        """
        Autoscale the z limit of both axes.
        """
        self.zautoscale = True
        self.zlim = self._calc_bounding_zlim()
        for _, _, kwargs in self.artist_args:
            kwargs.pop("zlim", None)
        self.redraw()

    def _calc_bounding_zlim(self) -> tuple[float, float]:
        """
        Calculate the z limits that will bound all the z data for all artists.
        """
        zlim = (np.inf, -np.inf)
        for _, args, kwargs in self.artist_args:
            # Check if 'z' is in the keyword arguments or if there is a third
            # argument of the same shape as x
            if "z" in kwargs:
                z = kwargs["z"]
            elif len(args) > 0:
                z, *args = args
            zlim = (min(zlim[0], np.min(z)), max(zlim[1], np.max(z)))
        return zlim

    def redraw(self):
        """
        Redraw the plot.
        """
        # Remove all the artists
        for artist in self.artists_left + self.artists_right:
            artist.remove()
        self.artists_left = []
        self.artists_right = []

        # Plot the data again
        self.is_redrawing = True
        artist_args = self.artist_args
        self.artist_args = []  # will repopulate in the getattr calls below
        for name, args, kwargs in artist_args:
            getattr(self, name)(*args, **kwargs)
        self.is_redrawing = False

    def set_xscale(self, *args: Any, **kwargs: Any):
        """
        Set the x-axis scale (e.g. "log") on all of the underlying axes and
        redraw, so the parallax offsets are recomputed for the new scale. The
        depth data is always treated as linear regardless of the x scale.
        """
        # Check __dict__ directly: __getattr__ returns a method for any missing
        # name, so a plain getattr would not fall through. 2D has .axs, the
        # anaglyph has a single .ax.
        axes = self.axs if "axs" in self.__dict__ else (self.ax,)
        for ax in axes:
            ax.set_xscale(*args, **kwargs)
        if self.artist_args:
            self.redraw()

    def semilogx(self, *args: Any, **kwargs: Any):
        """
        Stereo equivalent of `matplotlib.axes.Axes.semilogx`: log-scale the
        x-axis (forwarding *base*, *subs*, *nonpositive*) and then plot. The
        parallax is computed in the log-scaled x space, so it stays correct.
        """
        dx = log_scale_kwargs(kwargs, "x")
        self.set_xscale("log", **dx)
        return self.plot(*args, **{k: v for k, v in kwargs.items() if k not in dx})

    def loglog(self, *args: Any, **kwargs: Any):
        """
        Stereo equivalent of `matplotlib.axes.Axes.loglog`: log-scale both axes
        and then plot. The parallax is computed in the log-scaled x space, so it
        stays correct.
        """
        dx = log_scale_kwargs(kwargs, "x")
        self.set_xscale("log", **dx)
        dy = log_scale_kwargs(kwargs, "y")
        self.set_yscale("log", **dy)
        return self.plot(*args, **{k: v for k, v in kwargs.items() if k not in {*dx, *dy}})

    def plot2d(
        self,
        ax_left: Axes,
        ax_right: Axes,
        name: str,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        args: Any,
        kwargs: dict[str, Any],
    ) -> tuple[Any, Any]:
        """
        Plot the data twice, once for each eye view. This happens either on
        two subplots (for AxesStereo2D), or on the same subplot with different
        colors (for AxesAnaglyph).

        Parameters
        ----------
        ax_left : matplotlib.axes.Axes
            The left axes object.
        ax_right : matplotlib.axes.Axes
            The right axes object.
        name : str
            The name of the plotting method.
        x : np.ndarray
            An array of x-coordinates for the data points.
        y : np.ndarray
            An array of y-coordinates for the data points.
        z : np.ndarray
            An array of z-coordinates for the data points.
        args : Any
            The arguments passed to the plotting method.
        kwargs : dict[str, Any]
            The keyword arguments passed to the plotting method.

        Returns
        -------
        res_left : Any
            The return value of the left axes method call.
        res_right : Any
            The return value of the right axes method call.
        """
        # for scatter plots, sort the data by z to not occlude improperly
        kwargs_original = copy.deepcopy(kwargs)
        kwargs_original["x"] = x
        kwargs_original["y"] = y
        kwargs_original["z"] = z
        if name == "scatter":
            x, y, z, kwargs = sort_by_z(x, y, z, kwargs)

        # Extract the zzero and zscale keyword arguments if they exist
        zzero = kwargs.pop("zzero", None)
        if zzero is None and self.zzero is not None:
            zzero = self.zzero
        zscale = kwargs.pop("zscale", None)
        if zscale is None and self.zscale is not None:
            zscale = self.zscale

        # Extract the zlim keyword argument if it exists and update limits
        zlim = kwargs.pop("zlim", None)
        if zlim is not None:
            self.set_zlim(zlim, zscale, redraw=False)

        # Calculate the x-offsets. The offsets and zscale live in the x-axis's
        # scaled (transformed) space so the parallax is correct on log and other
        # non-linear x scales, so the x-range is transformed to match.
        xtransform = ax_left.xaxis.get_transform()
        xlim_scaled = xtransform.transform(np.asarray(ax_left.get_xlim(), dtype=float))
        offset_left, offset_right, zlim, zscale = calc_2d_offsets(
            self.eye_balance,
            z,
            self.d,
            self.ipd,
            self.zautoscale,
            zscale=zscale,
            zlim=self.zlim,
            zzero=zzero,
            xlim=xlim_scaled,
        )
        self.zscale = zscale
        self.zlim = zlim
        if len(self.artist_args) > 0 and not self.is_redrawing:
            self.redraw()

        if isinstance(self, AxesStereo2D):
            # Plot the data twice, once for each subplot (2D case)
            res_left = getattr(ax_left, name)(
                apply_x_offset(ax_left, x, offset_left), y, *args, **kwargs
            )
            res_right = getattr(ax_right, name)(
                apply_x_offset(ax_right, x, -offset_right), y, *args, **kwargs
            )
        elif isinstance(self, AxesAnaglyph):
            # Clear all color arguments (anaglyph case)
            kwargs.pop("c", None)
            kwargs.pop("color", None)
            kwargs.pop("cmap", None)
            kwargs.pop("alpha", None)

            # Plot the data twice, once for each color
            res_left = plot_anaglyph_eye(
                ax_left, name, x, offset_left, y, self.colors[1], self.alpha, args, kwargs
            )
            res_right = plot_anaglyph_eye(
                ax_right, name, x, -offset_right, y, self.colors[0], self.alpha, args, kwargs
            )

        # Keep track of the artists
        self.log_artists(res_left, res_right, name, args, kwargs_original)

        return res_left, res_right


class AxesStereo2D(AxesStereoSideBySide, AxesStereo2DBase):
    def __init__(
        self,
        fig: Optional[Figure] = None,
        axs: Optional[tuple[Axes, Axes]] = None,
        eye_balance: float = -1,
        d: float = 350,
        ipd: float = 65,
        zscale: Optional[float] = None,
        zlim: Optional[tuple[float, float]] = None,
        zzero: Optional[float] = None,
    ):
        """
        A class for creating stereoscopic 2D plots.

        Parameters
        ----------
        fig : matplotlib.figure.Figure, optional
            The figure object to plot on.
        axs : tuple[matplotlib.axes.Axes, matplotlib.axes.Axes], optional
            The axes objects to plot on (ax_left, ax_right).
        eye_balance : float
            The eye balance parameter, from -1 to 1. A value of -1 means the
            left plot will have accurate x-axis labels, and a value of 1 means
            the right plot will. For any other value, both plots will have
            inaccurate x-axis labels.
            All inaccurate axis labels will have transparency applied.
        d : float
            Distance from the focal plane to the viewer (in millimeters).
        ipd : float
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        zscale : Optional[float]
            Scaling factor for the visual depth of the z-data (in x-axis units).
            If None (default), then will be set to 1/4 the x-axis range.
        zlim : Optional[tuple[float, float]]
            Z-axis limits. If None (default), then z range will be autoscaled to
            (min(z), max(z)).
        zzero : Optional[float]
            The z-coordinate of the focal plane. Set to min(z) to have all the
            data float above the page, or set to max(z) to have all the data
            sink into the page. If None (default), will be set to the midpoint
            of the z range.
        """
        super().__init__(
            fig=fig,
            axs=axs,
            eye_balance=eye_balance,
            d=d,
            ipd=ipd,
            zlim=zlim,
            zscale=zscale,
            zzero=zzero,
            is_3d=False,
        )
        self.known_methods = [
            "plot",
            "step",
            "scatter",
            "stem",
            "bar",
            "fill",
            "fill_between",
            "errorbar",
            "text",
        ]

        # Minimize whitespace between plots
        self.fig.subplots_adjust(wspace=0.01)
        self.ax_right.tick_params(axis="y", length=0, labelcolor=(0, 0, 0, 0))

        # Give the innaccurate x-axis labels some transparency
        self.set_axlabel_alphas(alpha=0.5)

    def __getattr__(self, name: str) -> Any:
        """
        Delegate method calls to the left and right axes if the method is not
        defined in AxesStereoSideBySide. If the method has 'x' and 'y' as
        arguments, and either there is a third argument or 'z' is a keyword
        argument, then the z data will be used to offset the x data for the
        left and right axes and create the stereoscopic effect.

        Parameters
        ----------
        name : str
            The name of the attribute.
        """

        def method(*args: Any, **kwargs: dict[str, Any]) -> tuple[Any, Any]:
            """
            The method that will be called on the left and right axes. If the
            method has 'x' and 'y' as arguments, and either there is a third
            argument or 'z' is a keyword argument, then the z data will be used
            to offset the x data for the left and right axes and create the
            stereoscopic effect.

            Parameters
            ----------
            *args : Any
                The positional arguments passed to the method.
            **kwargs : dict[str, Any]
                The keyword arguments passed to the method.
            """
            ax_method = getattr(self.ax_left, name, None)
            args_original = args
            x, y, z, args, kwargs = process_args(ax_method, self.known_methods, args, kwargs)

            if all(var is not None for var in [ax_method, x, y, z]):
                res_left, res_right = self.plot2d(
                    self.ax_left, self.ax_right, name, x, y, z, args, kwargs
                )
            else:
                # For methods that don't plot x-y data
                res_left = getattr(self.ax_left, name)(*args_original, **kwargs)
                res_right = getattr(self.ax_right, name)(*args_original, **kwargs)
            return (res_left, res_right)

        return method

    def set_axlabel_alphas(self, alpha: float):
        """
        For axis labels that are not accurate to the plotted data, set their
        alpha to a value less than 1.

        Parameters
        ----------
        alpha : float
            Alpha value for the inaccurate axis labels.
        """
        if self.eye_balance != -1:
            for label in self.ax_left.get_xticklabels():
                label.set_alpha(alpha)
        elif self.eye_balance != 1:
            for label in self.ax_right.get_xticklabels():
                label.set_alpha(alpha)

    def _plot_area_fig(self) -> Figure:
        """Return a copy of the figure cropped to the two side-by-side panels."""
        return crop_fig_to_plot_area(self._copy_fig(), [0, 1])


class AxesStereo3D(AxesStereoSideBySide):
    def __init__(
        self,
        fig: Optional[Figure] = None,
        axs: Optional[tuple[Axes3D, Axes3D]] = None,
        eye_balance: float = -1,
        d: float = 350,
        ipd: float = 65,
    ):
        """
        A class for creating stereoscopic 3D plots.

        Parameters
        ----------
        fig : matplotlib.figure.Figure, optional
            The figure object to plot on.
        axs : tuple[mpl_toolkits.mplot3d.axes3d.Axes3D,
                    mpl_toolkits.mplot3d.axes3d.Axes3D], optional
            The axes3d objects to plot on (ax_left, ax_right).
        eye_balance : float
            The eye balance parameter, from -1 to 1. A value of -1 means the
            left plot will have accurate x-axis labels, and a value of 1 means
            the right plot will. For any other value, both plots will have
            inaccurate x-axis labels.
        d : float
            Distance from the focal plane to the viewer (in millimeters).
        ipd : float
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        """
        super().__init__(
            fig=fig,
            axs=axs,
            eye_balance=eye_balance,
            d=d,
            ipd=ipd,
            zscale=None,
            zlim=None,
            zzero=None,
            is_3d=True,
        )
        self.known_methods = [
            "plot",
            "plot3D",
            "scatter",
            "stem",
            "voxels",
            "plot_wireframe",
            "plot_surface",
            "plot_trisurf",
            "contour",
            "contourf",
        ]

        self.ax_left.sharez(self.ax_right)

        # Override the view_init method to link the views of the left and right
        # 3D axes while maintaining the correct offset
        self.ax_left.shareview(self.ax_right)
        self.ax_left.view_init = MethodType(view_init, self.ax_left)
        self.ax_right.view_init = MethodType(view_init, self.ax_right)

    def __getattr__(self, name: str):
        """
        Delegate method calls to the left and right axes if the method is not
        defined in AxesStereoSideBySide. If the method has 'x' and 'y' as
        arguments, and either there is a third argument or 'z' is a keyword
        argument, then the z data will be used to offset the x data for the
        left and right axes and create the stereoscopic effect.

        Parameters
        ----------
        name : str
            The name of the attribute.
        """

        def method(*args, **kwargs):
            """
            The method that will be called on the left and right axes. If the
            method has 'x' and 'y' as arguments, and either there is a third
            argument or 'z' is a keyword argument, then the z data will be used
            to offset the x data for the left and right axes and create the
            stereoscopic effect.

            Parameters
            ----------
            *args : Any
                The positional arguments passed to the method.
            **kwargs : dict[str, Any]
                The keyword arguments passed to the method.
            """
            # Reflect the method to check if 'x' and 'y' are in the arguments
            ax_method = getattr(self.ax_left, name, None)
            parameters = inspect.signature(ax_method).parameters

            is_plottable = False
            keyword_groups = [["x", "y", "z"], ["xs", "ys", "zs"], ["X", "Y", "Z"]]
            for keyword_group in keyword_groups:
                if all([keyword in parameters for keyword in keyword_group]):
                    is_plottable = True
            if name in ("plot", "plot3D"):
                # 'zs' is implicit and not in the parameters
                is_plottable = True

            if ax_method and is_plottable:
                offset_left, offset_right = self.calc_3d_offsets()

                # Set the views for both subplots
                self.ax_left.stereo_offset = -offset_left
                self.ax_right.stereo_offset = offset_right
                self.ax_left.view_init(azim=self.ax_left.azim - offset_left)
                self.ax_right.view_init(azim=self.ax_right.azim + offset_right)

                # Plot the data twice, once for each subplot
                res_left, res_right = self._plot_left_right(name, args, kwargs)

                # Keep track of the artists
                self.log_artists(res_left, res_right, name, args, kwargs)

            else:
                # For methods that do not involve 'x' and 'y'
                res_left = getattr(self.ax_left, name)(*args, **kwargs)
                res_right = getattr(self.ax_right, name)(*args, **kwargs)
            return (res_left, res_right)

        return method

    def calc_3d_offsets(self) -> tuple[float, float]:
        """
        Calculate the angular view offsets for a 3D plot.

        Returns
        -------
        offset_left : float
            The offset for the left subplot [deg].
        offset_right : float
            The offset for the right subplot [deg].
        """
        ang = 90 - np.rad2deg(np.arctan(2 * self.d / abs(self.ipd)))
        offset = ang / 2 * np.sign(self.ipd)
        offset_left = (self.eye_balance + 1) / 2 * offset
        offset_right = (1 - self.eye_balance) / 2 * offset
        return offset_left, offset_right

    def _plot_left_right(self, name: str, args: Any, kwargs: dict[str, Any]):
        """
        Draw the data on the left and right axes. Overridable so subclasses
        (e.g. AxesAnaglyph3D) can inject per-eye styling.
        """
        return (
            getattr(self.ax_left, name)(*args, **kwargs),
            getattr(self.ax_right, name)(*args, **kwargs),
        )


class AxesAnaglyph3D(AxesStereo3D):
    def __init__(
        self,
        fig: Optional[Figure] = None,
        axs: Optional[tuple[Axes3D, Axes3D]] = None,
        eye_balance: float = -1,
        d: float = 350,
        ipd: float = 65,
        colors: list[str, str] = ["red", "cyan"],
        alpha: float = 0.5,
    ):
        """
        A class for creating stereoscopic 3D anaglyph plots, viewed with red-cyan
        "3D glasses". Two 3D axes are stacked on top of each other at the correct
        left- and right-eye view offsets, and the data is drawn once per eye in
        each glasses-lens color. Any color arguments to plotting methods are
        ignored and replaced.

        Parameters
        ----------
        fig : matplotlib.figure.Figure, optional
            The figure object to plot on.
        axs : tuple[mpl_toolkits.mplot3d.axes3d.Axes3D,
                    mpl_toolkits.mplot3d.axes3d.Axes3D], optional
            The axes3d objects to plot on (ax_left, ax_right).
        eye_balance : float
            The eye balance parameter, from -1 to 1. A value of -1 means the
            left view will have accurate axis labels, and a value of 1 means the
            right view will.
        d : float
            Distance from the focal plane to the viewer (in millimeters).
        ipd : float
            Interpupillary distance (in millimeters). Default is 65. Note that
            for anaglyphs there is no cross-eyed viewing, so any negative ipds
            passed in will be changed to their absolute value.
        colors : list[str, str]
            Colors for the left and right glasses lenses. Default is
            ['red', 'cyan']. Because a lens color prevents that eye from seeing
            that color of data, each eye sees the opposite color. Eg.
            ['red', 'cyan'] means the left eye has a red lens and sees cyan.
        alpha : float
            The transparency of each eye's data, so the two overlaid views
            blend. Default is 0.5.
        """
        ipd = abs(ipd)  # anaglyphs are not cross-view
        super().__init__(fig=fig, axs=axs, eye_balance=eye_balance, d=d, ipd=ipd)
        self.colors = colors
        self.alpha = alpha

        # Stack the two eye axes on top of each other. When created standalone,
        # overlap them filling the figure. When handed existing axes (e.g. a
        # StereoSquare3D cell), the caller is responsible for overlapping them.
        # We leave their positions alone so 3D aspect handling stays consistent.
        if axs is None:
            for ax in (self.ax_left, self.ax_right):
                ax.set_position([0, 0, 1, 1])

        # Make both axes see-through (so the figure and the other eye show
        # through) and color every decoration in that eye's color, so the panes,
        # gridlines, ticks, and labels all become part of the anaglyph too.
        self.ax_left.patch.set_alpha(0)
        self.ax_right.patch.set_alpha(0)
        self._color_decorations(self.ax_left, self.colors[1])
        self._color_decorations(self.ax_right, self.colors[0])

    def _color_decorations(self, ax: Axes3D, color: str):
        """
        Color all of a 3D axis's decorations (panes, gridlines, axis lines,
        tick marks, tick labels, and axis labels) a single color, with the pane
        fill left transparent so the two overlaid views show through each other.
        """
        pane_color = mpl.colors.to_rgba(color, 0.08)  # faint pane shading
        grid_color = mpl.colors.to_rgba(color, 0.3)  # transparent so it recedes
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            # set_pane_color (not pane.set_facecolor) so the chosen alpha
            # survives the draw. A direct facecolor has its alpha reset to 0.5.
            axis.set_pane_color(pane_color)  # slight shading, mostly see-through
            axis.pane.set_edgecolor(color)  # colored box edges
            axis.line.set_color(color)
            axis.label.set_color(color)
            axis._axinfo["grid"]["color"] = grid_color
        ax.tick_params(colors=color)  # tick marks and tick labels, all three axes

    def _plot_left_right(self, name: str, args: Any, kwargs: dict[str, Any]):
        """
        Draw the data on each eye's axis in that eye's anaglyph color. Mirrors
        the color handling of AxesAnaglyph.plot2d: any user-supplied color
        arguments are dropped and replaced with the glasses-lens colors.
        """
        for key in ("c", "color", "cmap", "alpha"):
            kwargs.pop(key, None)
        res_left = getattr(self.ax_left, name)(
            *args, color=self.colors[1], alpha=self.alpha, **kwargs
        )
        res_right = getattr(self.ax_right, name)(
            *args, color=self.colors[0], alpha=self.alpha, **kwargs
        )
        return (res_left, res_right)


class AxesAnaglyph(AxesStereoBase, AxesStereo2DBase):
    def __init__(
        self,
        fig: Optional[Figure] = None,
        ax: Optional[Axes] = None,
        eye_balance: float = -1,
        d: float = 350,
        ipd: float = 65,
        zscale: Optional[float] = None,
        zlim: Optional[tuple[float, float]] = None,
        zzero: Optional[float] = None,
        colors: list[str, str] = ["red", "cyan"],
    ):
        """
        A class for creating anaglyph plots, that are viewed with red-cyan
        "3D glasses". Note that anaglyph plots need to be 2D. Also, any color
        arguments to plotting methods will be ignored and replaced.

        Parameters
        ----------
        fig : matplotlib.figure.Figure, optional
            The figure object to plot on.
        ax : matplotlib.axes.Axes, optional
            The axes object to plot on.
        eye_balance : float
            The eye balance parameter, from -1 to 1. A value of -1 means the
            left plot will have accurate x-axis labels, and a value of 1 means
            the right plot will. For any other value, both plots will have
            inaccurate x-axis labels.
            The x-axis will be colored to indicate which data is correct.
        d : float
            Distance from the focal plane to the viewer (in millimeters).
        ipd : float
            Interpupillary distance (in millimeters). Default is 65. Note that
            for anaglyphs, there is no cross-eyed viewing and coloring is
            determined by the `colors` argument. So any negative ipds passed in
            will be changed to their absolute value.
        zscale : Optional[float]
            Scaling factor for the visual depth of the z-data (in x-axis units).
            If None (default), then will be set to 1/4 the x-axis range.
        zlim : Optional[tuple[float, float]]
            Z-axis limits. If None (default), then z range will be autoscaled to
            (min(z), max(z)).
        zzero : Optional[float]
            The z-coordinate of the focal plane. Set to min(z) to have all the
            data float above the page, or set to max(z) to have all the data
            sink into the page. If None (default), will be set to the midpoint
            of the z range.
        colors : list[str, str]
            Colors for the left and right axes. Default is ['red', 'cyan'].
            The color ordering refers to the left and right glasses lens colors.
            Because that color prevents that eye from seeing that color data,
            the eye will see the opposite color data. Eg. ['red', 'cyan'] means
            the left eye has a red lens and sees cyan, and the right eye has a
            cyan lens and sees red.
        """
        ipd = abs(ipd)
        super().__init__(
            eye_balance=eye_balance,
            d=d,
            ipd=ipd,
            zscale=zscale,
            zlim=zlim,
            zzero=zzero,
            is_3d=False,
        )

        if fig is None and ax is None:
            self.fig, self.ax = plt.subplots()
        elif ax is None:
            self.fig = fig
            self.ax = fig.add_subplot(111)
        else:
            self.fig = ax.figure
            self.ax = ax

        self.known_methods = [
            "plot",
            "step",
            "scatter",
            "stem",
            "bar",
            "fill",
            "fill_between",
            "errorbar",
            "text",
        ]
        self.colors = colors
        self.alpha = 0.5

    def __getattr__(self, name: str):
        """
        Intercept plotting functions. If the method has 'x' and 'y' as
        arguments, and either there is a third argument or 'z' is a keyword
        argument, then the z data will be used to offset the x data for the two
        colors and create the stereoscopic effect.

        Parameters
        ----------
        name : str
            The name of the attribute.
        """

        def method(*args: Any, **kwargs: dict[str, Any]) -> Union[tuple[Any, Any], Any]:
            """
            The method that will be called on the left and right axes. If the
            method has 'x' and 'y' as arguments, and either there is a third
            argument or 'z' is a keyword argument, then the z data will be used
            to offset the x data for the left and right axes and create the
            stereoscopic effect.

            Parameters
            ----------
            *args : Any
                The positional arguments passed to the method.
            **kwargs : dict[str, Any]
                The keyword arguments passed to the method.

            Returns
            -------
            result : Union[tuple[Any, Any], Any]
                The result of the method call. If the method does not plot x-y
                data, then the result will be a single object. If the method
                does plot x-y data, then the result will be a tuple of two
                objects, one for each subplot.
            """
            ax_method = getattr(self.ax, name, None)
            args_original = args
            x, y, z, args, kwargs = process_args(ax_method, self.known_methods, args, kwargs)

            if all(var is not None for var in [ax_method, x, y, z]):
                res_left, res_right = self.plot2d(self.ax, self.ax, name, x, y, z, args, kwargs)
                result = (res_left, res_right)
                # Set the xlabel color to the right color
                self.set_axlabel_colors()
            else:
                # For methods that don't plot x-y data
                result = getattr(self.ax, name)(*args_original, **kwargs)
            return result

        return method

    def imshow_stereo(
        self,
        images: list[np.ndarray],
        method: Optional[str] = None,
        cmap: Optional[str] = None,
        crop: bool = False,
        *args: Any,
        **kwargs: dict[str, Any],
    ):
        """
        From existing stereo image data, combine into an anaglyph. Any further
        args or kwargs will be passed on to the `imshow()` function.

        A list of two or more images may be passed; the first (left-eye view)
        and the last (right-eye view) are used to form the anaglyph.

        The data can be of shape (M, N) or (M, N, 3) or (M, N, 4). If the data
        is (M, N), then it will be converted to (M, N, 3) by stacking the data
        three times, and so the baseline image will appear in grayscale.
        If the data is (M, N, 4), then the fourth alpha channel will be passed
        through to the anaglyph.

        The methods are based on the following paper:
        Sanders, William R., and David F. McAllister. "Producing anaglyphs from
          synthetic images." Stereoscopic displays and virtual reality systems
          X. Vol. 5006. SPIE, 2003.

        Parameters
        ----------
        images : list[numpy.ndarray]
            The stereo images, ordered from the left-eye view to the right-eye
            view. The first and last are combined into the anaglyph.
        method : str
            The method used to create the anaglyph. Options:
            'dubois', default when cmap is None, will always be red-cyan
            'photoshop'
            'photoshop2', default when cmap is not None
        cmap : str
            The colormap to use, default None. See above for default behavior.
            Only recommended for 1-D shape (M, N) scalar data rather than color
            images.
        crop : bool
            Whether to crop the image to the minimum size of the two images,
            keeping the images centered. Default is False.
        """
        if len(images) < 2:
            raise ValueError("imshow_stereo requires at least two images")
        data_left, data_right = images[0], images[-1]

        # Check that the method is valid
        if method is None:
            if cmap is None:
                method = "dubois"
            else:
                method = "photoshop2"
        if method.lower() not in ["photoshop", "photoshop2", "dubois"]:
            raise ValueError(
                f"method={method} must be one of 'photoshop', ", "'photoshop2', 'dubois'"
            )
        method = method.lower()

        data_left, data_right = sanitize_data_left_right(data_left, data_right, crop)

        # Map grayscale to rgb
        if len(data_left.shape) == 2:
            data_left = np.stack((data_left,) * 3, axis=-1)
            data_right = np.stack((data_right,) * 3, axis=-1)
        elif cmap is not None:
            warnings.warn("cmap is not recommended for color images")

        # Luma conversion
        luma_map = np.array([0.2126, 0.7152, 0.0722])

        # Color tuples (only used for 'photoshop', 'photoshop2' methods)
        color_left = mpl.colors.to_rgb(self.colors[0])
        color_right = mpl.colors.to_rgb(self.colors[1])

        if method == "dubois":
            M_left = np.array(
                [
                    [0.4561, 0.500484, 0.176381],
                    [-0.0400822, -0.0378246, -0.0157589],
                    [-0.0152161, -0.0205971, -0.00546856],
                ]
            )
            M_right = np.array(
                [
                    [-0.0434706, -0.0879388, -0.00155529],
                    [0.378476, 0.73364, -0.0184503],
                    [-0.0721527, -0.112961, 1.2264],
                ]
            )

        elif method == "photoshop":
            M_left = np.diag(color_left)
            M_right = np.diag(color_right)

        elif method == "photoshop2":
            # Known as "modified photoshop" in the paper
            M_luma = np.array([luma_map, [0, 1, 0], [0, 0, 1]])
            M_left = np.diag(color_left) @ M_luma
            M_right = np.diag(color_right)

        # Modify matrices to pass through alpha channel if it exists
        if data_left.shape[-1] == 4:
            M_left = np.block([[M_left, np.zeros((3, 1))], [np.zeros((1, 3)), 1]])
            M_right = np.block([[M_right, np.zeros((3, 1))], [np.zeros((1, 3)), 1]])
            luma_map = np.append(luma_map, 0)

        if cmap is None:
            # Matrix multiply every colored pixel
            data_colored = np.einsum("ij,klj->kli", M_left, data_left) + np.einsum(
                "ij,klj->kli", M_right, data_right
            )

        else:  # camp is not None
            cmap_left = copy.deepcopy(plt.get_cmap(cmap))
            cmap_right = copy.deepcopy(cmap_left)
            if isinstance(cmap_left, mpl.colors.ListedColormap):
                colors_left = []
                colors_right = []
                for color in cmap_left.colors:
                    colors_left.append(np.clip(M_left @ color, 0, 1).tolist())
                    colors_right.append(np.clip(M_right @ color, 0, 1).tolist())
                cmap_left.colors = colors_left
                cmap_right.colors = colors_right

            elif isinstance(cmap_left, mpl.colors.LinearSegmentedColormap):
                colors_left = []
                colors_right = []
                for i in range(len(cmap_left._segmentdata["red"])):
                    color = [
                        cmap_left._segmentdata[channel][i][1]
                        for channel in ["red", "green", "blue"]
                    ]
                    colors_left.append(np.clip(M_left @ color, 0, 1).tolist())
                    colors_right.append(np.clip(M_right @ color, 0, 1).tolist())
                for i, channel in enumerate(["red", "green", "blue"]):
                    entries_left = []
                    entries_right = []
                    for j in range(len(cmap_left._segmentdata[channel])):
                        entries_left.append(
                            (
                                cmap_left._segmentdata[channel][j][0],
                                colors_left[j][i],
                                colors_left[j][i],
                            )
                        )
                        entries_right.append(
                            (
                                cmap_right._segmentdata[channel][j][0],
                                colors_right[j][i],
                                colors_right[j][i],
                            )
                        )
                    cmap_left._segmentdata[channel] = entries_left
                    cmap_right._segmentdata[channel] = entries_right

            data_left_flattened = np.clip(data_left @ luma_map, 0, 1)
            data_right_flattened = np.clip(data_right @ luma_map, 0, 1)
            data_colored = cmap_left(data_left_flattened) + cmap_right(data_right_flattened)

        data_colored = np.clip(data_colored, 0, 1)
        self.ax.imshow(data_colored, *args, **kwargs)

    def set_axlabel_colors(self):
        """
        If either the first or second plotted data correctly matches the axis
        labels, then color the labels to indicate which data is correct.
        """
        if self.eye_balance == -1:
            for label in self.ax.get_xticklabels():
                label.set_color(self.colors[1])
        elif self.eye_balance == 1:
            for label in self.ax.get_xticklabels():
                label.set_color(self.colors[0])

    def _plot_area_fig(self) -> Figure:
        """Return a copy of the figure cropped to the anaglyph plot area."""
        return crop_fig_to_plot_area(self._copy_fig(), [0])
