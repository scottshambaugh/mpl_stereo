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
from matplotlib import _api
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d.axes3d import Axes3D


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
    if 'c' in kwargs and np.array(kwargs['c']).shape == np.array(z).shape:
        c = kwargs.pop('c')
        c = c[sort_idx]
        kwargs['c'] = c
    return x, y, z, kwargs


def process_args(ax_method: Any,
                 known_methods: list[str],
                 args: Any,
                 kwargs: dict[str, Any]):
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
    parameters = inspect.signature(ax_method).parameters
    if 'x' in kwargs:
        x = kwargs.pop('x')
    elif 'x' in parameters or ax_method.__name__ in known_methods:
        x, *args = args
    if 'y' in kwargs:
        y = kwargs.pop('y')
    elif 'y' in parameters or ax_method.__name__ in known_methods:
        y, *args = args

    # Check if 'z' is in the keyword arguments or if there is a third
    # argument of the same shape as x
    if 'z' in kwargs:
        z = kwargs.pop('z')
    elif len(args) > 0 and (np.array(args[0]).shape == np.array(x).shape):
        z, *args = args
    return x, y, z, args, kwargs


def calc_2d_offsets(eye_balance: float,
                    z: np.ndarray,
                    d: float,
                    ipd: float,
                    zautoscale: bool = True,
                    zscale: Optional[float] = None,
                    zlim: Optional[tuple[float, float]] = None,
                    zzero: Optional[float] = None,
                    xlim: Optional[tuple[float, float]] = None):
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
    z_midpoint = zlim[0] + zrange/2
    if zzero is None:
        zzero = z_midpoint

    zscaled = (z + z_midpoint - zzero) / zrange * zscale
    offset = ipd * zscaled / (d + zscaled)
    offset_left = (eye_balance + 1)/2 * offset
    offset_right = (1 - eye_balance)/2 * offset
    return offset_left, offset_right, zlim, zscale


def view_init(self, elev=None, azim=None, roll=None, vertical_axis="z",
              share=False):
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
    vertical_axis = _api.check_getitem(
        dict(x=0, y=1, z=2), vertical_axis=vertical_axis
    )

    if share:
        axes = {sibling for sibling
                in self._shared_axes['view'].get_siblings(self)}
    else:
        axes = [self]

    for ax in axes:
        ax.elev = elev
        ax.azim = azim
        ax.roll = roll
        if hasattr(ax, "stereo_offset") and hasattr(self, "stereo_offset") and ax is not self:
            ax.azim += (ax.stereo_offset - self.stereo_offset)
        ax._vertical_axis = vertical_axis


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
        raise ValueError(f'Crop shape ({cropx}, {cropy}) is larger than data shape ({x}, {y})')
    startx = x//2 - (cropx//2)
    starty = y//2 - (cropy//2)
    cropped_data = data[startx:startx + cropx, starty:starty + cropy, :]
    if ndim == 2:
        cropped_data = cropped_data[:, :, 0]
    return cropped_data


def sanitize_data_left_right(data_left: np.ndarray,
                             data_right: np.ndarray,
                             crop: bool) -> tuple[np.ndarray, np.ndarray]:
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
        data_left = data_left.astype(float)/255
        data_right = data_right.astype(float)/255
    return data_left, data_right


## Classes
class AxesStereoBase(ABC):
    def __init__(self,
                 eye_balance: float = -1,
                 d: float = 350,
                 ipd: float = 65,
                 zscale: Optional[float] = None,
                 zlim: Optional[tuple[float, float]] = None,
                 zzero: Optional[float] = None,
                 is_3d: bool = False):
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

    def log_artists(self,
                    res_left: Any,
                    res_right: Any,
                    name: str,
                    args: Any,
                    kwargs: dict[str, Any]):
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
    def __init__(self,
                 fig: Optional[Figure] = None,
                 axs: Optional[Union[tuple[Axes, Axes], tuple[Axes3D, Axes3D]]] = None,
                 eye_balance: float = -1,
                 d: float = 350,
                 ipd: float = 65,
                 zscale: Optional[float] = None,
                 zlim: Optional[tuple[float, float]] = None,
                 zzero: Optional[float] = None,
                 is_3d: bool = False):
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
        super().__init__(eye_balance=eye_balance, d=d, ipd=ipd, zscale=zscale,
                         zlim=zlim, zzero=zzero, is_3d=is_3d)

        # Generate two side-by-side subplots
        if fig is None and axs is None:
            if not is_3d:
                fig, axs = plt.subplots(1, 2)
            else:
                fig, axs = plt.subplots(1, 2, subplot_kw={'projection': '3d'})
            self.ax_left = axs[0]
            self.ax_right = axs[1]
        elif axs is None:
            if not is_3d:
                self.ax_left = fig.add_subplot(121)
                self.ax_right = fig.add_subplot(122)
            else:
                self.ax_left = fig.add_subplot(121, projection='3d')
                self.ax_right = fig.add_subplot(122, projection='3d')
        else:
            fig = axs[0].figure
            self.ax_left = axs[0]
            self.ax_right = axs[1]

        self.ax_left.sharex(self.ax_right)
        self.ax_left.sharey(self.ax_right)

        self.fig = fig
        self.axs = (self.ax_left, self.ax_right)

    def wiggle(self, filepath: Union[str, Path], interval: float = 125,
               ax: Optional[Axes] = None, yaxis_off: bool = False,
               *args: Any, **kwargs: dict[str, Any]):
        """
        Save the figure as a wiggle stereogram.

        Parameters
        ----------
        filepath : str | pathlib.Path
            The filepath to save the figure to.
        interval : float
            The interval between frames in milliseconds, default 125.
        ax : matplotlib.axes.Axes, optional
            The target axes to plot the wiggle stereogram on. If None (default),
            then will plot on the axes of a new Figure.
        yaxis_off : bool
            Whether to hide the y-axis. Default is False.
        *args : Any
            Additional arguments passed to `animation.save`.
        **kwargs : dict[str, Any]
            Additional keyword arguments passed to `animation.save`.
        """
        filepath = Path(filepath)
        if ax is None:
            fig, ax_target = plt.subplots()
        else:
            ax_target = ax
            fig = ax_target.figure
        pos = ax_target.get_position()

        # Duplicate the figure
        buf = io.BytesIO()
        pickle.dump(self.fig, buf)
        buf.seek(0)
        fig_buffer = pickle.load(buf)

        fig.delaxes(ax_target)

        # Set up the axes to fill the figure
        axs = fig_buffer.axes[0:2]
        for ax in axs:
            ax.figure = fig
            fig.axes.append(ax)
            fig.add_axes(ax)
            ax.set_position(pos)
            ax.set_visible(False)  # Hide initially
            if yaxis_off:
                ax.yaxis.set_visible(False)

            # workaround for https://github.com/matplotlib/matplotlib/issues/28448
            artists = ax.get_children()
            for artist in artists:
                if isinstance(artist, mpl.image.AxesImage):
                    array = artist.get_array()
                    artist.set_array(array.data.astype('float64'))

        def update(frame):
            axs[frame].set_visible(True)

            # Don't hide the underlying left axis for 2D plots, since the right
            # one has the y axis hidden
            if self.is_3d and frame == 1:
                axs[0].set_visible(False)
            return ax,

        ani = FuncAnimation(fig, update, frames=2, interval=interval)

        ani.save(filepath, *args, **kwargs)


class AxesStereo2DBase(ABC):
    def set_zlim(self,
                 zlim: tuple[float, float],
                 zscale: Optional[float] = None,
                 zautoscale: bool = False,
                 redraw: Optional[bool] = None):
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
        if (redraw is None
            and (self.zlim != zlim
                 or (self.zscale != zscale and zscale is not None))):
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
            kwargs.pop('zlim', None)
        self.redraw()

    def _calc_bounding_zlim(self) -> tuple[float, float]:
        """
        Calculate the z limits that will bound all the z data for all artists.
        """
        zlim = (np.inf, -np.inf)
        for _, args, kwargs in self.artist_args:
            # Check if 'z' is in the keyword arguments or if there is a third
            # argument of the same shape as x
            if 'z' in kwargs:
                z = kwargs['z']
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

    def plot2d(self,
               ax_left: Axes,
               ax_right: Axes,
               name: str,
               x: np.ndarray,
               y: np.ndarray,
               z: np.ndarray,
               args: Any,
               kwargs: dict[str, Any]) -> tuple[Any, Any]:
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
        kwargs_original['x'] = x
        kwargs_original['y'] = y
        kwargs_original['z'] = z
        if name == 'scatter':
            x, y, z, kwargs = sort_by_z(x, y, z, kwargs)

        # Extract the zzero and zscale keyword arguments if they exist
        zzero = kwargs.pop('zzero', None)
        if zzero is None and self.zzero is not None:
            zzero = self.zzero
        zscale = kwargs.pop('zscale', None)
        if zscale is None and self.zscale is not None:
            zscale = self.zscale

        # Extract the zlim keyword argument if it exists and update limits
        zlim = kwargs.pop('zlim', None)
        if zlim is not None:
            self.set_zlim(zlim, zscale, redraw=False)

        # Calculate the x-offsets
        offset_left, offset_right, zlim, zscale  = calc_2d_offsets(self.eye_balance,
                                                                   z,
                                                                   self.d,
                                                                   self.ipd,
                                                                   self.zautoscale,
                                                                   zscale=zscale,
                                                                   zlim=self.zlim,
                                                                   zzero=zzero,
                                                                   xlim=ax_left.get_xlim())
        self.zscale = zscale
        self.zlim = zlim
        if len(self.artist_args) > 0 and not self.is_redrawing:
            self.redraw()

        if isinstance(self, AxesStereo2D):
            # Plot the data twice, once for each subplot (2D case)
            res_left = getattr(ax_left, name)(x + offset_left, y, *args, **kwargs)
            res_right = getattr(ax_right, name)(x - offset_right, y, *args, **kwargs)
        elif isinstance(self, AxesAnaglyph):
            # Clear all color arguments (anaglyph case)
            kwargs.pop('c', None)
            kwargs.pop('color', None)
            kwargs.pop('cmap', None)
            kwargs.pop('alpha', None)

            # Plot the data twice, once for each color
            res_left = getattr(ax_left, name)(x + offset_left, y,
                                              color=self.colors[1], alpha=self.alpha,
                                              *args, **kwargs)
            res_right = getattr(ax_right, name)(x - offset_right, y,
                                                color=self.colors[0], alpha=self.alpha,
                                                *args, **kwargs)

        # Keep track of the artists
        self.log_artists(res_left, res_right, name, args, kwargs_original)

        return res_left, res_right


class AxesStereo2D(AxesStereoSideBySide, AxesStereo2DBase):
    def __init__(self,
                 fig: Optional[Figure] = None,
                 axs: Optional[tuple[Axes, Axes]] = None,
                 eye_balance: float = -1,
                 d: float = 350,
                 ipd: float = 65,
                 zscale: Optional[float] = None,
                 zlim: Optional[tuple[float, float]] = None,
                 zzero: Optional[float] = None):
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
        super().__init__(fig=fig, axs=axs, eye_balance=eye_balance, d=d, ipd=ipd,
                         zlim=zlim, zscale=zscale, zzero=zzero, is_3d=False)
        self.known_methods = ['plot', 'scatter', 'stem', 'bar', 'text']

        # Minimize whitespace between plots
        self.fig.subplots_adjust(wspace=0.01)
        self.ax_right.tick_params(axis='y', length=0, labelcolor=(0, 0, 0, 0))

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
                res_left, res_right = self.plot2d(self.ax_left, self.ax_right,
                                                  name, x, y, z, args, kwargs)
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


class AxesStereo3D(AxesStereoSideBySide):
    def __init__(self,
                 fig: Optional[Figure] = None,
                 axs: Optional[tuple[Axes3D, Axes3D]] = None,
                 eye_balance: float = -1,
                 d: float = 350,
                 ipd: float = 65):
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
        super().__init__(fig=fig, axs=axs, eye_balance=eye_balance, d=d, ipd=ipd,
                         zscale=None, zlim=None, zzero=None, is_3d=True)
        self.known_methods = ['plot', 'plot3D', 'scatter', 'stem', 'voxels',
                              'plot_wireframe', 'plot_surface', 'plot_trisurf',
                              'contour', 'contourf']

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
            keyword_groups = [['x', 'y', 'z'], ['xs', 'ys', 'zs'], ['X', 'Y', 'Z']]
            for keyword_group in keyword_groups:
                if all([keyword in parameters for keyword in keyword_group]):
                    is_plottable = True
            if name in ('plot', 'plot3D'):
                # 'zs' is implicit and not in the parameters
                is_plottable = True

            if (ax_method and is_plottable):
                offset_left, offset_right = self.calc_3d_offsets()

                # Set the views for both subplots
                self.ax_left.view_init(azim=self.ax_left.azim - offset_left)
                self.ax_right.view_init(azim=self.ax_right.azim + offset_right)
                self.ax_left.stereo_offset = -offset_left
                self.ax_right.stereo_offset = offset_right

                # Plot the data twice, once for each subplot
                res_left = getattr(self.ax_left, name)(*args, **kwargs)
                res_right = getattr(self.ax_right, name)(*args, **kwargs)

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
        offset_left = (self.eye_balance + 1)/2 * offset
        offset_right = (1 - self.eye_balance)/2 * offset
        return offset_left, offset_right


class AxesAnaglyph(AxesStereoBase, AxesStereo2DBase):
    def __init__(self,
                 fig: Optional[Figure] = None,
                 ax: Optional[Axes] = None,
                 eye_balance: float = -1,
                 d: float = 350,
                 ipd: float = 65,
                 zscale: Optional[float] = None,
                 zlim: Optional[tuple[float, float]] = None,
                 zzero: Optional[float] = None,
                 colors: list[str, str] = ['red', 'cyan']):
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
        super().__init__(eye_balance=eye_balance, d=d, ipd=ipd,
                         zscale=zscale, zlim=zlim, zzero=zzero, is_3d=False)

        if fig is None and ax is None:
            self.fig, self.ax = plt.subplots()
        elif ax is None:
            self.fig = fig
            self.ax = fig.add_subplot(111)
        else:
            self.fig = ax.figure
            self.ax = ax

        self.known_methods = ['plot', 'scatter', 'bar', 'text']
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
                res_left, res_right = self.plot2d(self.ax, self.ax,
                                                  name, x, y, z, args, kwargs)
                result = (res_left, res_right)
                # Set the xlabel color to the right color
                self.set_axlabel_colors()
            else:
                # For methods that don't plot x-y data
                result = getattr(self.ax, name)(*args_original, **kwargs)
            return result

        return method

    def imshow_stereo(self, data_left: np.ndarray, data_right: np.ndarray,
                      method: Optional[str] = None, cmap: Optional[str] = None,
                      crop: bool = False, *args: Any, **kwargs: dict[str, Any]):
        """
        From existing stereo image data, combine into an anaglyph. Any further
        args or kwargs will be passed on to the `imshow()` function.

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
        data_left : numpy.ndarray
            The data from the left image.
        data_right : numpy.ndarray
            The data from the right image.
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
        # Check that the method is valid
        if method is None:
            if cmap is None:
                method = 'dubois'
            else:
                method = 'photoshop2'
        if method.lower() not in ['photoshop', 'photoshop2', 'dubois']:
            raise ValueError(f"method={method} must be one of 'photoshop', ",
                             "'photoshop2', 'dubois'")
        method = method.lower()

        data_left, data_right = sanitize_data_left_right(data_left, data_right, crop)

        # Map grayscale to rgb
        if len(data_left.shape) == 2:
            data_left = np.stack((data_left,)*3, axis=-1)
            data_right = np.stack((data_right,)*3, axis=-1)
        elif cmap is not None:
            warnings.warn("cmap is not recommended for color images")

        # Luma conversion
        luma_map = np.array([0.2126, 0.7152, 0.0722])

        # Color tuples (only used for 'photoshop', 'photoshop2' methods)
        color_left = mpl.colors.to_rgb(self.colors[0])
        color_right = mpl.colors.to_rgb(self.colors[1])

        if method == 'dubois':
            M_left = np.array([[ 0.4561,  0.500484,  0.176381],
                               [-0.0400822, -0.0378246, -0.0157589],
                               [-0.0152161, -0.0205971, -0.00546856]])
            M_right = np.array([[ -0.0434706, -0.0879388, -0.00155529],
                                [0.378476,  0.73364, -0.0184503],
                                [-0.0721527, -0.112961,  1.2264]])

        elif method == 'photoshop':
            M_left = np.diag(color_left)
            M_right = np.diag(color_right)

        elif method == 'photoshop2':
            # Known as "modified photoshop" in the paper
            M_luma = np.array([luma_map,
                               [0, 1, 0],
                               [0, 0, 1]])
            M_left = np.diag(color_left) @ M_luma
            M_right = np.diag(color_right)

        # Modify matrices to pass through alpha channel if it exists
        if data_left.shape[-1] == 4:
            M_left = np.block([[M_left,           np.zeros((3, 1))],
                               [np.zeros((1, 3)), 1]])
            M_right = np.block([[M_right,         np.zeros((3, 1))],
                                [np.zeros((1, 3)), 1]])
            luma_map = np.append(luma_map, 0)

        if cmap is None:
            # Matrix multiply every colored pixel
            data_colored = (np.einsum('ij,klj->kli', M_left, data_left) +
                            np.einsum('ij,klj->kli', M_right, data_right))

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
                for i in range(len(cmap_left._segmentdata['red'])):
                    color = [cmap_left._segmentdata[channel][i][1]
                                for channel in ['red', 'green', 'blue']]
                    colors_left.append(np.clip(M_left @ color, 0, 1).tolist())
                    colors_right.append(np.clip(M_right @ color, 0, 1).tolist())
                for i, channel in enumerate(['red', 'green', 'blue']):
                    entries_left = []
                    entries_right = []
                    for j in range(len(cmap_left._segmentdata[channel])):
                        entries_left.append((cmap_left._segmentdata[channel][j][0],
                                             colors_left[j][i],
                                             colors_left[j][i]))
                        entries_right.append((cmap_right._segmentdata[channel][j][0],
                                              colors_right[j][i],
                                              colors_right[j][i]))
                    cmap_left._segmentdata[channel] = entries_left
                    cmap_right._segmentdata[channel] = entries_right

            data_left_flattened = np.clip(data_left @ luma_map, 0, 1)
            data_right_flattened = np.clip(data_right @ luma_map, 0, 1)
            data_colored = (cmap_left(data_left_flattened) + cmap_right(data_right_flattened))

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
