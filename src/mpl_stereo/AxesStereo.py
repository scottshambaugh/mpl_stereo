import matplotlib.pyplot as plt
import numpy as np
import inspect

from abc import ABC
from typing import Optional, Any
from types import MethodType
from matplotlib import _api
from matplotlib.figure import Figure

## Functions
def sort_by_z(x: np.ndarray, y: np.ndarray, z: np.ndarray, kwargs: dict[str, Any]):
    """
    Sort the data by z to not occlude improperly
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


def process_args(ax_method: Any, known_methods: list[str], args: Any, kwargs: dict[str, Any]):
    """
    Process the arguments to a method call to determine if the method is
    plotting x-y data and if there is a z argument or keyword argument.
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


def calc_2d_offsets(focal_plane: float, z: np.ndarray, z_scale: float, d: float, ipd: float):
    """
    Calculate the x-offsets for a 2D plot
    """
    z_scaled = z / np.ptp(z) * z_scale
    offset = ipd * z_scaled / (d + z_scaled)
    offset_left = (focal_plane + 1)/2 * offset
    offset_right = (1 - focal_plane)/2 * offset
    return offset_left, offset_right


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


## Classes
class AxesStereoBase(ABC):
    def __init__(self,
                 focal_plane: float = -1,
                 z_scale: float = 2,
                 d: float = 350,
                 ipd: float = 65,
                 is_3d: bool = False):
        """
        Parameters
        ----------
        - focal_plane : float
            Location of the focal plane, from -1 to 1. A value of -1 means all
            data will float above the plane. A value of 1 means all data will
            float below the plane. A value of 0 puts the focal plane on the
            page.
        - z_scale : float
            Scaling factor for the z-data (in millimeters). Default is 2.
        - d : float
            Distance from the focal plane to the viewer (in millimeters).
        - ipd : float
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        - is_3d : bool
            Whether the axes are 3D. Default is False.
        """
        self.focal_plane = focal_plane
        self.z_scale = z_scale
        self.d = d
        self.ipd = ipd
        self.is_3d = is_3d
        self.known_methods: list[str] = []


class AxesStereo(AxesStereoBase):
    def __init__(self,
                 fig: Optional[Figure] = None,
                 focal_plane: float = -1,
                 z_scale: float = 2,
                 d: float = 350,
                 ipd: float = 65,
                 is_3d: bool = False):
        """
        Parameters
        ----------
        - fig : matplotlib.figure.Figure, optional
            The figure object to which these axes belong.
        - focal_plane : float
            Location of the focal plane, from -1 to 1. A value of -1 means all
            data will float above the plane. A value of 1 means all data will
            float below the plane. A value of 0 puts the focal plane on the
            page.
        - z_scale : float
            Scaling factor for the z-data (in millimeters). Default is 2.
        - d : float
            Distance from the focal plane to the viewer (in millimeters).
        - ipd : float
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        - is_3d : bool
            Whether the axes are 3D. Default is False.
        """
        super().__init__(focal_plane=focal_plane, z_scale=z_scale, d=d, ipd=ipd, is_3d=is_3d)

        # Generate two side-by-side subplots
        if fig is None:
            if not is_3d:
                fig, axs = plt.subplots(1, 2)
            else:
                fig, axs = plt.subplots(1, 2, subplot_kw={'projection': '3d'})
            self.ax_left = axs[0]
            self.ax_right = axs[1]
        else:
            if not is_3d:
                self.ax_left = fig.add_subplot(121)
                self.ax_right = fig.add_subplot(122)
            else:
                self.ax_left = fig.add_subplot(121, projection='3d')
                self.ax_right = fig.add_subplot(122, projection='3d')

        self.ax_left.sharex(self.ax_right)
        self.ax_left.sharey(self.ax_right)

        self.fig = fig
        self.axs = (self.ax_left, self.ax_right)


class AxesStereo2D(AxesStereo):
    def __init__(self,
                 fig: Optional[Figure] = None,
                 focal_plane: float = -1,
                 z_scale: float = 2,
                 d: float = 350,
                 ipd: float = 65):
        """
        A class for creating stereoscopic 2D plots.

        Parameters
        ----------
        - fig : matplotlib.figure.Figure, optional
            The figure object to which these axes belong.
        - focal_plane : float
            Location of the focal plane, from -1 to 1.
            A value of -1 means all data will float above the plane,
            and only the left axis labels are accurate. (The right labels will
            have transparancy applied)
            A value of 1 means all data will float below the plane,
            and only the right axis labels are accurate. (The left labels will
            have transparancy applied)
            A value of 0 puts the focal plane on the page and neither axes'
            labels are accurate. (Both will have transparancy applied)
        - z_scale : float
            Scaling factor for the z-data (in millimeters). Default is 2.
        - d : float
            Distance from the focal plane to the viewer (in millimeters).
        - ipd : float
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        """
        super().__init__(fig=fig, focal_plane=focal_plane, z_scale=z_scale, d=d, ipd=ipd,
                         is_3d=False)
        self.known_methods = ['plot', 'scatter', 'stem', 'bar']

        # Give the innacurate x-axis labels some transparency
        self.set_axlabel_alphas(alpha=0.5)

    def __getattr__(self, name: str):
        """
        Delegate method calls to the left and right axes if the method is not
        defined in AxesStereo. If the method has 'x' and 'y' as arguments, and
        either there is a third argument or 'z' is a keyword argument, then the
        z data will be used to offset the x data for the left and right axes
        and create the stereoscopic effect.
        """
        def method(*args, **kwargs):
            ax_method = getattr(self.ax_left, name, None)
            args_original = args
            x, y, z, args, kwargs = process_args(ax_method, self.known_methods, args, kwargs)

            if all(var is not None for var in [ax_method, x, y, z]):
                # for scatter plots, sort the data by z to not occlude improperly
                if name == 'scatter':
                    x, y, z, kwargs = sort_by_z(x, y, z, kwargs)

                # Calculate the x-offsets
                offset_left, offset_right = calc_2d_offsets(self.focal_plane, z, self.z_scale,
                                                            self.d, self.ipd)

                # Plot the data twice, once for each subplot
                res_left = getattr(self.ax_left, name)(x + offset_left, y, *args, **kwargs)
                res_right = getattr(self.ax_right, name)(x - offset_right, y, *args, **kwargs)
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
        """
        if self.focal_plane != -1:
            for label in self.ax_left.get_xticklabels():
                label.set_alpha(alpha)
        elif self.focal_plane != 1:
            for label in self.ax_right.get_xticklabels():
                label.set_alpha(alpha)


class AxesStereo3D(AxesStereo):
    def __init__(self,
                 fig: Optional[Figure] = None,
                 focal_plane: float = -1,
                 z_scale: float = 2,
                 d: float = 350,
                 ipd: float = 65):
        """
        A class for creating stereoscopic 3D plots.

        Parameters
        ----------
        - fig : matplotlib.figure.Figure, optional
            The figure object to which these axes belong.
        - focal_plane : float
            Location of the focal plane, from -1 to 1.
            A value of -1 means all data will float above the plane.
            A value of 1 means all data will float below the plane.
            A value of 0 puts the focal plane on the page.
        - z_scale : float
            Scaling factor for the z-data (in millimeters). Default is 2.
        - d : float
            Distance from the focal plane to the viewer (in millimeters).
        - ipd : float
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        """
        super().__init__(fig=fig, focal_plane=focal_plane, z_scale=z_scale, d=d, ipd=ipd,
                         is_3d=True)
        self.known_methods = ['plot', 'scatter', 'stem', 'voxels', 'plot_wireframe',
                              'plot_surface', 'plot_trisurf', 'contour', 'contourf']

        self.ax_left.sharez(self.ax_right)

        # Override the view_init method to link the views of the left and right
        # 3D axes while maintaining the correct offset
        self.ax_left.shareview(self.ax_right)
        self.ax_left.view_init = MethodType(view_init, self.ax_left)
        self.ax_right.view_init = MethodType(view_init, self.ax_right)

    def __getattr__(self, name: str):
        """
        Delegate method calls to the left and right axes if the method is not
        defined in AxesStereo. If the method has 'x' and 'y' as arguments, and
        either there is a third argument or 'z' is a keyword argument, then the
        z data will be used to offset the x data for the left and right axes
        and create the stereoscopic effect.
        """

        def method(*args, **kwargs):
            # Reflect the method to check if 'x' and 'y' are in the arguments
            ax_method = getattr(self.ax_left, name, None)
            parameters = inspect.signature(ax_method).parameters

            is_plottable = False
            keyword_groups = [['x', 'y', 'z'], ['xs', 'ys', 'zs'], ['X', 'Y', 'Z']]
            for keyword_group in keyword_groups:
                if all([keyword in parameters for keyword in keyword_group]):
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

            else:
                # For methods that do not involve 'x' and 'y'
                res_left = getattr(self.ax_left, name)(*args, **kwargs)
                res_right = getattr(self.ax_right, name)(*args, **kwargs)
            return (res_left, res_right)

        return method

    def calc_3d_offsets(self):
        """
        Calculate the angular view offsets for a 3D plot
        """
        ang = 90 - np.rad2deg(np.arctan(self.d / self.ipd))
        offset = ang * self.z_scale / self.ax_left._dist
        offset_left = (self.focal_plane + 1)/2 * offset
        offset_right = (1 - self.focal_plane)/2 * offset
        return offset_left, offset_right


class AxesAnaglyph(AxesStereoBase):
    def __init__(self,
                 fig: Optional[Figure] = None,
                 focal_plane: float = -1,
                 z_scale: float = 2,
                 d: float = 350,
                 ipd: float = 65,
                 colors: list[str] = ['red', 'cyan']):
        """
        A class for creating anaglyph plots, that are viewed with red-cyan
        "3D glasses". Note that anaglyph plots need to be 2D. Also, any color
        arguments to plotting methods will be ignored and replaced.

        Parameters
        ----------
        - fig : matplotlib.figure.Figure, optional
            The figure object to which these axes belong.
        - focal_plane : float
            Location of the focal plane, from -1 to 1. A value of -1 means all
            data will float above the plane. A value of 1 means all data will
            float below the plane. A value of 0 puts the focal plane on the
            page.
        - z_scale : float
            Scaling factor for the z-data (in millimeters). Default is 2.
        - d : float
            Distance from the focal plane to the viewer (in millimeters).
        - ipd : float
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        - colors : list[str]
            Colors for the left and right axes. Default is ['red', 'cyan'].
            The color ordering refers to the left and right glasses lens colors.
            Because that color prevents that eye from seeing that color data,
            the eye will see the opposite color data. Eg. ['red', 'cyan'] means
            the left eye has a red lens and sees cyan, and the right eye has a
            cyan lens and sees red.
        """
        super().__init__(focal_plane=focal_plane, z_scale=z_scale, d=d, ipd=ipd, is_3d=False)

        if fig is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.ax = fig.add_subplot(111)

        self.known_methods = ['plot', 'scatter', 'bar']
        self.colors = colors
        self.alpha = 0.5

    def __getattr__(self, name: str):
        """
        Intercept plotting functions. If the method has 'x' and 'y' as
        arguments, and either there is a third argument or 'z' is a keyword
        argument, then the z data will be used to offset the x data for the two
        colors and create the stereoscopic effect.
        """
        def method(*args, **kwargs):
            ax_method = getattr(self.ax, name, None)
            args_original = args
            x, y, z, args, kwargs = process_args(ax_method, self.known_methods, args, kwargs)

            if all(var is not None for var in [ax_method, x, y, z]):
                offset_left, offset_right = calc_2d_offsets(self.focal_plane, z, self.z_scale,
                                                            self.d, self.ipd)
                # Delete any color arguments
                kwargs.pop('c', None)
                kwargs.pop('color', None)
                kwargs.pop('cmap', None)
                kwargs.pop('alpha', None)

                # Set the xlabel color to the right color
                self.set_axlabel_colors()

                # Plot the data twice, once for each color
                res_left = getattr(self.ax, name)(x + offset_left, y,
                                                  color=self.colors[1], alpha=self.alpha,
                                                  *args, **kwargs)
                res_right = getattr(self.ax, name)(x - offset_right, y,
                                                   color=self.colors[0], alpha=self.alpha,
                                                   *args, **kwargs)
                result = (res_left, res_right)
            else:
                # For methods that don't plot x-y data
                result = getattr(self.ax, name)(*args_original, **kwargs)
            return result

        return method

    def set_axlabel_colors(self):
        """
        If either the first or second plotted data correctly matches the axis
        labels, then color the labels to indicate which data is correct.
        """
        if self.focal_plane == -1:
            for label in self.ax.get_xticklabels():
                label.set_color(self.colors[1])
        elif self.focal_plane == 1:
            for label in self.ax.get_xticklabels():
                label.set_color(self.colors[0])
