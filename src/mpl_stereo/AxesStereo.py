import matplotlib.pyplot as plt
import numpy as np
import inspect

from abc import ABC
from typing import Optional, Any
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
        - focal_plane : float, optional
            Location of the focal plane, from -1 to 1. A value of -1 means all
            data will float above the plane. A value of 1 means all data will
            float below the plane. A value of 0 puts the focal plane on the
            page.
        - z_scale : float, optional
            Scaling factor for the z-data (in millimeters). Default is 2.
        - d : float, optional
            Distance from the focal plane to the viewer (in millimeters).
        - ipd : float, optional
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        - is_3d : bool, optional
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
        - fig : matplotlib.figure.Figure
            The figure object to which these axes belong.
        - focal_plane : float, optional
            Location of the focal plane, from -1 to 1. A value of -1 means all
            data will float above the plane. A value of 1 means all data will
            float below the plane. A value of 0 puts the focal plane on the
            page.
        - z_scale : float, optional
            Scaling factor for the z-data (in millimeters). Default is 2.
        - d : float, optional
            Distance from the focal plane to the viewer (in millimeters).
        - ipd : float, optional
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        - is_3d : bool, optional
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
        - fig : matplotlib.figure.Figure
            The figure object to which these axes belong.
        - focal_plane : float, optional
            Location of the focal plane, from -1 to 1.
            A value of -1 means all data will float above the plane,
            and only the left axis labels are accurate. (The right labels will
            have transparancy applied)
            A value of 1 means all data will float below the plane,
            and only the right axis labels are accurate. (The left labels will
            have transparancy applied)
            A value of 0 puts the focal plane on the page and neither axes'
            labels are accurate. (Both will have transparancy applied)
        - z_scale : float, optional
            Scaling factor for the z-data (in millimeters). Default is 2.
        - d : float, optional
            Distance from the focal plane to the viewer (in millimeters).
        - ipd : float, optional
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        """
        super().__init__(fig=fig, focal_plane=focal_plane, z_scale=z_scale, d=d, ipd=ipd,
                         is_3d=False)
        self.known_methods = ['plot', 'scatter', 'stem', 'bar']

        # Turn on all the axis labels then make semitransparent not accurate at the edges
        self.ax_right.yaxis.set_tick_params(labelleft=True)
        self.set_axlabel_alphas(0.5)

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
                offset_left, offset_right = calc_2d_offsets(self.focal_plane, z, self.z_scale,
                                                            self.d, self.ipd)

                getattr(self.ax_left, name)(x + offset_left, y, *args, **kwargs)
                return getattr(self.ax_right, name)(x - offset_right, y, *args, **kwargs)
            else:
                # For methods that don't plot x-y data
                getattr(self.ax_left, name)(*args_original, **kwargs)
                return getattr(self.ax_right, name)(*args_original, **kwargs)

        return method

    def set_axlabel_alphas(self, alpha: float):
        if self.focal_plane != -1:
            for label in self.ax_left.get_xticklabels():
                label.set_alpha(alpha)
            for label in self.ax_left.get_yticklabels():
                label.set_alpha(alpha)
        elif self.focal_plane != 1:
            for label in self.ax_right.get_xticklabels():
                label.set_alpha(alpha)
            for label in self.ax_right.get_yticklabels():
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
        - fig : matplotlib.figure.Figure
            The figure object to which these axes belong.
        - focal_plane : float, optional
            Location of the focal plane, from -1 to 1.
            A value of -1 means all data will float above the plane.
            A value of 1 means all data will float below the plane.
            A value of 0 puts the focal plane on the page.
        - z_scale : float, optional
            Scaling factor for the z-data (in millimeters). Default is 2.
        - d : float, optional
            Distance from the focal plane to the viewer (in millimeters).
        - ipd : float, optional
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        """
        super().__init__(fig=fig, focal_plane=focal_plane, z_scale=z_scale, d=d, ipd=ipd,
                         is_3d=True)
        self.known_methods = ['plot', 'scatter', 'stem', 'voxels', 'plot_wireframe',
                              'plot_surface', 'plot_trisurf', 'contour', 'contourf']

        self.ax_left.sharez(self.ax_right)

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
                ang = 90 - np.rad2deg(np.arctan(self.d / self.ipd))
                offset = ang * self.z_scale / self.ax_left._dist
                offset_left = (self.focal_plane + 1)/2 * offset
                offset_right = (1 - self.focal_plane)/2 * offset

                azim_init = self.ax_left.azim
                elev_init = self.ax_left.elev
                roll_init = self.ax_left.roll
                self.ax_left.view_init(elev=elev_init,
                                       azim=azim_init - offset_left,
                                       roll=roll_init)
                self.ax_right.view_init(elev=elev_init,
                                        azim=azim_init + offset_right,
                                        roll=roll_init)

                getattr(self.ax_left, name)(*args, **kwargs)
                return getattr(self.ax_right, name)(*args, **kwargs)

            else:
                # For methods that do not involve 'x' and 'y'
                getattr(self.ax_left, name)(*args, **kwargs)
                return getattr(self.ax_right, name)(*args, **kwargs)

        return method
