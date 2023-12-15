import matplotlib.pyplot as plt
import numpy as np
import inspect
from abc import ABC

class AxesStereo(ABC):
    def __init__(self, fig=None, focal_plane=-1, ipd=65, is_3d=False):
        """
        - fig : matplotlib.figure.Figure
            The figure object to which these axes belong.
        - focal_plane : float, optional
            Location of the focal plane, from -1 to 1. A value of -1 means all
            data will float above the plane. A value of 1 means all data will
            float below the plane. A value of 0 puts the focal plane on the
            page.
        - ipd : float, optional
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        - is_3d : bool, optional
            Whether the axes are 3D. Default is False.
        """
        # Generate two side-by-side subplots
        if fig is None:
            if not is_3d:
                fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
            else:
                fig, axs = plt.subplots(1, 2, subplot_kw={'projection': '3d'},
                                        sharex=True, sharey=True)
            self.ax_left = axs[0]
            self.ax_right = axs[1]
        else:
            if not is_3d:
                self.ax_left = fig.add_subplot(121)
                self.ax_right = fig.add_subplot(122, sharex=self.ax_left, sharey=self.ax_left)
            else:
                self.ax_left = fig.add_subplot(121, projection='3d')
                self.ax_right = fig.add_subplot(122, projection='3d',
                                                sharex=self.ax_left, sharey=self.ax_left)
        self.focal_plane = focal_plane
        self.ipd = ipd
        self.is_3d = is_3d


class AxesStereo2D(AxesStereo):
    def __init__(self, fig=None, focal_plane=-1, ipd=65):
        """
        - fig : matplotlib.figure.Figure
            The figure object to which these axes belong.
        - focal_plane : float, optional
            Location of the focal plane, from -1 to 1. A value of -1 means all
            data will float above the plane. A value of 1 means all data will
            float below the plane. A value of 0 puts the focal plane on the
            page.
        - ipd : float, optional
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        """
        super().__init__(fig=fig, focal_plane=focal_plane, ipd=ipd, is_3d=False)

    def __getattr__(self, name):
        """
        Delegate method calls to the left and right axes if the method is not defined in AxesStereo.
        If the method has 'x' and 'y' as arguments, and either there is a third argument or 'z' is a keyword argument,
        then the z data will be used to offset the x data for the left and right axes and create the stereoscopic effect.
        """

        def method(*args, **kwargs):
            # Reflect the method to check if 'x' and 'y' are in the arguments
            ax_method = getattr(self.ax_left, name, None)
            parameters = inspect.signature(ax_method).parameters

            x = y = z = None
            if x in kwargs:
                x = kwargs.pop('x')
            elif 'x' in parameters:
                x, *args = args
            if y in kwargs:
                y = kwargs.pop('y')
            elif 'y' in parameters:
                y, *args = args

            # Check if 'z' is in the keyword arguments or if there is a third argument of the same shape as x
            if 'z' in kwargs:
                z = kwargs.pop('z')
            elif len(args) > 0 and (np.array(args[0]).shape == np.array(x).shape):
                z, *args = args

            if (ax_method and x is not None and y is not None and z is not None):
                offset = z / np.ptp(z) * self.ipd
                offset_left = (self.focal_plane + 1)/2 * offset
                offset_right = (1 - self.focal_plane)/2 * offset

                getattr(self.ax_left, name)(x - offset_left, y, *args, **kwargs)
                return getattr(self.ax_right, name)(x + offset_right, y, *args, **kwargs)
            else:
                # For methods that do not involve 'x' and 'y'
                getattr(self.ax_left, name)(*args, **kwargs)
                return getattr(self.ax_right, name)(*args, **kwargs)

        return method
    

class AxesStereo3D(AxesStereo):
    def __init__(self, fig=None, focal_plane=-1, ipd=65):
        """
        - fig : matplotlib.figure.Figure
            The figure object to which these axes belong.
        - focal_plane : float, optional
            Location of the focal plane, from -1 to 1. A value of -1 means all
            data will float above the plane. A value of 1 means all data will
            float below the plane. A value of 0 puts the focal plane on the
            page.
        - ipd : float, optional
            Interpupillary distance (in millimeters). Default is 65. Negative
            values for cross-view.
        """
        super().__init__(fig=fig, focal_plane=focal_plane, ipd=ipd, is_3d=True)


    def __getattr__(self, name):
        """
        Delegate method calls to the left and right axes if the method is not defined in AxesStereo.
        If the method has 'x' and 'y' as arguments, and either there is a third argument or 'z' is a keyword argument,
        then the z data will be used to offset the x data for the left and right axes and create the stereoscopic effect.
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
                offset = 5 # [deg]
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
