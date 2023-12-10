import matplotlib.pyplot as plt
import numpy as np

class AxesStereo():
    def __init__(self, fig=None, focal_plane=-1, ipd=65):
        """
        - fig : matplotlib.figure.Figure
            The figure object to which these axes belong.
        - focal_plane : float, optional
            Location of the focal plane, from -1 to 1. A value of -1 means all data will float above the plane. A value of 1 means all data will float below the plane. A value of 0 puts the focal plane on the page.
        - ipd : float, optional
            Interpupillary distance (in millimeters). Default is 65. Negative values for cross-view.
        """
        # Generate two side-by-side subplots
        if fig is None:
            fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
            self.ax_left = axs[0]
            self.ax_right = axs[1]
        else:
            self.ax_left = fig.add_subplot(121)
            self.ax_right = fig.add_subplot(122, sharex=self.ax_left, sharey=self.ax_left)
        
        self.focal_plane = focal_plane
        self.ipd = ipd


    def __getattr__(self, name):
        """
        Delegate method calls to the left and right axes if the method is not defined in AxesStereo.
        """

        def method(*args, **kwargs):
            # Call the method on both the left and right axes
            getattr(self.left_ax, name)(*args, **kwargs)
            return getattr(self.right_ax, name)(*args, **kwargs)

        return method

    # Would have methods for all the 2D plotting functions
    def scatter(self, x, y, z, *args, **kwargs):
        x, y, z = np.array(x).ravel(), np.array(y).ravel(), np.array(z).ravel()
        if 'c' in kwargs and np.array(kwargs['c']).ravel().shape == z.shape:
            sorted_indices = np.argsort(z)[::-1]
            x, y, z = x[sorted_indices], y[sorted_indices], z[sorted_indices]
            kwargs['c'] = z  # Update the color argument with sorted z

        # this math is definitely not right, but is illustrative
        offset = z / np.ptp(z) * self.ipd
        offset_left = (self.focal_plane + 1)/2 * offset
        offset_right = (1 - self.focal_plane)/2 * offset

        self.ax_left.scatter(x - offset_left, y, *args, **kwargs)
        self.ax_right.scatter(x + offset_right, y, *args, **kwargs)
