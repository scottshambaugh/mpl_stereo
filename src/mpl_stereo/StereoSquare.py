import matplotlib.pyplot as plt
import numpy as np

from abc import ABC
from typing import Optional, Union, Any
from pathlib import Path

from mpl_stereo.AxesStereo import (AxesStereo2D, AxesStereo3D, AxesStereoSideBySide,
                                   AxesAnaglyph, sanitize_data_left_right)

class StereoSquareBase(ABC):
    """
    Base class for StereoSquare2D and StereoSquare3D.
    """
    def __init__(self):
        self.axs = None
        self.fig = None

        self.axesstereo: AxesStereoSideBySide = None
        self.axesanaglyph: AxesAnaglyph = None

    def __getattr__(self, name: str) -> Any:
        """
        Delegate method calls to the different axes.

        Parameters
        ----------
        name : str
            The name of the attribute.
        """
        def method(*args: Any, **kwargs: dict[str, Any]) -> tuple[Any, Any, Any]:
            """
            The method that will be called on the stereo and anaglyph axes.

            Parameters
            ----------
            *args : Any
                The positional arguments passed to the method.
            **kwargs : dict[str, Any]
                The keyword arguments passed to the method.
            """
            res_left, res_right = getattr(self.axesstereo, name)(*args, **kwargs)
            res_anaglyph = None
            if self.axesanaglyph is not None:
                res_anaglyph = getattr(self.axesanaglyph, name)(*args, **kwargs)
            return (res_left, res_right, res_anaglyph)

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
        data_left, data_right = sanitize_data_left_right(data_left, data_right, crop)
        res_left = self.axesstereo.ax_left.imshow(data_left, *args, **kwargs)
        res_right = self.axesstereo.ax_right.imshow(data_right, *args, **kwargs)
        if self.axesanaglyph is not None:
            res_anaglyph = self.axesanaglyph.imshow_stereo(data_left, data_right,
                                                           method=method, cmap=cmap, crop=crop,
                                                           *args, **kwargs)

        return (res_left, res_right, res_anaglyph)

    def wiggle(self, filepath: Union[str, Path], interval: float = 125,
               *args: Any, **kwargs: dict[str, Any]):
        """
        Save the figure with a wiggle stereogram in the bottom right plot.

        Parameters
        ----------
        filepath : str | pathlib.Path
            The filepath to save the figure to.
        interval : float
            The interval between frames in milliseconds, default 125.
        *args : Any
            Additional arguments passed to `animation.save`.
        **kwargs : dict[str, Any]
            Additional keyword arguments passed to `animation.save`.
        """
        self.axesstereo.wiggle(filepath=filepath, interval=interval,
                               ax=self.axs[1, 1], yaxis_off=True,
                               *args, **kwargs)


class StereoSquare2D(StereoSquareBase):
    """
    A stereo square plot for 2D data.

    A stereo square is a 2x2 grid of plots for showing all viewing methods at
    once. Along the top, the stereo plots are shown side-by-side. Along the
    bottom, the left plot is the anaglyph and the right plot is a placeholder
    for a wiggle stereogram animation. Create the figure with the wiggle
    stereogram animation by calling the `wiggle` method.

    """
    def __init__(self):
        super().__init__()
        self.fig, self.axs = plt.subplots(2, 2)

        self.axesstereo = AxesStereo2D(fig=self.fig, axs=self.axs[0, :])
        self.axesanaglyph = AxesAnaglyph(fig=self.fig, ax=self.axs[1, 0])


class StereoSquare3D(StereoSquareBase):
    """
    A stereo square plot for 3D data.

    A stereo square is a 2x2 grid of plots for showing all viewing methods at
    once. Along the top, the stereo plots are shown side-by-side. Along the
    bottom, the left plot is the anaglyph (currently unsupported for 3D plots),
    and the right plot is a placeholder for a wiggle stereogram animation.
    Create the figure with the wiggle stereogram animation by calling the
    `wiggle` method.
    """
    def __init__(self):
        super().__init__()
        self.fig, self.axs = plt.subplots(2, 2, subplot_kw={'projection': '3d'})

        self.axesstereo = AxesStereo3D(fig=self.fig, axs=self.axs[0, :])
        self.axesanaglyph = None  # No anaglyph for 3D plots
        self.fig.delaxes(self.axs[1, 0])
