import matplotlib.pyplot as plt
import numpy as np

from abc import ABC
from typing import Optional, Union, Any
from pathlib import Path

from mpl_stereo.AxesStereo import (
    AxesStereo2D,
    AxesStereo3D,
    AxesStereoSideBySide,
    AxesAnaglyph,
    AxesAnaglyph3D,
    sanitize_images,
    ANIMATION_SUFFIXES,
)


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

        A list of two or more images may be passed. The first is shown on the
        left axes and the last on the right axes, and the first and last are
        combined into the anaglyph. All of the images are kept for use by
        `wiggle`.

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
            view.
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
            Whether to crop the images to the minimum common size, keeping them
            centered. Default is False.
        """
        if len(images) < 2:
            raise ValueError("imshow_stereo requires a list of at least two images")

        # Sanitize once over all images so the side-by-side and anaglyph views
        # share a consistent crop.
        images = sanitize_images(list(images), crop)
        res_left, res_right = self.axesstereo.imshow_stereo(images, cmap=cmap, *args, **kwargs)
        res_anaglyph = None
        if self.axesanaglyph is not None:
            res_anaglyph = self.axesanaglyph.imshow_stereo(
                [images[0], images[-1]], method=method, cmap=cmap, *args, **kwargs
            )

        return (res_left, res_right, res_anaglyph)

    def save(
        self,
        filepath: Union[str, Path],
        plot_area: bool = False,
        animate: Optional[bool] = None,
        interval: float = 125,
        frames: Optional[int] = 2,
        *args: Any,
        **kwargs: dict[str, Any],
    ):
        """
        Save the stereo square to a file. This is the single entry point for both
        a static image and an animated version with a wiggle in the bottom right
        plot. Whether to animate is chosen from the file extension (e.g. ``.gif``
        animates) or forced with the ``animate`` argument.

        Parameters
        ----------
        filepath : str | pathlib.Path
            The filepath to save to. The extension selects the format.
        plot_area : bool
            If True, strip all spines, ticks, and labels from every panel and
            remove the surrounding padding, so the 2x2 grid of plot areas fills
            the frame. Default is False.
        animate : Optional[bool]
            Force (True) or forbid (False) the wiggle animation. If None
            (default), inferred from the file extension.
        interval : float
            The interval between wiggle frames in milliseconds, default 125.
            Only used when animating.
        frames : Optional[int]
            The number of distinct viewpoints to sample across the stereo
            baseline, default 2. See `AxesStereoSideBySide.wiggle`. Only used
            when animating.
        *args : Any
            Additional positional arguments passed to the save call.
        **kwargs : dict[str, Any]
            Additional keyword arguments. For animations, passed to
            `animation.save`; for static images, passed to `savefig`.
        """
        filepath = Path(filepath)
        if animate is None:
            animate = filepath.suffix.lower() in ANIMATION_SUFFIXES

        if plot_area:
            # Strip the static cells and collapse the surrounding padding so the
            # grid of plot areas fills the frame. The wiggle cell is stripped by
            # passing plot_area through to the inner wiggle below.
            self.axesstereo.ax_left.set_axis_off()
            self.axesstereo.ax_right.set_axis_off()
            if self.axesanaglyph is not None:
                # AxesAnaglyph has a single .ax, AxesAnaglyph3D has overlaid .axs
                anaglyph_axes = getattr(self.axesanaglyph, "axs", (self.axesanaglyph.ax,))
                for ax in anaglyph_axes:
                    ax.set_axis_off()
            self.fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)

        if animate:
            self.axesstereo._run_wiggle(
                filepath,
                interval=interval,
                frames=frames,
                ax=self.axs[1, 1],
                yaxis_off=True,
                plot_area=plot_area,
                *args,
                **kwargs,
            )
        else:
            self.fig.savefig(filepath, *args, **kwargs)

    def wiggle(
        self,
        filepath: Union[str, Path],
        interval: float = 125,
        frames: Optional[int] = 2,
        plot_area: bool = False,
        *args: Any,
        **kwargs: dict[str, Any],
    ):
        """
        Save the figure with a wiggle stereogram in the bottom right plot.
        Convenience wrapper, equivalent to ``save`` with ``animate=True``. See
        `save` for the meaning of each argument.
        """
        self.save(filepath, plot_area, True, interval, frames, *args, **kwargs)


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
    bottom, the left plot is the anaglyph and the right plot is a placeholder
    for a wiggle stereogram animation. Create the figure with the wiggle
    stereogram animation by calling the `wiggle` method.
    """

    def __init__(self):
        super().__init__()
        self.fig, self.axs = plt.subplots(2, 2, subplot_kw={"projection": "3d"})

        self.axesstereo = AxesStereo3D(fig=self.fig, axs=self.axs[0, :])

        # The anaglyph overlays two 3D axes, so add a second one on the
        # bottom-left cell. Reuse that cell's subplotspec (with a distinct label)
        # so both axes share identical 3D aspect handling and line up exactly.
        overlay = self.fig.add_subplot(
            self.axs[1, 0].get_subplotspec(), projection="3d", label="anaglyph_overlay"
        )
        self.axesanaglyph = AxesAnaglyph3D(fig=self.fig, axs=(self.axs[1, 0], overlay))
