![Release](https://img.shields.io/github/v/release/scottshambaugh/mpl_stereo?sort=semver)
![Builds](https://github.com/scottshambaugh/mpl_stereo/actions/workflows/builds.yml/badge.svg)
![Tests](https://github.com/scottshambaugh/mpl_stereo/actions/workflows/tests.yml/badge.svg)
[![codecov](https://codecov.io/gh/scottshambaugh/mpl_stereo/graph/badge.svg?token=V2ZSLFUK03)](https://codecov.io/gh/scottshambaugh/mpl_stereo)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mpl_stereo)

# mpl_stereo
Matplotlib add-on to make [stereograms](https://en.wikipedia.org/wiki/Stereoscopy) and [anaglyphs](https://en.wikipedia.org/wiki/Anaglyph_3D).

Stereographic images can significantly enhance the interpretability of 3D data by leveraging human binocular vision. Instead of looking at a flat projection on a page, stereograms give us "3D glasses" for 2D data with just our eyes.

It takes some practice to be able to view the stereoscopic effect for the first time, but the effort is well worth it!

## Usage
**Installation**
```
pip install mpl_stereo
```
**Setup**
```python
import numpy as np
from mpl_stereo import AxesStereo2D, AxesStereo3D, AxesAnaglyph
# Generate some data, we'll make a trefoil knot
t = np.linspace(0, 2*np.pi, 100)
x = np.cos(2*t) * (3 + np.cos(3*t))
y = np.sin(2*t) * (3 + np.cos(3*t))
z = np.sin(3*t)
```

### 2D Stereogram plots
Currently, only a subset of matplotlib's 2D plots are officially supported. See the list in `axstereo.known_methods`.
```python
axstereo = AxesStereo2D()
axstereo.plot(x, y, z, c='k', alpha=0.2)
axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
```
<p float="left" align="center">
<img width="500" height="250" src="https://raw.githubusercontent.com/scottshambaugh/mpl_stereo/main/docs/trefoil_2d.png">
</p>
When you are viewing the stereogram properly, you should see the knot weave in and out of the page!

*Warning*: Please note that for 2D plots the stereoscopic effect requires shifting data, so the data will *not* necessarily line up with the x-axis labels! Right now this is controlled with the `focal_plane` parameter. Calling `AxesStereo2D(focal_plane=-1)` (the default) will ensure that the left x-axis data is not shifted, whereas `AxesStereo2D(focal_plane=1)` will lock down the right x-axis data. The tick labels for x-axes where the data is not aligned will have transparency applied. So in the plot above, the right side labels being lighter gray indicates that you should not trust that x-axis data to be positioned correctly, but the left subplot with its black labeling is accurate.

### 3D Stereogram plots
The stereoscopic effect in 3D is made just by rotating the plot view, so all of matplotlib's 3D plot types are supported, and there are no concerns about data not lining up with the axis labels.
```python
axstereo = AxesStereo3D()
axstereo.plot(x, y, z, c='k', alpha=0.2)
axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
```
<p float="left" align="center">
<img width="500" height="250" src="https://raw.githubusercontent.com/scottshambaugh/mpl_stereo/main/docs/trefoil_3d.png">
</p>

### Red-Cyan Anaglyphs
Some 2D plots can also be made into anaglyphs, which are stereograms that can be viewed with red-cyan 3D glasses. While this allows for seeing the stereoscopic effect without training your eyes, it also means that the data cannot be otherwise colored.

The same warning as for the 2D stereo plots about shifting data applies here as well. If `focal_plane` is -1 or +1 such that the data for one of the colors is not shifted, then that color will be applied to the x-axis tick labels to show that they are accurate.
```python
axstereo = AxesAnaglyph()
axstereo.plot(x, y, z)
axstereo.scatter(x, y, z, s=10)
```
<p float="left" align="center">
<img width="250" height="250" src="https://raw.githubusercontent.com/scottshambaugh/mpl_stereo/main/docs/trefoil_anaglyph.png">
</p>

### Working With Plots
The figure and subplot axes can be accessed with the following:
```python
axstereo.fig
axstereo.axs  # (ax_left, ax_right), for AxesStereo2D and AxesStereo3D
axstereo.ax  # for AxesAnaglyph
```

Calling any method on `axstereo` will pass that method call onto all the subplot axes. In the 2D cases, the plotting methods which take in `x` and `y` arguments are intercepted and the additional `z` data is processed to obtain appropriate horizontal offsets.

```python
# For example instead of:
for ax in axstereo.axs:
    ax.set_xlabel('X Label')

# You can do the equivalent:
axstereo.set_xlabel('X Label')

# When applicable, the return values from the two method calls are returned as a tuple
(line_left, line_right) = axstereo.plot(x, y, z)
```

## Viewing Stereograms

These are not [*auto*stereograms](https://en.wikipedia.org/wiki/Autostereogram), like the "Magic Eye" books that were popular in the 1990's. However, they use the same viewing technique. Below is ChatGPT's how-to guide on viewing these, but I'll try to find a better beginner-friendly resource to put here.

1) **Position the Stereogram**: Place it at arm's length and ensure it's level with your eyes.
2) **Relax Your Focus**: Look through the image, as if focusing on something distant, rather than the stereogram itself.
3) **Parallel Viewing**: Try to view the image with your eyes parallel, similar to how you would look at a distant object.
4) **Align and Overlap**: Adjust the angle and distance of the stereogram until the two images begin to overlap.
5) **Perceive the 3D Image**: As the images overlap, a 3D image should emerge. Keep your focus steady to maintain the illusion.
6) **Practice**: If initially unsuccessful, take breaks and try again. It might require some practice to get used to this method.

### Parallel vs Cross-Eyed Viewing
By default, the stereograms are set up for "parallel" viewing method as described above. For "cross-eyed" viewing, initialize with a negative `ipd` parameter. An ipd (Inter-Pupilary Distance) of 65 millimeters is the default, so call `AxesStereo2D(ipd=-65)` for the default cross-eyed viewing.

### Derivation of Geometry
TODO
