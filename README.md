![Release](https://img.shields.io/github/v/release/scottshambaugh/mpl_stereo?sort=semver)
![Builds](https://github.com/scottshambaugh/mpl_stereo/actions/workflows/builds.yml/badge.svg)
![Tests](https://github.com/scottshambaugh/mpl_stereo/actions/workflows/tests.yml/badge.svg)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mpl_stereo)

# mpl_stereo
Matplotlib add-on to make [stereograms](https://en.wikipedia.org/wiki/Stereoscopy).

Stereograms can significantly enhance the interpretability of 3D data by leveraging human binocular vision. Instead of trying to imagine how the flat projection of 3D data on a page would look in real life, stereograms give us "3D glasses" for 2D data with just our eyes.

It takes some practice to be able to view the stereoscopic effect for the first time, but the effort is well worth it!

## Usage
### Installation
```
pip install mpl_stereo
```
### Setup
```python
import numpy as np
from mpl_stereo import AxesStereo2D, AxesStereo3D

# Generate some data, here a (3,2) trefoil knot
t = np.linspace(0, 2*np.pi, 100)
x = np.cos(2*t) * (3 + np.cos(3*t))
y = np.sin(2*t) * (3 + np.cos(3*t))
z = np.sin(3*t)
```

### 2D plots
Currently, only a subset of matplotlib's 2D plots are officially supported. See the list by calling `axstereo.known_methods`.
```python
axstereo = AxesStereo2D()
axstereo.plot(x, y, z, c='k', alpha=0.2)
axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
```
<p float="left" align="center">
<img width="500" height="250" src="https://raw.githubusercontent.com/scottshambaugh/mpl_stereo/main/docs/trefoil_2d.png">
</p>
If you can see view the stereogram properly, you should see the knot weave in and out of the page!

### 3D plots
The stereoscopic effect in 3D can be made just by rotating the plot view, so all of matplotlib's 3D plot types are supported.
```python
axstereo = AxesStereo2D()
axstereo.plot(x, y, z, c='k', alpha=0.2)
axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
```
<p float="left" align="center">
<img width="500" height="250" src="https://raw.githubusercontent.com/scottshambaugh/mpl_stereo/main/docs/trefoil_3d.png">
</p>

## Viewing Stereograms

These are not [*auto*stereograms](https://en.wikipedia.org/wiki/Autostereogram), like the "Magic Eye" books that were popular in the 1990's. However, they use the same viewing technique. Below is ChatGPT's how-to guide on viewing these, but I'll try to find a better beginner-friendly resource to put here.

1) **Position the Stereogram**: Place it at arm's length and ensure it's level with your eyes.
2) **Relax Your Focus**: Look through the image, as if focusing on something distant, rather than the stereogram itself.
3) **Parallel Viewing**: Try to view the image with your eyes parallel, similar to how you would look at a distant object.
4) **Align and Overlap**: Adjust the angle and distance of the stereogram until the two images begin to overlap.
5) **Perceive the 3D Image**: As the images overlap, a 3D image should emerge. Keep your focus steady to maintain the illusion.
6) **Practice**: If initially unsuccessful, take breaks and try again. It might require some practice to get used to this method.

### Parallel vs Cross-Eyed Viewing
By default, the stereograms are set up for "parallel" viewing method as described above. For "cross-eyed" viewing, initialize with a negative `ipd` "inter-pupilary distance". An ipd of 65 millimeters is the default, so call `AxesStereo2D(ipd=-65)` as a default cross-eyed method.
