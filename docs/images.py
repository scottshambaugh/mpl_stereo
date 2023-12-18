import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_stereo import AxesStereo2D, AxesStereo3D, AxesAnaglyph

def main():
    currdir = Path(__file__).parent.resolve()

    data = dict()
    # Parametric equations for a (3,2) trefoil knot
    t = np.linspace(0, 2*np.pi, 100)
    x = np.cos(2*t) * (3 + np.cos(3*t))
    y = np.sin(2*t) * (3 + np.cos(3*t))
    z = np.sin(3*t)
    data['trefoil'] = (x, y, z)

    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.fig.set_size_inches(6.0, 3)
    plt.savefig(currdir / 'trefoil_2d.png', dpi=100)

    axstereo = AxesStereo3D()
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.fig.set_size_inches(6.0, 3)
    plt.savefig(currdir / 'trefoil_3d.png', dpi=100)

    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z)
    axstereo.scatter(x, y, z, s=10)
    axstereo.fig.set_size_inches(3.0, 3)
    plt.savefig(currdir / 'trefoil_anaglyph.png', dpi=100)


if __name__ == '__main__':
    main()
