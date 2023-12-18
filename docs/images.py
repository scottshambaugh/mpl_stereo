import matplotlib
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_stereo import AxesStereo2D, AxesStereo3D, AxesAnaglyph, calc_2d_offsets

N_STEPS = 10

def generate_trefoil(i=0):
    dt = 2*np.pi*i/100/N_STEPS
    t = np.linspace(0, 2*np.pi, 100) + dt
    x = np.cos(2*t) * (3 + np.cos(3*t))
    y = np.sin(2*t) * (3 + np.cos(3*t))
    z = np.sin(3*t)
    return x, y, z

def plot_2d_trefoil(savedir):
    x, y, z = generate_trefoil()
    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.fig.set_size_inches(6.0, 3)
    plt.savefig(savedir / 'trefoil_2d.png', dpi=100)

def plot_3d_trefoil(savedir):
    x, y, z = generate_trefoil()
    axstereo = AxesStereo3D()
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.fig.set_size_inches(6.0, 3)
    plt.savefig(savedir / 'trefoil_3d.png', dpi=100)

def plot_anaglyph_trefoil(savedir):
    x, y, z = generate_trefoil()
    axstereo = AxesAnaglyph()
    axstereo.plot(x, y, z)
    axstereo.scatter(x, y, z, s=10)
    axstereo.fig.set_size_inches(3.0, 3)
    plt.savefig(savedir / 'trefoil_anaglyph.png', dpi=100)

def animate_2d_trefoil(savedir):
    x, y, z = generate_trefoil(0)
    cmap = matplotlib.colormaps['viridis']
    axstereo = AxesStereo2D()
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    scatter = axstereo.scatter(x, y, z, c=z, cmap=cmap, s=10)
    axstereo.fig.set_size_inches(6.0, 3)

    def animate(frame):
        x, y, z = generate_trefoil(frame)
        x, y, z = sort_by_z(x, y, z)
        offset_left, offset_right = calc_2d_offsets(axstereo.focal_plane, z, axstereo.z_scale,
                                                    axstereo.d, axstereo.ipd)
        scatter[0].set_offsets(np.stack([x + offset_left, y]).T)
        scatter[1].set_offsets(np.stack([x - offset_right, y]).T)
        for scat in scatter:
            scat.set_edgecolor(cmap(z))
            scat.set_facecolor(cmap(z))
        return scatter

    ani = animation.FuncAnimation(axstereo.fig, animate, frames=np.arange(1, N_STEPS), interval=20, repeat=False)
    ani.save(savedir / "trefoil_animation.gif", fps=10, dpi=100)


def main():
    currdir = Path(__file__).parent.resolve()
    plot_2d_trefoil(currdir)
    plot_3d_trefoil(currdir)
    plot_anaglyph_trefoil(currdir)
    animate_2d_trefoil(currdir)

if __name__ == '__main__':
    main()
