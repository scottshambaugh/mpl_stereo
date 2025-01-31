import matplotlib
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_stereo import (AxesStereo2D, AxesStereo3D, AxesAnaglyph,
                        StereoSquare2D, StereoSquare3D,
                        calc_2d_offsets, sort_by_z)
from mpl_stereo.example_data import trefoil, sun_left_right, church_left_right

N_STEPS = 30
IPD = 65  # 65 for parallel, -65 for cross-eye

def _cmap_norm(z):
    return (z - np.min(z)) / np.ptp(z)

def plot_2d_trefoil(savedir=None, show=True):
    x, y, z = trefoil()
    axstereo = AxesStereo2D(ipd=IPD)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.fig.set_size_inches(6, 3)
    if savedir is not None:
        plt.savefig(savedir / 'trefoil_2d.png', dpi=100)
    if show:
        plt.show()

def wiggle_2d_trefoil(savedir=None, show=True):
    x, y, z = trefoil()
    axstereo = AxesStereo2D(ipd=IPD)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.fig.set_size_inches(3, 3)
    if savedir is not None:
        axstereo.wiggle(savedir / 'trefoil_2d_wiggle.gif', dpi=100)

def plot_3d_trefoil(savedir=None, show=True):
    x, y, z = trefoil()
    axstereo = AxesStereo3D(ipd=IPD)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.fig.set_size_inches(6, 3)
    if savedir is not None:
        plt.savefig(savedir / 'trefoil_3d.png', dpi=100)
    if show:
        plt.show()

def wiggle_3d_trefoil(savedir=None, show=True):
    x, y, z = trefoil()
    axstereo = AxesStereo3D(ipd=IPD)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)
    axstereo.fig.set_size_inches(3, 3)
    if savedir is not None:
        axstereo.wiggle(savedir / 'trefoil_3d_wiggle.gif', dpi=100)

def plot_anaglyph_trefoil(savedir=None, show=True):
    x, y, z = trefoil()
    axstereo = AxesAnaglyph(ipd=IPD)
    axstereo.plot(x, y, z)
    axstereo.scatter(x, y, z, s=10)
    axstereo.fig.set_size_inches(3.0, 3)
    if savedir is not None:
        plt.savefig(savedir / 'trefoil_anaglyph.png', dpi=100)
    if show:
        plt.show()

def plot_anaglyph_trefoil_zzero(savedir=None, show=True):
    x, y, z = trefoil()
    fig, axs = plt.subplots(1, 3)
    zzeros = (min(z), None, max(z))
    titles = ('zzero = min(z)\ndata floats above page',
              'zzero = None\ndata midpoint at page',
              'zzero = max(z)\ndata sinks below page')
    for ax, zzero, title in zip(axs, zzeros, titles):
        axstereo = AxesAnaglyph(ax=ax, ipd=IPD, zzero=zzero)
        axstereo.plot(x, y, z)
        axstereo.scatter(x, y, z, s=10)
        axstereo.set_title(title)
        axstereo.grid(True)
    fig.set_size_inches(12, 4)
    if savedir is not None:
        plt.savefig(savedir / 'trefoil_anaglyph_zzero.png', dpi=100)
    if show:
        plt.show()

def animate_2d_trefoil(savedir=None, show=True):
    x, y, z = trefoil(0, n_steps=N_STEPS)
    cmap = matplotlib.colormaps['viridis']
    axstereo = AxesStereo2D(ipd=IPD)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    scatter = axstereo.scatter(x, y, z, s=10)
    colors = cmap(_cmap_norm(z))
    for scat in scatter:
        scat.set_edgecolor(colors)
        scat.set_facecolor(colors)
    axstereo.fig.set_size_inches(6, 3)

    def animate(frame):
        x, y, z = trefoil(frame, n_steps=N_STEPS)
        x, y, z, _ = sort_by_z(x, y, z, kwargs=dict())
        zautoscale = False
        offset_left, offset_right, zlim, zscale = calc_2d_offsets(axstereo.eye_balance, z,
                                                                  axstereo.d, axstereo.ipd,
                                                                  zautoscale,
                                                                  axstereo.zscale,
                                                                  axstereo.zlim,
                                                                  axstereo.zzero)
        scatter[0].set_offsets(np.stack([x + offset_left, y]).T)
        scatter[1].set_offsets(np.stack([x - offset_right, y]).T)
        colors = cmap(_cmap_norm(z))
        for scat in scatter:
            scat.set_edgecolor(colors)
            scat.set_facecolor(colors)
        return scatter

    ani = animation.FuncAnimation(axstereo.fig, animate, frames=np.arange(N_STEPS),
                                  interval=20, repeat=False)
    if savedir is not None:
        ani.save(savedir / "trefoil_2d_animation.gif", fps=30, dpi=100)
    if show:
        plt.show()

def animate_3d_trefoil(savedir=None, show=True):
    n_frames = 180*3
    dazim = 360/n_frames

    x, y, z = trefoil(0, n_steps=N_STEPS)
    cmap = matplotlib.colormaps['viridis']
    axstereo = AxesStereo3D(ipd=IPD)
    axstereo.plot(x, y, z, c='k', alpha=0.2)
    scatter = axstereo.scatter(x, y, z, s=10)
    colors = cmap(_cmap_norm(z))
    for scat in scatter:
        scat.set_edgecolor(colors)
        scat.set_facecolor(colors)
    axstereo.fig.set_size_inches(6, 3)

    def animate(frame):
        x, y, z = trefoil(frame, n_steps=N_STEPS)
        colors = cmap(_cmap_norm(z))
        for scat in scatter:
            scat._offsets3d = (x, y, z)
            scat.set_edgecolor(colors)
            scat.set_facecolor(colors)
        axstereo.ax_left.view_init(azim=(axstereo.ax_left.azim + dazim), share=True)
        return scatter

    ani = animation.FuncAnimation(axstereo.fig, animate, frames=np.arange(n_frames),
                                  interval=20, repeat=False)
    if savedir is not None:
        ani.save(savedir / "trefoil_3d_animation.gif", fps=30, dpi=100)
    if show:
        plt.show()

def gen_logo(savedir=None, show=True):
    # Roll the data so the end of the line is in a straight portion
    x, y, z = trefoil()
    n_roll = 15
    x = np.roll(x, n_roll)
    y = np.roll(y, n_roll)
    z = np.roll(z, n_roll)

    axstereo = AxesAnaglyph(ipd=150, zscale=1, zzero=min(z))
    axstereo.plot(x, y, z, linewidth=6)
    axstereo.set_xlim(-4.5, 4.5)
    axstereo.set_ylim(-4.5, 4.5)
    axstereo.set_aspect('equal', adjustable='box')
    axstereo.set_axis_off()
    axstereo.fig.set_size_inches(1, 1)
    if savedir is not None:
        plt.savefig(savedir / 'mpl_stereo_logo.png', dpi=640)
    if show:
        plt.show()
    return axstereo

def gen_logo_with_text(savedir=None, show=True):
    axstereo = gen_logo(savedir=None, show=False)
    axstereo.ax.text(6, 0, 'mpl_stereo', fontsize=30, weight='bold', ha='left', va='center')
    axstereo.set_xlim(-4.5, 4.5+9*3.5)
    axstereo.fig.set_size_inches(1+3, 1)
    if savedir is not None:
        plt.savefig(savedir / 'mpl_stereo_logo_with_text.png', bbox_inches='tight', dpi=640)
    if show:
        plt.show()

def plot_2d_sun(savedir=None, show=True):
    sun_left_data, sun_right_data = sun_left_right()

    axstereo = AxesStereo2D()
    if IPD > 0:
        axstereo.ax_left.imshow(sun_left_data, cmap='gray')
        axstereo.ax_right.imshow(sun_right_data, cmap='gray')
    else:
        axstereo.ax_right.imshow(sun_left_data, cmap='gray')
        axstereo.ax_left.imshow(sun_right_data, cmap='gray')
    axstereo.fig.set_size_inches(6, 3)
    if savedir is not None:
        plt.savefig(savedir / 'sun_2d.png', bbox_inches='tight', dpi=640)
    if show:
        plt.show()

def wiggle_sun(savedir=None, show=True):
    sun_left_data, sun_right_data = sun_left_right()

    axstereo = AxesStereo2D()
    axstereo.ax_left.imshow(sun_left_data, cmap='gray')
    axstereo.ax_right.imshow(sun_right_data, cmap='gray')
    fig, ax = plt.subplots()
    fig.set_size_inches(3, 3)
    if savedir is not None:
        axstereo.wiggle(savedir / 'sun_wiggle.gif', ax=ax, dpi=640)

def plot_anaglyph_sun(savedir=None, show=True):
    sun_left_data, sun_right_data = sun_left_right()
    axstereo = AxesAnaglyph()
    axstereo.imshow_stereo(sun_left_data, sun_right_data, cmap='gray')
    axstereo.fig.set_size_inches(3, 3)
    if savedir is not None:
        plt.savefig(savedir / 'sun_anaglyph.png', bbox_inches='tight', dpi=640)
    if show:
        plt.show()

def plot_2d_church(savedir=None, show=True):
    church_left_data, church_right_data = church_left_right()

    axstereo = AxesStereo2D()
    if IPD > 0:
        axstereo.ax_left.imshow(church_left_data)
        axstereo.ax_right.imshow(church_right_data)
    else:
        axstereo.ax_right.imshow(church_left_data)
        axstereo.ax_left.imshow(church_right_data)
    axstereo.fig.set_size_inches(8, 3)
    if savedir is not None:
        plt.savefig(savedir / 'church_2d.png', bbox_inches='tight', dpi=640)
    if show:
        plt.show()

def wiggle_church(savedir=None, show=True):
    church_left_data, church_right_data = church_left_right()

    axstereo = AxesStereo2D()
    axstereo.ax_left.imshow(church_left_data)
    axstereo.ax_right.imshow(church_right_data)
    axstereo.fig.set_size_inches(4, 3)
    fig, ax = plt.subplots()
    if savedir is not None:
        axstereo.wiggle(savedir / 'church_wiggle.gif', ax=ax, dpi=640)

def plot_anaglyph_church(savedir=None, show=True):
    church_left_data, church_right_data = church_left_right()
    axstereo = AxesAnaglyph()
    axstereo.imshow_stereo(church_left_data, church_right_data)
    axstereo.fig.set_size_inches(4, 3)
    if savedir is not None:
        plt.savefig(savedir / 'church_anaglyph.png', bbox_inches='tight', dpi=640)
    if show:
        plt.show()

def stereo_square_2d_trefoil(savedir=None, show=True):
    x, y, z = trefoil()
    stereosquare = StereoSquare2D()
    stereosquare.plot(x, y, z, c='k', alpha=0.2)
    stereosquare.scatter(x, y, z, c=z, cmap='viridis', s=10)
    stereosquare.fig.set_size_inches(6, 6)
    if savedir is not None:
        stereosquare.wiggle(savedir / 'trefoil_2d_square.gif', dpi=100)
    if show:
        plt.show()

def stereo_square_3d_trefoil(savedir=None, show=True):
    x, y, z = trefoil()
    stereosquare = StereoSquare3D()
    stereosquare.plot(x, y, z, c='k', alpha=0.2)
    stereosquare.scatter(x, y, z, c=z, cmap='viridis', s=10)
    stereosquare.fig.set_size_inches(6, 6)
    if savedir is not None:
        stereosquare.wiggle(savedir / 'trefoil_3d_square.gif', dpi=100)
    if show:
        plt.show()

def stereo_square_church(savedir=None, show=True):
    church_left_data, church_right_data = church_left_right()
    stereosquare = StereoSquare2D()
    stereosquare.imshow_stereo(church_left_data, church_right_data, crop=True)
    stereosquare.fig.set_size_inches(8, 6)
    stereosquare.fig.set_dpi(320)
    if savedir is not None:
        stereosquare.wiggle(savedir / 'church_2d_square.gif')
    if show:
        plt.show()


def main():
    currdir = Path(__file__).parent.resolve()
    savedir = currdir
    show = False

    plot_2d_trefoil(savedir, show)
    plot_2d_sun(savedir, show)
    plot_2d_church(savedir, show)
    plot_3d_trefoil(savedir, show)
    plot_anaglyph_trefoil(savedir, show)
    plot_anaglyph_trefoil_zzero(savedir, show)
    plot_anaglyph_sun(savedir, show)
    plot_anaglyph_church(savedir, show)
    wiggle_2d_trefoil(savedir, show)
    wiggle_3d_trefoil(savedir, show)
    wiggle_sun(savedir, show)
    wiggle_church(savedir, show)
    animate_2d_trefoil(savedir, show)
    animate_3d_trefoil(savedir, show)
    stereo_square_2d_trefoil(savedir, show)
    stereo_square_3d_trefoil(savedir, show)
    stereo_square_church(savedir, show)
    gen_logo(savedir, show)
    gen_logo_with_text(savedir, show)


if __name__ == '__main__':
    main()
