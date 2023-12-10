import matplotlib.pyplot as plt
import numpy as np

from mpl_stereo import AxesStereo

def test_AxesStereo():
    # Smoke test
    axstereo = AxesStereo(ipd=0.2, focal_plane=0)
    assert True


def plotting_tests():
    # Plot up a sphere
    n = 30  # Number of points
    radius = 2  # Radius of the sphere
    offset = np.radians(6)

    # Create theta and phi angles for the sphere
    theta = np.linspace(0, np.pi, n)
    phi = np.linspace(0, 2*np.pi, n) + offset
    theta, phi = np.meshgrid(theta, phi)

    # Parametric equations for the sphere
    x = radius * np.sin(theta) * np.sin(phi)
    y = radius * np.cos(theta)
    z = radius * np.sin(theta) * np.cos(phi)
    # Make the stereoscopic plots
    axstereo = AxesStereo(ipd=0.2, focal_plane=0)
    axstereo.scatter(x, y, z, c=z, cmap='viridis', s=10)


if __name__ == '__main__':
    plotting_tests()