"""
Typical usage
=================

To demonstrate the use of the Hankel Transform class, we will give an example
of propagating a radially-symmetric beam using the beam propagation method.

In this case, it will be a simple Gaussian beam propagating way from focus and
diverging.

First we will use a loop over :math:`z` position, and then we will demonstrate
the vectorisation of the :func:`.HankelTransforms.iqdht` (and
:func:`~.HankelTransforms.qdht`) functions.

"""

from srxraylib.plot.gol import set_qt, plot, plot_image
from pyhank import HankelTransform
import numpy as np
import matplotlib.pyplot as plt

# 1D Gaussian function
def gauss1d(x, x0, fwhm):
    return np.exp(-2 * np.log(2) * ((x - x0) / fwhm) ** 2)

def slit1d(x, x0, fwhm):
    y = np.ones_like(x)
    y[ (np.abs(x - x0) > fwhm)] = 0
    return y

def get_initial_efield(
                        Dr=100e-6 / 2,  # Beam radius (100um)
                        lambda_ = 1.5e-10,  # 488e-9  # wavelength 488nm
                        nr = 1024 * 5,  # Number of sample points
                        r_max = 5e-3 * 5, # Maximum radius (5mm)
                        ):

    # Initialise radius grid

    r = np.linspace(0, r_max, nr)

    # Set up beam parameters


    # Set up the electric field profile at :math:`z = 0`, and resample onto the correct radial grid
    # (``transformer.r``) as required for the QDHT.
    # Er = gauss1d(r, 0, Dr)  + 0j  # Initial field
    Er = slit1d(r, 0, Dr)  + 0j  # Initial field

    return Er, r, lambda_

def hankel_propagate(Er, r, lambda_, z_max=100):
    # Initialise :math:`z` grid

    Nz = 5  # Number of z positions
    z = np.linspace(0, z_max, Nz)  # Propagation axis


    H = HankelTransform(order=0, radial_grid=r)  # Set up an object, telling it the order (``0``) and the radial grid.
    ErH = H.to_transform_r(Er)  # Resampled field


    # Perform Hankel Transform
    # Convert from physical field to physical wavevector
    EkrH = H.qdht(ErH)

    # Propagate the beam - loop
    # Do the propagation in a loop over :math:`z`

    # Pre-allocate an array for field as a function of r and z
    # Erz = np.zeros((r.size, Nz), dtype=complex)
    k0 = 2 * np.pi / lambda_  # Vacuum k vector
    kz = np.sqrt(k0 ** 2 - H.kr ** 2)

    # for i, z_loop in enumerate(z):
    phi_z = kz * z_max  # Propagation phase
    EkrHz = EkrH * np.exp(1j * phi_z)  # Apply propagation
    ErHz = H.iqdht(EkrHz)  # iQDHT
    Erz = H.to_original_r(ErHz)  # Interpolate output

    return Erz


if __name__ == "__main__":
    set_qt()

    Dr = 100e-6 / 2  # Beam radius (100um)
    lambda_ = 1.5e-10  # 488e-9  # wavelength 488nm
    nr = 1024 * 5  # Number of sample points
    r_max = 5e-3 * 5  # Maximum radius (5mm)
    z_max = 100.0
    Er, r, lambda1 = get_initial_efield(Dr=Dr, lambda_=lambda_, nr=nr, r_max=r_max)

    # H = HankelTransform(order=0, radial_grid=r)  # Set up an object, telling it the order (``0``) and the radial grid.
    Erz = hankel_propagate(Er, r, lambda1, z_max=z_max)

    Irz = np.abs(Erz) ** 2

    # theory
    from scipy.special import jv

    sin_theta_array = np.sin(r / z_max)
    aperture_diameter = 2 * Dr
    x = (2 * np.pi / lambda_) * (aperture_diameter / 2) * sin_theta_array
    x_over_pi = x / np.pi
    electric_field = 2 * jv(1, x) / x
    intensity = electric_field ** 2

    plot(np.degrees(r / z_max), np.abs(Er) ** 2,
         np.degrees(r / z_max), Irz / Irz.max(),
         np.degrees(r / z_max), intensity,
         xrange=[0., 0.0003], legend=['data', 'Hankel normalized', 'theory'])



