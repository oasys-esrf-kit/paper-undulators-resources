import numpy

#
# far field
#
a1d = numpy.loadtxt("fit1Dwofry1D/id06_farfield_y.dat")
a2d = numpy.loadtxt("fit1Dwofry2D/id06_farfield_y.dat")
aSrw =         numpy.loadtxt(    "fit1Dsrw/id06_farfield_y.dat")
aTheory =          numpy.loadtxt("fit1D/id06_farfield.dat")

import scipy.constants as codata
u_nperiods = 111.111
u_period = 0.018
wavelength = codata.h * codata.c / codata.e / 10000.
sigmap = 0.690 * numpy.sqrt(wavelength / (u_period * u_nperiods))
x = numpy.linspace(-0.003, 0.003, 100) / 100.0
gaussian_div = numpy.exp(-x**2 / (2 * sigmap**2))


from srxraylib.plot.gol import plot

plot(a1d[:, 0], a1d[:, 1] / a1d[:, 1].max(),
     a2d[:, 0], a2d[:, 1] / a2d[:, 1].max(),
     aSrw[:, 0], aSrw[:, 1] / aSrw[:, 1].max(),
     aTheory[:, 0], aTheory[:, 1] / aTheory[:, 1].max(),
     x * 100, gaussian_div,
     legend=['1d', '2d-v', 'srw-v', 'theory', 'Gaussian'], title="farfield",
     xtitle="Size @ 100 m [m]", ytitle="Normalized intensity")

#
# backpropagatopm
#
a1d = numpy.loadtxt("fit1Dwofry1D/id06_backpropagated_y.dat")
a2d = numpy.loadtxt("fit1Dwofry2D/id06_backpropagated_y.dat")
aSrw = numpy.loadtxt("fit1Dsrw/id06_backpropagated_y.dat")
aSpectra = numpy.loadtxt("spectra_windows/spectra_results/id06_backpropagated_spectra_1D_hor_zero_emitt_zero_spread.txt")

aTheory = numpy.loadtxt("fit1D/id06_backpropagated.dat")

import scipy.constants as codata
u_nperiods = 111.111
u_period = 0.018
wavelength = codata.h * codata.c / codata.e / 10000.
sigma = 2.704 / 4 / numpy.pi * numpy.sqrt(u_period * u_nperiods * wavelength)


from srxraylib.plot.gol import plot

plot(1e6 * a1d[:, 0], a1d[:, 1] / a1d[:, 1].max(),
     1e6 * a2d[:, 0], a2d[:, 1] / a2d[:, 1].max(),
     1e6 * aSrw[:, 0], aSrw[:, 1] / aSrw[:, 1].max(),
     1e3 * aSpectra[:, 0], aSpectra[:, 1] / aSpectra[:, 1].max(),
     1e6 * aTheory[:, 0], aTheory[:, 1] / aTheory[:, 1].max(),
     1e6 * aTheory[:, 0], numpy.exp(- aTheory[:, 0]**2 / (2 * sigma**2)),
     legend=['1d', '2d-v', 'srw-v','spectra v', 'theory (4 pi)', 'Gaussian'], title='backpropagation',
     xtitle="Size @ 0 m [um]", ytitle="Normalized intensity")