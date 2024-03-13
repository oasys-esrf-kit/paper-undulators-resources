import numpy
from scipy.optimize import curve_fit



def gauss(x, height, center, sigma):
    return height * numpy.exp(-(x - center)**2 / (2 * sigma**2))

def fitgauss(x, y, p0 = [1.0, 0, 0.5]):
    popt, pcov = curve_fit(gauss, x, y, p0=p0, method='trf')  # ‘lm’, ‘trf’, ‘dogbox’
    y1 = gauss(x, *popt)
    print("FWHM = ", 2.355 * popt[2])
    return popt[2], popt, pcov


#
# far field
#
if False:
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
          xtitle="Size @ 100 m [m]", ytitle="Normalized intensity", show=0)

#
# backpropagatopm
#
a1d = numpy.loadtxt("fit1Dwofry1D/id06_backpropagated_y.dat")
a2d = numpy.loadtxt("fit1Dwofry2D/id06_backpropagated_y.dat")
aSrw = numpy.loadtxt("fit1Dsrw/id06_backpropagated_y.dat")
aSpectra = numpy.loadtxt("spectra_windows/spectra_results/id06_backpropagated_spectra_1D_hor_zero_emitt_zero_spread.txt")

aTheory = numpy.loadtxt("fit1D/id06_backpropagated.dat")
aTheory2 = numpy.loadtxt("fit1D/undulator_size.dat")
aTheory2 = numpy.loadtxt("fit1D/undulator_size_twice_interval.dat")

import scipy.constants as codata
u_nperiods = 111.111
u_period = 0.018
wavelength = codata.h * codata.c / codata.e / 10000.
sigma = 2.704 / 4 / numpy.pi * numpy.sqrt(u_period * u_nperiods * wavelength)


from srxraylib.plot.gol import plot
cte = numpy.sqrt(wavelength * u_nperiods * u_period)

if True:
     plot(1e6 * a1d[:, 0], a1d[:, 1] / a1d[:, 1].max(),
          1e6 * a2d[:, 0], a2d[:, 1] / a2d[:, 1].max(),
          1e6 * aSrw[:, 0], aSrw[:, 1] / aSrw[:, 1].max(),
          1e3 * aSpectra[:, 0], aSpectra[:, 1] / aSpectra[:, 1].max(),
          1e6 * aTheory[:, 0], aTheory[:, 1] / aTheory[:, 1].max(),
          1e6 * aTheory[:, 0], numpy.exp(- aTheory[:, 0]**2 / (2 * sigma**2)),
          1e6 * cte * aTheory2[:, 0], aTheory2[:, 1] / aTheory2[:, 1].max(),
          legend=['1d', '2d-v', 'srw-v','spectra v', 'theory (4 pi)', 'Gaussian', 'Theory raw nb'], title='backpropagation',
          xtitle="Size @ 0 m [um]", ytitle="Normalized intensity", xrange=[-20,20], show=1)

     sigma2 = sigma # 0.40994 * cte
     # plot()
     plot(
          1e6 * aTheory[:, 0], aTheory[:, 1] / aTheory[:, 1].max(),
          1e6 * aTheory[:, 0], numpy.exp(- aTheory[:, 0]**2 / (2 * sigma2**2)),
          1e6 * cte * aTheory2[:, 0], aTheory2[:, 1] / aTheory2[:, 1].max(),
          legend=['theory (4 pi)', 'Gaussian', 'Theory raw nb'], title='backpropagation',
          marker=[None, None, '+'], linestyle=[None,None,''],
          xtitle="Size @ 0 m [um]", ytitle="Normalized intensity", xrange=[-20,20], show=1)

if False: # for Chubar
     plot(
          1e6 * aSrw[:, 0], aSrw[:, 1] / aSrw[:, 1].max(),
          1e3 * aSpectra[:, 0], aSpectra[:, 1] / aSpectra[:, 1].max(),
          1e6 * aTheory[:, 0], numpy.exp(- aTheory[:, 0]**2 / (2 * sigma**2)),
          legend=['srw-v','spectra v', 'Gaussian'], title='backpropagation',
          xtitle="Size @ 0 m [um]", ytitle="Normalized intensity", xrange=[-20,20], show=1)

#
# plot(
#      aTheory2[:, 0], aTheory2[:, 1] / aTheory2[:, 1].max(),
#      aTheory2[:, 0], numpy.exp(- aTheory2[:, 0]**2 / (2 * sigma2**2)),
#      legend=['Theory raw nb', 'fit'], title='backpropagation',
#      # marker=[None, None, '+'], linestyle=[None,None,''],xrange=[-20,20],
#      xtitle="Size @ 0 m [um]", ytitle="Normalized intensity", show=1)

if False:
     x = aTheory2[:, 0]
     y = aTheory2[:, 1]
     y = y / numpy.trapz(y, x)

     plot(x, y)

     # ss = numpy.sqrt( (y * x**2).sum() * (x[1] - x[0]))
     # ss = numpy.sqrt( numpy.trapz(x**2 * y, x)) # RMS numeric 0.262849885337894 BAD
     ss = 0.409924 / 2  # Mathematica divided by 2.. Why?? NOT BAD
     ss = 0.218042 #Elleume # NOT BAD
     ss, _, _ = fitgauss(x, y)  # fit 0.21034769697537756 NOT BAD
     print("stdev", ss)

     plot(x, y / y.max(),
          x, numpy.exp(- x**2 / (2 * ss**2)))


     # cte = numpy.sqrt(wavelength * u_nperiods * u_period)
     # plot(1e6 * cte * aTheory2[:, 0], aTheory[:, 1] / aTheory[:, 1].max(),)