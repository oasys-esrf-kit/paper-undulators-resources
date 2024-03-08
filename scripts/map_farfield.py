
#
# Import section
#
import numpy
from scipy.special import erf


from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from srxraylib.plot.gol import plot


def get_fwhm(histogram, bins, ret0=None):
    fwhm = ret0
    quote = ret0
    coordinates = None

    if histogram.size > 1:
        quote = numpy.max(histogram)*0.5
        cursor = numpy.where(histogram >= quote)

        if histogram[cursor].size > 1:
            bin_size    = bins[1]-bins[0]
            fwhm        = bin_size*(cursor[0][-1]-cursor[0][0])
            coordinates = (bins[cursor[0][0]], bins[cursor[0][-1]])

    return fwhm, quote, coordinates

def get_stdev(x, y):
    delta = (x[1] - x[0])
    Y = y.copy()
    Y /= y.sum() * delta
    m1 = (x ** 1 * Y).sum() * delta
    m2 = (x ** 2 * Y).sum() * delta
    return numpy.sqrt(m2 - m1**2)


def wofry_results(harmonic_number=1, shift=0.0, do_plot=0):

    #
    # Wofry results
    #
    Electron_energy = numpy.linspace(5.9, 6.1, 201)
    dir = "/scisoft/data/srio/paper-undulator/enerSpreadScanWofry1D/"
    if shift == 0.0:
        filename = dir + "wfr1D_elec_energy_scan_farfield_n%d.h5" % harmonic_number
    else:
        filename = dir + "wfr1D_elec_energy_scan_farfield_n%d_shift%d.h5" % (harmonic_number, shift)

    Y = numpy.zeros_like(Electron_energy, dtype=float)
    SD = numpy.zeros_like(Electron_energy, dtype=float)

    for i, electron_energy in enumerate(Electron_energy):
        subgroupname = "wfr_%5.3f" % electron_energy

        output_wavefront = GenericWavefront1D.load_h5_file(filename, subgroupname)
        x = output_wavefront.get_abscissas()
        y = output_wavefront.get_intensity()
        fwhm, _, _ = get_fwhm(y, x)

        # plot(x, y, title="electron energy = %f GeV FWHM = %f" % (electron_energy, fwhm))

        Y[i] = fwhm
        SD[i] = get_stdev(x, y)

    if do_plot: plot(Electron_energy, Y,
         Electron_energy, SD*2.355,
         xtitle="electron energy [GeV]", ytitle='WOFRY FWHM % 100m [m]', legend=['FWHM','SD*2.355'], show=1)


    #
    # average
    #
    npoints = Electron_energy.size
    e0 = 6.0

    e_energies = Electron_energy
    spread = numpy.linspace(0, 0.005, npoints)  # (e_energies - e0) / e0  #

    Weights = numpy.zeros((npoints, npoints)) # spread, e energy

    print(">>>>", e_energies.shape, spread.shape)
    for i in range(npoints):
        Weights[i, :] = numpy.exp(-(e_energies - e0)**2 / 2 / (e0 * spread[i])**2 )

    if do_plot: plot_image(Weights, spread, e_energies, xtitle="sigma spread", ytitle="e energy [GeV]", title="",
                           aspect='auto', show=1)

    # store all results in an image
    for j, electron_energy in enumerate(Electron_energy):
        subgroupname = "wfr_%5.3f" % electron_energy
        output_wavefront = GenericWavefront1D.load_h5_file(filename, subgroupname)
        x = output_wavefront.get_abscissas()
        y = output_wavefront.get_intensity()  #* Weights[i, j]
        if j == 0:
            Y_IMG = numpy.zeros((Electron_energy.size, y.size))
        Y_IMG[j, :] = y

    if do_plot: plot_image(Y_IMG, Electron_energy, x * 1e3, ytitle="x [mm] @ far field (100m)", xtitle="e energy [GeV]", title="",
                           yrange=[-2,2], aspect='auto', show=1)

    FWHM = numpy.zeros(npoints)
    SD = numpy.zeros(npoints)
    for i in range(npoints):  # scan in spread values
        ys = Weights[i, :] # spread, e energy
        ys = ys / ys.sum()
        Y_IMG_weighted = Y_IMG * numpy.outer(ys, numpy.ones_like(x))

        yy = Y_IMG_weighted.sum(axis=0)
        fwhm, _, _ = get_fwhm(yy, x)
        SD[i] = get_stdev(x, yy)
        FWHM[i] = fwhm
        if i == 0:
            Y_IMG_WEIGHTED = numpy.zeros((spread.size, yy.size))
        Y_IMG_WEIGHTED[i, :] = yy

    if do_plot: plot_image(Y_IMG_WEIGHTED, spread, x * 1e3,
                           ytitle="x [mm] @ far field (100m) averaged", xtitle="e energy spread",
                           title="WEIGHTED shift=%d" % shift,
                           yrange=[-2,2], aspect='auto', show=1)

    return spread.copy(), (FWHM / FWHM[1]).copy(), (SD / SD[1]).copy()


if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image, plot_show

    # WOFRY on resonance results
    if False:
        Spread1, Fwhm1, Sd1 = wofry_results(harmonic_number=1, shift=0.0, do_plot=0)


        plot(Spread1, Fwhm1,
             Spread1, Sd1,
             xtitle="spread ", ytitle="correction for angular width", legend =[
                                                     "n=1 FWHM", "n=1 SD",
                                                     ], title="pySRU", yrange=[0,5], show=0)




    Shift = numpy.arange(-300, 501, 1)

    # Shift = [-300, -200, -100, 0, 100, 200, 300, 400]
    # Shift = [-50, -25, 0, 25, 50]

    Mfwhm_a = numpy.zeros( (201, Shift.size) )
    Mstdev_a = numpy.zeros_like(Mfwhm_a)

    for i, shift in enumerate(Shift):
        Spread1, Fwhm1, Sd1 = wofry_results(harmonic_number=1, shift=shift, do_plot=0)


        # plot(Spread1, Fwhm1,
        #      Spread1, Sd1,
        #      xtitle="spread ", ytitle="correction for angular width", legend =[
        #                                              "n=1 FWHM", "n=1 SD",
        #                                              ], title="pySRU shift=%d" % shift, show=1)

        Mfwhm_a[:, i] = Fwhm1
        Mstdev_a[:, i] = Sd1

    plot_image(Mfwhm_a, show=0)
    plot_image(Mstdev_a, show=0)
    plot_show()

    plot_image(Mfwhm_a,  Shift, Spread1, show=0)
    plot_image(Mstdev_a, Shift, Spread1,  show=0)

    plot_show()