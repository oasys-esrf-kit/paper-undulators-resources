
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

def q_a(x):
    # Tanaka & Kitamura 2009 equation (17), forzed to give 1 in the limit close to zero.
    if x > 1e-5:
        f_1 = -1 + numpy.exp(-2 * numpy.square(x)) + numpy.sqrt(2 * numpy.pi) * x * erf(numpy.sqrt(2) * x)
        value = numpy.sqrt(2 * numpy.square(x) / f_1)
    elif x < 1e-5 and x >= 0:
        value = 1.0
    else:
        raise RuntimeError('ERROR: Please provide a positive energy spread')

    return value


def q_s(x, factor=0.5):
    # Tanaka & Kitamura 2009 equation (24), please noticed the correction factor
    # which in our case of using Onuki&Elleaume should be factor = 0.5
    return 2 * numpy.power(q_a(x / 4), 2 / 3) * factor


def norm_energ_spr(energy_spread, harmonic_number=1, number_of_periods=111.11, verbose=1):
    out = 2 * numpy.pi * harmonic_number * number_of_periods * energy_spread
    if verbose:
        print("n=%d, N=%f, sigma_delta=%g; Normalized energy spread = %g" %
              (harmonic_number, number_of_periods, energy_spread, out))
    return out

def wofry_results(harmonic_number=1, shift=0.0, do_plot=0):

    #
    # Wofry results
    #
    Electron_energy = numpy.linspace(5.9, 6.1, 201)
    if shift == 0.0:
        filename = "enerSpreadScanWofry1D/wfr1D_elec_energy_scan_farfield_n%d.h5" % harmonic_number
    else:
        filename = "enerSpreadScanWofry1D/wfr1D_elec_energy_scan_farfield_n%d_shift%d.h5" % (harmonic_number, shift)

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

    for i in range(npoints): # scan in spread values
        Qa1 = numpy.zeros_like(spread)

        for i in range(spread.size):
            Qa1[i] = q_a(norm_energ_spr(spread[i], harmonic_number=harmonic_number, verbose=0))

    return spread.copy(), (FWHM / FWHM[1]).copy(), (SD / SD[1]).copy(), Qa1.copy()

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image


    if False:
        Spread1, Fwhm1, Sd1, Qa1 = wofry_results(harmonic_number=1, do_plot=0)
        Spread3, Fwhm3, Sd3, Qa3 = wofry_results(harmonic_number=3, do_plot=0)
        Spread5, Fwhm5, Sd5, Qa5 = wofry_results(harmonic_number=5, do_plot=0)
        Spread7, Fwhm7, Sd7, Qa7 = wofry_results(harmonic_number=7, do_plot=0)

        plot(Spread1, Fwhm1,
             Spread1, Sd1,
             Spread1, Qa1,
             Spread3, Fwhm3,
             Spread3, Sd3,
             Spread3, Qa3,
             Spread5, Fwhm5,
             Spread5, Sd5,
             Spread5, Qa5,
             Spread7, Fwhm7,
             Spread7, Sd7,
             Spread7, Qa7,
             xtitle="spread ", ytitle="correction for angular width", legend =[
                                                     "n=1 FWHM", "n=1 SD", "n=1 Qa",
                                                     "n=3 FWHM", "n=3 SD", "n=3 Qa",
                                                     "n=5 FWHM", "n=5 SD", "n=5 Qa",
                                                     "n=7 FWHM", "n=7 SD", "n=7 Qa",
                                                     ], title="", show=1)

    Spread1,   Fwhm1, Sd1, Qa1 = wofry_results(harmonic_number=1, do_plot=0)
    Spread1s,  Fwhm1s, Sd1s, Qa1s = wofry_results(harmonic_number=1, shift=-25, do_plot=0)
    Spread1s2, Fwhm1s2, Sd1s2, Qa1s2 = wofry_results(harmonic_number=1, shift=25, do_plot=0)
    Spread1s3, Fwhm1s3, Sd1s3, Qa1s3 = wofry_results(harmonic_number=1, shift=-50, do_plot=0)
    Spread1s4, Fwhm1s4, Sd1s4, Qa1s4 = wofry_results(harmonic_number=1, shift=+50, do_plot=0)

    plot(Spread1, Fwhm1,
         Spread1, Sd1,
         Spread1, Qa1,
         # Spread1s, Fwhm1s,
         # Spread1s, Sd1s,
         # Spread1s, Qa1s,
         # Spread1s2, Fwhm1s2,
         # Spread1s2, Sd1s2,
         # Spread1s2, Qa1s2,
         Spread1s3, Fwhm1s3,
         Spread1s3, Sd1s3,
         Spread1s3, Qa1s3,
         Spread1s4, Fwhm1s4,
         Spread1s4, Sd1s4,
         Spread1s4, Qa1s4,
         xtitle="spread ", ytitle="correction for angular width", legend=[
            "n=1 FWHM", "n=1 SD", "n=1 Qa",
            # "n=1 shifted -25 FWHM", "n=1 shifted -25 SD", "n=1 shifted -25 Qa",
            # "n=1 shifted +25 FWHM", "n=1 shifted +25 SD", "n=1 shifted +25 Qa",
            "n=1 shifted -50 FWHM", "n=1 shifted -50 SD", "n=1 shifted -50 Qa",
            "n=1 shifted +50 FWHM", "n=1 shifted +50 SD", "n=1 shifted +50 Qa",
        ], title="", show=1)
    #
    # spectra results
    #

    # if False:
    #
    #     datadir = "/scisoft/data/srio/paper-undulator/spectra_results/"
    #     if case == 1:
    #         Shift = numpy.arange(-400, 100, 1) # case 1
    #         filename = datadir + "Spectra_Divergence_Off_Res_ESRF_ID06_EBS_CPMU18_1_0_emitt_0_spread_case1.hdf5"
    #     elif case == 2:
    #         Shift = numpy.arange(-1000, 1000, 2) # case 2
    #         filename = datadir + "Spectra_Divergence_Off_Res_ESRF_ID06_EBS_CPMU18_1_0_emitt_0_spread_case2.hdf5"
    #
    #     import h5py
    #
    #     file = h5py.File(filename, 'r')
    #     shift =     file['/Scan/harmonic 01/shift'][()]
    #     flux_dens = file["/Scan/harmonic 01/flux_dens"][()]
    #     x_div =         file["/Scan/harmonic 01/x_div"][()]
    #     y_div =         file["/Scan/harmonic 01/y_div"][()]
    #
    #     print(shift.shape, flux_dens.shape, x_div.shape, y_div.shape)
    #
    #     nx = flux_dens.shape[1]
    #
    #     # plot(y_div, flux_dens[0, nx // 2, :], xtitle="Angle in mrad")
    #
    #     YS = numpy.zeros_like(Shift, dtype=float)
    #     SD = numpy.zeros_like(Shift, dtype=float)
    #     for i, shift in enumerate(Shift):
    #         ys = flux_dens[i, nx // 2, :]
    #         fwhm, _, _ = get_fwhm(ys, y_div)
    #         # plot(x_div, ys, title="shift = %d eV FWHM = %f" % (shift, fwhm), xtitle="Angle in mrad")
    #         YS[i] = fwhm * 1e-3 * 100 # in m @ 100m
    #         SD[i] = get_stdev(y_div * 1e-3 * 100, ys)
    #
    #     plot(Shift, YS,
    #          Shift, SD * 2.355,
    #          xtitle="Shift [eV]", ytitle='SPECTRA FWHM % 100m [m]', legend=['FWHM','SD*2.355'], show=0)
    #     Delta = 1 * 111.111 * (Shift) / 10000.0
    #     plot(Delta, YS,
    #          Delta, SD * 2.355,
    #          xtitle="Delta", ytitle='SPECTRA FWHM % 100m [m]', legend=['FWHM','SD*2.355'], show=0)
    #
    #     plot(Delta, YS,
    #          DeltaW, YW,
    #          Delta,  SD,
    #          DeltaW, SDW,
    #          xtitle="Delta", legend=['SPECTRA FWHM','WOFRY FWHM','SPECTRA SD*2.355','WOFRY SD*2.355'])