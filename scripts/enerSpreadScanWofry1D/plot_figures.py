
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

def wofry_results(harmonic_number=1, shift=0.0, do_plot=0, figfile=""):

    #
    # Wofry results
    #
    Electron_energy = numpy.linspace(5.9, 6.1, 201)
    if shift == 0.0:
        filename = "wfr1D_elec_energy_scan_farfield_n%d.h5" % harmonic_number
    else:
        filename = "wfr1D_elec_energy_scan_farfield_n%d_shift%d.h5" % (harmonic_number, shift)

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

    if do_plot >= 3: plot(Electron_energy, Y,
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

    if do_plot >= 2: plot_image(Weights, spread, e_energies,
                           xtitle="sigma spread", ytitle=r"Electron energy $E_e$ [GeV]", title="",
                           figsize=(9,7), aspect='auto', show=1)

    # store all results in an image
    for j, electron_energy in enumerate(Electron_energy):
        subgroupname = "wfr_%5.3f" % electron_energy
        output_wavefront = GenericWavefront1D.load_h5_file(filename, subgroupname)
        x = output_wavefront.get_abscissas()
        y = output_wavefront.get_intensity()  #* Weights[i, j]
        if j == 0:
            Y_IMG = numpy.zeros((Electron_energy.size, y.size))
        Y_IMG[j, :] = y

    if do_plot >= 2:
        plot_image(Y_IMG, Electron_energy, x * 1e3,
                           ytitle="y [mm] @ far field (100m)", xtitle=r"Electron energy $E_e$ [GeV]", title="",
                           yrange=[-2,2], figsize=(9,7), aspect='auto', show=0)
        # plt.savefig("detuning-map-electron.png")
        # print("File written to disk: detuning-map-electron.png")
        plot_show()

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

    if do_plot >=1:
        plot_image(Y_IMG_WEIGHTED, spread, x * 1e3,
                           ytitle="y [mm] @ far field (100m) averaged",
                           xtitle=r"Electron energy spread $\delta_\mathcal{E}$",
                           title="", #"WEIGHTED shift=%d" % shift,
                           yrange=[-2,2], figsize=(9,7), aspect='auto', show=0)

        if figfile != "":
            plt.savefig(figfile)
            print("File written to disk: %s" % figfile)
            plot_show()

    if do_plot >= 3:
        plot(x * 1e3, Y_IMG_WEIGHTED[0, :]          / Y_IMG_WEIGHTED[0, :]         .max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*1, :] / Y_IMG_WEIGHTED[(200//5)*1, :].max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*2, :] / Y_IMG_WEIGHTED[(200//5)*2, :].max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*3, :] / Y_IMG_WEIGHTED[(200//5)*3, :].max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*4, :] / Y_IMG_WEIGHTED[(200//5)*4, :].max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*5, :] / Y_IMG_WEIGHTED[(200//5)*5, :].max(),
             legend=[r'$\delta_\mathcal{E}$=0.000',r'$\delta_\mathcal{E}$=0.001',r'$\delta_\mathcal{E}$=0.002',r'$\delta_\mathcal{E}$=0.003',r'$\delta_\mathcal{E}$=0.004',r'$\delta_\mathcal{E}$=0.005'])

    for i in range(npoints): # scan in spread values
        Qa1 = numpy.zeros_like(spread)

        for i in range(spread.size):
            Qa1[i] = q_a(norm_energ_spr(spread[i], harmonic_number=harmonic_number, verbose=0))

    return spread.copy(), (FWHM / FWHM[1]).copy(), (SD / SD[1]).copy(), Qa1.copy()



def wofry_results_backpropagated(harmonic_number=1, shift=0.0, do_plot=0, figfile=""):

    #
    # Wofry results
    #
    Electron_energy = numpy.linspace(5.9, 6.1, 201)
    if shift == 0.0:
        filename = "wfr1D_elec_energy_scan_backpropagated_n%d.h5" % harmonic_number
    else:
        filename = "wfr1D_elec_energy_scan_backpropagated_n%d_shift%d.h5" % (harmonic_number, shift)

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

    if do_plot >= 3: plot(Electron_energy, Y,
         Electron_energy, SD*2.355,
         xtitle="electron energy [GeV]", ytitle='BACKPROPAGATED WOFRY FWHM % 100m [m]', legend=['FWHM','SD*2.355'], show=1)


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

    if do_plot >= 2: plot_image(Weights, spread, e_energies,
                           xtitle="sigma spread", ytitle=r"Electron energy $E_e$ [GeV]", title="",
                           figsize=(9,7), aspect='auto', show=1)

    # store all results in an image
    for j, electron_energy in enumerate(Electron_energy):
        subgroupname = "wfr_%5.3f" % electron_energy
        output_wavefront = GenericWavefront1D.load_h5_file(filename, subgroupname)
        x = output_wavefront.get_abscissas()
        y = output_wavefront.get_intensity()  #* Weights[i, j]
        if j == 0:
            Y_IMG = numpy.zeros((Electron_energy.size, y.size))
        Y_IMG[j, :] = y

    if do_plot >= 2:
        plot_image(Y_IMG, Electron_energy, x * 1e3,
                           ytitle="y [mm] @ backpropagated", xtitle=r"Electron energy $E_e$ [GeV]", title="",
                           yrange=[-0.05,0.05], figsize=(9,7), aspect='auto', show=0)
        plot_show()

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

    if do_plot >= 1:
        plot_image(Y_IMG_WEIGHTED, spread, x * 1e3,
                           ytitle="y [mm] @ backpropagated averaged",
                           xtitle=r"Electron energy spread $\delta_\mathcal{E}$",
                           title="", #"WEIGHTED shift=%d" % shift,
                           yrange=[-0.05,0.05], figsize=(9,7), aspect='auto', show=0)
        if figfile != "":
            plt.savefig(figfile)
            print("File written to disk: %s" % figfile)
        plot_show()

    if do_plot >= 3:

        for i in range(30):
            print(i, Y_IMG_WEIGHTED[i, :].shape, Y_IMG_WEIGHTED[i, :])
        # plot(x, Y_IMG_WEIGHTED[0, :], marker='+', show=0)

        plot(x * 1e3, Y_IMG_WEIGHTED[1, :]          / Y_IMG_WEIGHTED[1, :].max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*1, :] / Y_IMG_WEIGHTED[(200//5)*1, :].max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*2, :] / Y_IMG_WEIGHTED[(200//5)*2, :].max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*3, :] / Y_IMG_WEIGHTED[(200//5)*3, :].max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*4, :] / Y_IMG_WEIGHTED[(200//5)*4, :].max(),
             x * 1e3, Y_IMG_WEIGHTED[(200//5)*5, :] / Y_IMG_WEIGHTED[(200//5)*5, :].max(),
             legend=[r'$\delta_\mathcal{E}$=0.000',r'$\delta_\mathcal{E}$=0.001',r'$\delta_\mathcal{E}$=0.002',r'$\delta_\mathcal{E}$=0.003',r'$\delta_\mathcal{E}$=0.004',r'$\delta_\mathcal{E}$=0.005'],
             xrange=[-0.025, 0.025], marker=['+',None,None,None,None,None])

    for i in range(npoints): # scan in spread values
        Qs1 = numpy.zeros_like(spread)

        for i in range(spread.size):
            Qs1[i] = q_s(norm_energ_spr(spread[i], harmonic_number=harmonic_number, verbose=0))

    return spread.copy(), (FWHM / FWHM[1]).copy(), (SD / SD[1]).copy(), Qs1.copy()

def spectra_results(harmonic_number=1, shift=0.0, do_plot=0):
    datadir = "/scisoft/data/srio/paper-undulator/spectra_results/"
    filename = datadir + "Spectra_div_ESRF_ID06_EBS_CPMU18_1_0_emitt_0_spread_diff_elec_energy_harm%s.hdf5" % harmonic_number

    import h5py

    print(filename)
    file = h5py.File(filename, 'r')
    # //scisoft/data/srio/paper-undulator/spectra_results/Spectra_div_ESRF_ID06_EBS_CPMU18_1_0_emitt_0_spread_diff_elec_energy_harm1.hdf5::/Scan/harmonic 01
    Electron_energy =     file['/Scan/harmonic 0%d/electron_energy'% harmonic_number][()]
    flux_dens =           file["/Scan/harmonic 0%d/flux_dens"      % harmonic_number][()]
    x_div =               file["/Scan/harmonic 0%d/x_div"          % harmonic_number][()]
    y_div =               file["/Scan/harmonic 0%d/y_div"          % harmonic_number][()]

    print(Electron_energy.shape, flux_dens.shape, x_div.shape, y_div.shape)

    nx = flux_dens.shape[1]

    # plot(y_div, flux_dens[0, nx // 2, :], xtitle="Angle in mrad")

    # Electron_energy = numpy.linspace(5.9, 6.1, 201)

    YS = numpy.zeros_like(Electron_energy, dtype=float)
    SD = numpy.zeros_like(Electron_energy, dtype=float)
    for i, electron_energy in enumerate(Electron_energy):
        ys = flux_dens[i, nx // 2, :]
        fwhm, _, _ = get_fwhm(ys, y_div)
        # plot(x_div, ys, title="shift = %d eV FWHM = %f" % (shift, fwhm), xtitle="Angle in mrad")
        if fwhm is not None: YS[i] = fwhm * 1e-3 * 100 # in m @ 100m
        SD[i] = get_stdev(y_div * 1e-3 * 100, ys)

    if do_plot: plot(Electron_energy, YS,
         Electron_energy, SD * 2.355,
         xtitle="Electron_energy [GeV]", ytitle='SPECTRA FWHM % 100m [m]', legend=['FWHM','SD*2.355'], show=1)


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
        # subgroupname = "wfr_%5.3f" % electron_energy
        # output_wavefront = GenericWavefront1D.load_h5_file(filename, subgroupname)
        # x = output_wavefront.get_abscissas()
        # y = output_wavefront.get_intensity()  #* Weights[i, j]

        y = flux_dens[j, nx // 2, :]
        x = y_div * 1e-3 * 100

        if j == 0:
            Y_IMG = numpy.zeros((Electron_energy.size, y.size))
        Y_IMG[j, :] = y

    if do_plot: plot_image(Y_IMG, Electron_energy, x * 1e3, ytitle="y [mm] @ far field (100m)", xtitle="e energy [GeV]", title="",
                           yrange=[-2,2],
                           aspect='auto', show=1)

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
                           ytitle="y [mm] @ far field (100m) averaged", xtitle="e energy spread",
                           title="WEIGHTED shift=%d" % shift,
                           yrange=[-2,2], aspect='auto', show=1)

    for i in range(npoints): # scan in spread values
        Qa1 = numpy.zeros_like(spread)

        for i in range(spread.size):
            Qa1[i] = q_a(norm_energ_spr(spread[i], harmonic_number=harmonic_number, verbose=0))

    return spread.copy(), (FWHM / FWHM[1]).copy(), (SD / SD[1]).copy(), Qa1.copy()


def plot_farfield(figfile=""):
    # FIG 4 : WOFRY on resonance results
    Spread1, Fwhm1, Sd1, Qa1 = wofry_results(harmonic_number=1, do_plot=0)
    Spread3, Fwhm3, Sd3, Qa3 = wofry_results(harmonic_number=3, do_plot=0)
    Spread5, Fwhm5, Sd5, Qa5 = wofry_results(harmonic_number=5, do_plot=0)

    plot(Spread1, Fwhm1,
         Spread1, Sd1,
         Spread1, Qa1,
         Spread3, Fwhm3,
         Spread3, Sd3,
         Spread3, Qa3,
         Spread5, Fwhm5,
         Spread5, Sd5,
         Spread5, Qa5,
         xtitle=r"$\delta_\mathcal{E}$", ytitle="correction factor for angular width", legend =[
                                                 "n=1 FWHM", "n=1 SD", "n=1 Qa",
                                                 "n=3 FWHM", "n=3 SD", "n=3 Qa",
                                                 "n=5 FWHM", "n=5 SD", "n=5 Qa",
                                                 ],
         linestyle=['dashed', 'dotted', None,
                    'dashed', 'dotted', None,
                    'dashed', 'dotted', None,
                    ],
         color=['black', 'black', 'black',
                'blue', 'blue', 'blue',
                'red', 'red', 'red',
                ],
         title="", yrange=[0,5], figsize=(12,8), show=0)

    if figfile != "":
        plt.savefig(figfile)
        print("File writte to disk: %s" % figfile)

    plot_show()


def plot_backpropagated(figfile=""):
    # FIG 4 : WOFRY on resonance results
    Spread1, Fwhm1, Sd1, Qs1 = wofry_results_backpropagated(harmonic_number=1, do_plot=0)
    Spread3, Fwhm3, Sd3, Qs3 = wofry_results_backpropagated(harmonic_number=3, do_plot=0)
    Spread5, Fwhm5, Sd5, Qs5 = wofry_results_backpropagated(harmonic_number=5, do_plot=0)

    plot(
         Spread1, Fwhm1,
         Spread1, Sd1,
         Spread1, Qs1,
         Spread3, Fwhm3,
         Spread3, Sd3,
         Spread3, Qs3,
         Spread5, Fwhm5,
         Spread5, Sd5,
         Spread5, Qs5,
         xtitle=r"$\delta_\mathcal{E}$", ytitle="correction factor for width", legend =[
                                                 "n=1 FWHM", "n=1 SD", "n=1 Qs",
                                                 "n=3 FWHM", "n=3 SD", "n=3 Qs",
                                                 "n=5 FWHM", "n=5 SD", "n=5 Qs",
                                                 ],
        linestyle=['dashed', 'dotted', None,
                   'dashed', 'dotted', None,
                   'dashed', 'dotted', None,
                   ],
        color=['black', 'black', 'black',
               'blue', 'blue', 'blue',
               'red', 'red', 'red',
               ],
        xrange=[0, 0.0025],
        title="", yrange=[0,3], figsize=(12,8), show=0)

    if figfile != "":
        plt.savefig(figfile)
        print("File writte to disk: %s" % figfile)

    # plot(
    #      Spread1, Qs1,
    #      Spread3, Qs3,
    #      Spread5, Qs5,
    #      xtitle=r"$\delta_\mathcal{E}$", ytitle="correction for angular width", legend =[
    #                                              "n=1 Qs",
    #                                              "n=3 Qs",
    #                                              "n=5 Qs",
    #                                              ],
    #     title="", yrange=[0,5], figsize=(12,8), show=0)


    plot_show()

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_image, plot_show
    import matplotlib.pylab as plt
    import matplotlib

    matplotlib.rcParams.update({'font.size': 18})


    Spread1, Fwhm1, Sd1, Qa1 = wofry_results(harmonic_number=1,
                                             do_plot=1, # 0=None, 1=image for paper, 2=all images, 3=plots+images,
                                             figfile="spread-map-farfield.png")
    plot_farfield(figfile="spread-farfield.pdf")


    plot_backpropagated(figfile="spread-backpropagated.pdf")

    Spread1, Fwhm1, Sd1, Qa1 = wofry_results_backpropagated(harmonic_number=1,
                                             do_plot=3, # 0=None, 1=image for paper, 2=all images, 3=plots+images,
                                             figfile="spread-map-backpropagated.png")

