import numpy
from scipy.special import erf

from srxraylib.plot.gol import plot
from srxraylib.plot.gol import plot, plot_image, plot_show

import matplotlib.pylab as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 18})


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

if __name__ == "__main__":

    # Load numerical arrays for universal functions
    Fflux = numpy.loadtxt("gwSrwBrilUndHarmUnivFlux.txt")
    Fdiv =  numpy.loadtxt("gwSrwBrilUndHarmUnivDiv.txt")
    Fsize = numpy.loadtxt("gwSrwBrilUndHarmUnivSize.txt")



    nDelta = Fflux.shape[0]
    nepsilon = Fflux.shape[1]
    Delta = numpy.linspace(-10, 10, nDelta)
    epsilon = numpy.linspace(0, 5, nepsilon)
    print(Fflux.shape, Fdiv.shape, Fsize.shape, nDelta, nepsilon )


    if False: plot_image(Fflux, Delta, epsilon, xtitle=r'$\Delta$', ytitle=r'$\varepsilon$', title="Fflux", aspect='auto', show=0)
    if False: plot_image(Fsize, Delta, epsilon, xtitle=r'$\Delta$', ytitle=r'$\varepsilon$', title="Fsize", aspect='auto', show=0)
    if False: plot_image(Fdiv, Delta, epsilon, xtitle=r'$\Delta$', ytitle=r'$\varepsilon$', title="Fdiv", aspect='auto', show=0)

    # scan vs epsilon
    if True:
        Qa = numpy.zeros_like(epsilon)
        Qs = numpy.zeros_like(epsilon)

        Gs0 = Fsize[ nDelta // 2, 0]
        Ga0 = Fdiv [nDelta // 2, 0]

        print("Gs0: ", Gs0)
        print("Ga0: ", Ga0)

        for i, eps in enumerate(epsilon):
            Qa[i] = q_a(2 * numpy.pi * eps)
            Qs[i] = q_s(2 * numpy.pi * eps)

        if True:
            Fs = Fsize[nDelta // 2, :]
            plot(epsilon, Fs / Gs0,
                 epsilon, Qs,
                 xtitle="epsilon", title="cut at Delta=%d" % (Delta[nDelta//2]),
                 legend=['Size Sirepo (norm)', 'Size Tanaka', ], figsize=(9,7), show=0)

            Fa = Fdiv[ nDelta // 2, :]
            plot(
                epsilon, Fa / Ga0,
                epsilon, Qa,
                xtitle="epsilon", title="", # cut at Delta=%d" % (Delta[nDelta // 2, ]),
                legend=['Div Sirepo (norm) Delta=%d eV' % (Delta[nDelta // 2]),
                        'Div Tanaka', ], figsize=(9,7), show=0)


            N = 111.11
            n = 1
            e0 = 10000.0
            E = Delta * e0 / N / n
            spread = epsilon / N / n # / 2 / numpy.pi
            plot(
                spread, Fa / Ga0,
                spread, Qa,
                xtitle="spread", title="", # cut at Delta=%d" % (Delta[nDelta // 2, ]),
                legend=['Div Sirepo (norm) Delta=%d eV' % (Delta[nDelta // 2]),
                        'Div Tanaka', ], figsize=(9,7), show=0)


        if False:
            plot(
                epsilon, Fdiv[0, :]                  , #/ Fdiv[0, 0],
                epsilon, Fdiv[int(nDelta * 0.25), :] , #/ Fdiv[int(nDelta * 0.25), 0],
                epsilon, Fdiv[int(nDelta * 0.50), :] , #/ Fdiv[int(nDelta * 0.50), 0],
                epsilon, Fdiv[int(nDelta * 0.75), :] , #/ Fdiv[int(nDelta * 0.75), 0],
                epsilon, Fdiv[-1, :]                 , #/ Fdiv[-1, 0],
                epsilon, Qa,
                xtitle="epsilon", title="", # cut at Delta=%d" % (Delta[nDelta // 2, ]),
                legend=['Div Sirepo Delta=%d eV' % (Delta[0]),
                        'Div Sirepo Delta=%d eV' % (Delta[int(nDelta * 0.25)] ),
                        'Div Sirepo Delta=%d eV' % (Delta[int(nDelta * 0.50)] ),
                        'Div Sirepo Delta=%d eV' % (Delta[int(nDelta * 0.75)] ),
                        'Div Sirepo Delta=%d eV' % (Delta[-1]),
                        'Div Tanaka', ], show=0)

        FdivN = Fdiv.copy()
        for i, eps in enumerate(epsilon):
            FdivN[:, i] = FdivN[:, i] / Fdiv[:, 0]

        FsizeN = Fsize.copy()
        for i, eps in enumerate(epsilon):
            FsizeN[:, i] = FsizeN[:, i] / Fsize[:, 0]

        if True:
            plot_image(FdivN, Delta, epsilon, xtitle=r'$\Delta$', ytitle=r'$\varepsilon$', title="Fdiv NORMALIZED",
                                aspect='auto', figsize=(9,7), show=0)

            plot_image(FsizeN, Delta, epsilon, xtitle=r'$\Delta$', ytitle=r'$\varepsilon$', title="Fsize NORMALIZED",
                                aspect='auto', figsize=(9,7), show=0)

            N = 111.11
            n = 1
            e0 = 10000.0
            E = Delta * e0 / N / n
            spread = epsilon / N / n # / 2 / numpy.pi

            ibad = numpy.argwhere(spread > 0.005)
            FdivN[:, ibad] = 1
            FsizeN[:, ibad] = 1

            f, _ = plot_image(FdivN, E, spread, xtitle=r'photon energy $E-E_0$ [eV]',
                       ytitle=r'electron energy spread $\delta_\mathcal{E}$', title=r"$F_a$ NORMALIZED",
                       aspect='auto', xrange=[-300, 450], yrange=[0, 0.005], add_colorbar=True, figsize=(9,7), show=0)
            plt.savefig("FdivNash.png")
            print("File written to disk: FdivNash.png")

            plot_image(FsizeN, E, spread, xtitle=r'photon energy $E-E_0$ [eV]',
                       ytitle=r'electron energy spread $\delta_\mathcal{E}$', title=r"$F_s$ NORMALIZED",
                       aspect='auto',  xrange=[-300, 450], yrange=[0, 0.005], figsize=(9,7), show=0)
            plt.savefig("FsizeNash.png")
            print("File written to disk: FsizeNash.png")


        if False:
            plot(
            epsilon, FdivN[0, :],  # / Fdiv[0, 0],
            epsilon, FdivN[int(nDelta * 0.25), :],  # / Fdiv[int(nDelta * 0.25), 0],
            epsilon, FdivN[int(nDelta * 0.50), :],  # / Fdiv[int(nDelta * 0.50), 0],
            epsilon, FdivN[int(nDelta * 0.75), :],  # / Fdiv[int(nDelta * 0.75), 0],
            epsilon, FdivN[-1, :],
            epsilon, Qa,
            xtitle="epsilon", title="",  # cut at Delta=%d" % (Delta[nDelta // 2, ]),
            legend=['DivN Sirepo Delta=%d eV' % (Delta[0]),
                    'DivN Sirepo Delta=%d eV' % (Delta[int(nDelta * 0.25)]),
                    'DivN Sirepo Delta=%d eV' % (Delta[int(nDelta * 0.50)]),
                    'DivN Sirepo Delta=%d eV' % (Delta[int(nDelta * 0.75)]),
                    'DivN Sirepo Delta=%d eV' % (Delta[-1]),
                    'DivN Tanaka', ], show=0)


    # scan vs delta
    if False:
        Gs0 = Fsize[ nDelta // 2, 0]
        Ga0 = Fdiv [nDelta // 2, 0]

        print("Gs0: ", Gs0)
        print("Ga0: ", Ga0)

        Gs = Fsize[:, 0] / Gs0
        Ga =  Fdiv[:, 0] / Ga0
        plot(Delta, Gs,
             Delta, Ga,
             # Delta, Gs * Ga,
             xtitle="Delta", title="cut at epsilon=%f" % (epsilon[0]),
             legend=['Size Sirepo (norm)', 'Div Sirepo (norm)', 'prod Sirepo',], show=0)




    plot_show()
