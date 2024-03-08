import numpy
from scipy.special import erf

from srxraylib.plot.gol import plot
from srxraylib.plot.gol import plot, plot_image, plot_show



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

    if False:
        plot_image(Fflux, Delta, epsilon, xtitle=r'$\Delta$', ytitle=r'$\varepsilon$', title="Fflux", aspect='auto', show=0)
        plot_image(Fsize, Delta, epsilon, xtitle=r'$\Delta$', ytitle=r'$\varepsilon$', title="Fsize", aspect='auto', show=0)
        plot_image(Fdiv, Delta, epsilon, xtitle=r'$\Delta$', ytitle=r'$\varepsilon$', title="Fdiv", aspect='auto', show=0)

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

        Fs = Fsize[nDelta // 2, :]
        Fa = Fdiv[ nDelta // 2, :]
        plot(epsilon, Fs / Gs0,
             epsilon, Fa / Ga0,
             epsilon, Qs,
             epsilon, Qa,
             epsilon, Fs / Gs0 * Fa / Ga0,
             epsilon, Qa * Qs,
             xtitle="epsilon", title="cut at Delta=%d" % (Delta[nDelta//2]),
             legend=['Size Sirepo (norm)', 'Div Sirepo (norm)', 'Size Tanaka', 'Div Tanaka', 'prod Sirepo', 'prod Tanaka'], show=0)


    # scan vs delta
    if True:


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
