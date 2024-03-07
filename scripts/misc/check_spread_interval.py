import numpy
import scipy.constants as codata
from scipy.special import erf

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
    from srxraylib.plot.gol import plot, plot_image
    e0 = 6.0 # GeV


    for i in range(1,7):
        print("0.5%% %dsigma" % i, e0 - e0 * 0.005 * i, e0 + e0 * 0.005 * i)

    #
    #
    #
    npoints = 201
    e_energies = numpy.linspace(5.9, 6.1, npoints)
    spread = numpy.linspace(0, 0.005, npoints)  # (e_energies - e0) / e0  #

    Weights = numpy.zeros((npoints, npoints)) # spread, e energy

    for i in range(npoints):
        Weights[i, :] = numpy.exp(-(e_energies - e0)**2 / 2 / (e0 * spread[i])**2 )

    plot_image(Weights, spread, e_energies, xtitle="sigma spread", ytitle="e energy [GeV]", aspect='auto', show=0)

    #
    # check: recalculate spread
    FWHM = numpy.zeros(npoints)

    for i in range(npoints):
        ys = Weights[i, :]
        # plot(e_energies, ys, title="sigma = %s" % (spread[i]))
        fwhm, _, _ = get_fwhm(ys, e_energies)
        FWHM[i] = fwhm

    plot(e_energies, FWHM / 2.355 / e0,
         e_energies, spread, xtitle="e energy", ytitle="spread", legend=['recalculated from weights', 'original'])

    # Spread1 = numpy.abs(e_energies - e0) / e0
    # plot(e_energies, Spread1, xtitle="e energy [GeV]", ytitle="Spread1")
    # spread = (e_energies - e0) / e0 #  numpy.linspace(0, 0.005, 100)

    Qs1 = numpy.zeros_like(spread)
    Qs5 = numpy.zeros_like(spread)
    Qs7 = numpy.zeros_like(spread)
    Qa1 = numpy.zeros_like(spread)
    Qa5 = numpy.zeros_like(spread)
    Qa7 = numpy.zeros_like(spread)

    for i in range(spread.size):
        Qs1[i] = q_s(norm_energ_spr( numpy.abs(spread[i]), harmonic_number=1, verbose=0))
        Qs5[i] = q_s(norm_energ_spr( numpy.abs(spread[i]), harmonic_number=5, verbose=0))
        Qs7[i] = q_s(norm_energ_spr( numpy.abs(spread[i]), harmonic_number=7, verbose=0))
        Qa1[i] = q_a(norm_energ_spr( numpy.abs(spread[i]), harmonic_number=1, verbose=0))
        Qa5[i] = q_a(norm_energ_spr( numpy.abs(spread[i]), harmonic_number=5, verbose=0))
        Qa7[i] = q_a(norm_energ_spr( numpy.abs(spread[i]), harmonic_number=7, verbose=0))

    plot(spread, Qs1,
         spread, Qs5,
         spread, Qs7,
         spread, Qa1,
         spread, Qa5,
         spread, Qa7,
         legend=['Qs n=1','Qs n=5','Qs n=7','Qa n=1','Qa n=5','Qa n=7'], xtitle="Electron energy spread", ytitle='Qs or Qa')

