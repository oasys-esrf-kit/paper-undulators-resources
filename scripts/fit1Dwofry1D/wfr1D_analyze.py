
#
# Import section
#
import numpy


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

if __name__ == "__main__":
    # delete  wfr1D_farfield.h5
    # Shift = numpy.arange(-400, 100, 1)  # RUN1
    # Shift = numpy.arange(-60, 60, 5)
    Shift = numpy.arange(-1000, 1000, 2)
    print(Shift)
    Y = numpy.zeros_like(Shift, dtype=float)
    for i, shift in enumerate(Shift):

        if shift < 0:
            subgroupname = "wfr_%04d" % shift
        else:
            subgroupname = "wfr_%03d" % shift

        output_wavefront = GenericWavefront1D.load_h5_file("wfr1D_farfield_run2.h5", subgroupname)
        x = output_wavefront.get_abscissas()
        y = output_wavefront.get_intensity()
        fwhm, _, _ = get_fwhm(y, x)

        # plot(x, y, title="shift = %d eV FWHM = %f" % (shift, fwhm))

        Y[i] = fwhm

    print(Y)
    plot(Shift, Y, xtitle="Shift [eV]", ytitle='FWHM [eV]', show=0)
    Delta = 1 * 111.111 * (Shift) / 10000.0
    plot(Delta, Y, xtitle="Delta", ytitle='FWHM [eV]')