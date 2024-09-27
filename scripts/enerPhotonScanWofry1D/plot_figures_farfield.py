
#
# Import section
#
import numpy

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D
from srxraylib.plot.gol import plot

import matplotlib.pylab as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 18})

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

if __name__ == "__main__":

    from srxraylib.plot.gol import plot, plot_image, plot_show

    import matplotlib.pylab as plt
    import matplotlib
    matplotlib.rcParams.update({'font.size': 18})


    harmonic_number = 1


    #
    # Wofry results
    #

    Shift = Shift = numpy.linspace(-300, 501, 201)
    filename = "wfr1D_photon_energy_scan_farfield_n%d.h5" % harmonic_number

    Y = numpy.zeros_like(Shift, dtype=float)
    SD = numpy.zeros_like(Shift, dtype=float)
    INTENSITY = numpy.zeros_like(Shift, dtype=float)
    for i, shift in enumerate(Shift):
        subgroupname = "wfr_%5.3f" % (10000.0 * harmonic_number + shift)
        output_wavefront = GenericWavefront1D.load_h5_file(filename, subgroupname)
        x = output_wavefront.get_abscissas()
        y = output_wavefront.get_intensity()
        fwhm, _, _ = get_fwhm(y, x)

        Y[i] = fwhm
        SD[i] = get_stdev(x, y)
        INTENSITY[i] = y.sum()
        if i == 0:
            YW_IMG = numpy.zeros((Shift.size, y.size))
        YW_IMG[i, :] = y

    plot_image(YW_IMG, Shift, x * 1e3,
               yrange=[-2, 2],
               xtitle=r"Photon energy $E-E_0$ [eV]", ytitle="y [mm] @ far field (100m)", title="",
               aspect='auto', figsize=(9,7), show=0)
    plt.savefig("detuning-map-farfield.png")
    print("File written to disk: detuning-map-resonance.png")


    Delta = 1 * 111.111 * (Shift) / 10000.0
    ShiftW = Shift.copy()
    DeltaW = Delta.copy()
    YW = Y.copy()
    SDW = SD.copy()

    plot(ShiftW, YW  * 1e3,
         ShiftW, SDW * 1e3,
         ShiftW, 5 * INTENSITY / INTENSITY.max(),
         xtitle=r"Photon energy $E-E_0$ [eV]", ytitle="",
         legend=[r'FWHM [mm]',r'SD [mm]', 'Intensity [a.u.]'],
         xrange=[-300,500],
         # yrange=[0,5],
         figsize=(12,8), show=0)
    plt.savefig("detuning-farfield.pdf")
    print("File written to disk: angular_width_vs_detuning.pdf")

    # plot(ShiftW, INTENSITY,
    #      xtitle=r"Photon energy $E-E_0$ [eV]", ytitle="Intensity [arbitrary units]",
    #      xrange=[-300,500], figsize=(12,8), show=0)

    plot_show()