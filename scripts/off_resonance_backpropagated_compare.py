
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

    case = 1


    #
    # Wofry results
    #

    if False:
        if case == 1:
            Shift = numpy.arange(-400, 100, 1) # case 1
            filename = "fit1Dwofry1D/wfr1D_farfield_run1.h5"
        elif case == 2:
            Shift = numpy.arange(-1000, 1000, 2) # case 2
            filename = "fit1Dwofry1D/wfr1D_farfield_run2.h5"

        Y = numpy.zeros_like(Shift, dtype=float)
        for i, shift in enumerate(Shift):

            if shift < 0:
                subgroupname = "wfr_%04d" % shift
            else:
                subgroupname = "wfr_%03d" % shift

            output_wavefront = GenericWavefront1D.load_h5_file(filename, subgroupname)
            x = output_wavefront.get_abscissas()
            y = output_wavefront.get_intensity()
            fwhm, _, _ = get_fwhm(y, x)

            # plot(x, y, title="shift = %d eV FWHM = %f" % (shift, fwhm))

            Y[i] = fwhm

        plot(Shift, Y, xtitle="Shift [eV]", ytitle='FWHM [eV]', show=0)
        Delta = 1 * 111.111 * (Shift) / 10000.0
        plot(Delta, Y, xtitle="Delta", ytitle='FWHM [eV]')

    #
    # spectra results
    #

    if True:
        if case == 1:
            Shift = numpy.arange(-400, 100, 1) # case 1
            filename = "/tmp_14_days/reyesher/to_Manolo/spectra_results/Spectra_charac_at_source_point_Off_Res_ESRF_ID06_EBS_CPMU18_1_0_emitt_0_spread_case1.hdf5"
            filename = "/tmp_14_days/reyesher/to_Manolo/spectra_results/Spectra_charac_at_source_point_Off_Res_ESRF_ID06_EBS_CPMU18_1_0_emitt_0_spread_case1.hdf5"
        elif case == 2:
            Shift = numpy.arange(-1000, 1000, 2) # case 2
            filename = "/tmp_14_days/reyesher/to_Manolo/spectra_results/Spectra_charac_at_source_point_Off_Res_ESRF_ID06_EBS_CPMU18_1_0_emitt_0_spread_case2.hdf5"

        import h5py

        file = h5py.File(filename, 'r')
        shift =     file['/Scan/harmonic 01/shift'][()]
        flux_dens = file["/Scan/harmonic 01/flux_dens"][()]
        x =         file["/Scan/harmonic 01/x"][()]
        y =         file["/Scan/harmonic 01/y"][()]

        print(shift.shape, flux_dens.shape, x.shape, y.shape)

        nx = flux_dens.shape[1]

        plot(y, flux_dens[0, nx // 2, :])

        YS = numpy.zeros_like(Shift, dtype=float)
        for i, shift in enumerate(Shift):
            ys = flux_dens[i, nx // 2, :]
            fwhm, _, _ = get_fwhm(ys, y)
            # plot(y, ys, title="shift = %d eV FWHM = %f" % (shift, fwhm))
            YS[i] = fwhm

        plot(Shift, YS, xtitle="Shift [eV]", ytitle='SPECTRA FWHM [mm]', show=0)
        Delta = 1 * 111.111 * (Shift) / 10000.0
        plot(Delta, YS, xtitle="Delta", ytitle='SPECTRA FWHM [mm]')