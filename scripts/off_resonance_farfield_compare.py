
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

def get_stdev(x, y):
    delta = (x[1] - x[0])
    Y = y.copy()
    Y /= y.sum() * delta
    m1 = (x ** 1 * Y).sum() * delta
    m2 = (x ** 2 * Y).sum() * delta
    return numpy.sqrt(m2 - m1**2)

if __name__ == "__main__":

    from srxraylib.plot.gol import plot, plot_image, plot_show

    case = 2
    harmonic_number = 1


    #
    # Wofry results
    #

    if True:
        if case == 1:
            Shift = numpy.arange(-400, 100, 1) # case 1
            filename = "fit1Dwofry1D/wfr1D_farfield_n%d_run1.h5" % harmonic_number
        elif case == 2:
            Shift = numpy.arange(-1000, 1000, 2) # case 2
            filename = "fit1Dwofry1D/wfr1D_farfield_n%d_run2.h5" % harmonic_number

        Y = numpy.zeros_like(Shift, dtype=float)
        SD = numpy.zeros_like(Shift, dtype=float)
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
            SD[i] = get_stdev(x, y)
            if i == 0:
                YW_IMG = numpy.zeros((Shift.size, y.size))
            YW_IMG[i, :] = y

        plot_image(YW_IMG, Shift, x * 1e3,
                   yrange=[-2, 2],
                   xtitle=r"Photon energy $E-E_0$ [eV]", ytitle="x [mm] @ far field (100m)", title="", aspect='auto', show=0)

        # plot(Shift, Y,
        #      Shift, SD*2.355,
        #      xtitle="Shift [eV]", ytitle='WOFRY FWHM % 100m [m]', legend=['FWHM','SD*2.355'], show=0)

        Delta = 1 * 111.111 * (Shift) / 10000.0

        # plot(Delta, Y,
        #      Delta, SD*2.355,
        #      xtitle="Delta", ytitle='WOFRY FWHM % 100m [m]', legend=['FWHM','SD*2.355'], show=0)

        ShiftW = Shift.copy()
        DeltaW = Delta.copy()
        YW = Y.copy()
        SDW = SD.copy()

        #
        #  check resonance shift-> e energy shift
        #
        if False:
            import scipy.constants as codata

            e = 10000.0 + Shift
            w1 = codata.c * codata.h / codata.e / e
            K = 1.341095
            lambda_u = 0.018
            gamma1 = numpy.sqrt((1 + K ** 2 / 2) * (lambda_u / 2 / w1))  # DG + gamma0
            Ee1 = gamma1 * 0.51099895e-3
            plot(gamma1, Ee1, show=0)


            plot_image(YW_IMG, Ee1, x * 1e3,
                       yrange=[-2, 2],
                       xtitle=r"CONVERTED e energy [GeV]", ytitle="x [mm] @ far field (100m)", title="", aspect='auto', show=0)

            YW_IMGi = YW_IMG.copy()
            Eei = numpy.linspace(Ee1.min(), Ee1.max(), Ee1.size)
            for i in range(YW_IMG.shape[1]):
                tmp = numpy.interp(Eei, Ee1, YW_IMG[:, i])
                YW_IMGi[:, i] = tmp

            plot_image(YW_IMGi, Eei, x * 1e3,
                       yrange=[-2, 2],
                       xtitle=r"INTERPOLATED e energy [GeV]", ytitle="x [mm] @ far field (100m)", title="", aspect='auto', show=0)

    #
    # spectra results
    #

    if True: # and harmonic_number < 2:

        datadir = "/scisoft/data/srio/paper-undulator/spectra_results/"
        if case == 1:
            Shift = numpy.arange(-400, 100, 1) # case 1
            filename = datadir + "Spectra_Divergence_Off_Res_ESRF_ID06_EBS_CPMU18_1_0_emitt_0_spread_case1.hdf5"
        elif case == 2:
            Shift = numpy.arange(-1000, 1000, 2) # case 2
            filename = datadir + "Spectra_Divergence_Off_Res_ESRF_ID06_EBS_CPMU18_1_0_emitt_0_spread_case2.hdf5"

        import h5py

        file = h5py.File(filename, 'r')
        shift =     file['/Scan/harmonic 01/shift'][()]
        flux_dens = file["/Scan/harmonic 01/flux_dens"][()]
        x_div =         file["/Scan/harmonic 01/x_div"][()]
        y_div =         file["/Scan/harmonic 01/y_div"][()]

        print(shift.shape, flux_dens.shape, x_div.shape, y_div.shape)

        nx = flux_dens.shape[1]

        # plot(y_div, flux_dens[0, nx // 2, :], xtitle="Angle in mrad")

        YS = numpy.zeros_like(Shift, dtype=float)
        SD = numpy.zeros_like(Shift, dtype=float)
        for i, shift in enumerate(Shift):
            ys = flux_dens[i, nx // 2, :]
            fwhm, _, _ = get_fwhm(ys, y_div)
            # plot(x_div, ys, title="shift = %d eV FWHM = %f" % (shift, fwhm), xtitle="Angle in mrad")
            YS[i] = fwhm * 1e-3 * 100 # in m @ 100m
            SD[i] = get_stdev(y_div * 1e-3 * 100, ys)
            if i == 0:
                Y_IMG = numpy.zeros((Shift.size, y_div.size))
            Y_IMG[i, :] = ys

        plot_image(Y_IMG, Shift, y_div * 1e-3 * 100,
                   yrange=[-0.002, 0.002],
                   xtitle="Shift [eV]", ytitle="Far field x[mm]", title="SPECTRA", aspect='auto', show=0)

        # plot(Shift, YS,
        #      Shift, SD * 2.355,
        #      xtitle="Shift [eV]", ytitle='SPECTRA FWHM % 100m [m]', legend=['FWHM','SD*2.355'], show=0)

        Delta = 1 * 111.111 * (Shift) / 10000.0

        # plot(Delta, YS,
        #      Delta, SD * 2.355,
        #      xtitle="Delta", ytitle='SPECTRA FWHM % 100m [m]', legend=['FWHM','SD*2.355'], show=0)

        plot(Delta, YS,
             DeltaW, YW,
             Delta,  SD,
             DeltaW, SDW,
             xtitle="Delta=N n DE/Eo", ytitle="width at far field (100 m) [m]",
             legend=['FWHM (SPECTRA)','FWHM (WOFRY)','SD*2.355 (SPECTRA)','SD*2.355 (WOFRY)'],
             xrange=[-2,5], yrange=[0,0.005], show=0)

        plot(Shift , YS,
             ShiftW, YW,
             Shift ,  SD,
             ShiftW, SDW,
             xtitle="DE [eV]", ytitle="width at far field (100 m) [m]",
             legend=['FWHM (SPECTRA)','FWHM (WOFRY)','SD*2.355 (SPECTRA)','SD*2.355 (WOFRY)'],
             xrange=[-300,500], yrange=[0,0.005], show=0)


    plot_show()