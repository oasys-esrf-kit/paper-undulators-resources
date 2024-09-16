

def run_shadow4(energy=10000, code_undul_phot='internal'):
    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam
    electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
    electron_beam.set_sigmas_all(sigma_x=3.01836e-05, sigma_y=3.63641e-06, sigma_xp=4.36821e-06, sigma_yp=1.37498e-06)

    # magnetic structure
    from shadow4.sources.undulator.s4_undulator import S4Undulator
    source = S4Undulator(
        K_vertical=1.341095,  # syned Undulator parameter
        period_length=0.018,  # syned Undulator parameter
        number_of_periods=111.111,  # syned Undulator parameter
        emin=energy,  # Photon energy scan from energy (in eV)
        emax=energy,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=1.63e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot=code_undul_phot,  # internal, pysru, srw
        flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=0,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.05,  # for SRW backpropagation, the range factor
        srw_resolution=50,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=0,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.05,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=0,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=0,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.5,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=25000, seed=5676561)
    beam = light_source.get_beam()

    # # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')

    return beam

if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_scatter, plot_image, plot_show
    import matplotlib.pylab as plt
    import matplotlib
    matplotlib.rcParams.update({'font.size': 18})

    xrange = [-20e-6, 20e-6]
    beam = run_shadow4(energy=10000.0, code_undul_phot='internal')
    tkt_zero_0 = beam.histo2(3,6, nbins=101, xrange=xrange, yrange=xrange, ref=23)
    tkt_zero_0_v = beam.histo1(6, nbins=101, xrange=xrange, ref=23)

    beam = run_shadow4(energy=9909.91, code_undul_phot='internal')
    tkt_red_0 = beam.histo2(3,6, nbins=101, xrange=xrange, yrange=xrange, ref=23)
    tkt_red_0_v = beam.histo1(6, nbins=101, xrange=xrange, ref=23)

    beam = run_shadow4(energy=10036.04, code_undul_phot='internal')
    tkt_blue_0 = beam.histo2(3,6, nbins=101, xrange=xrange, yrange=xrange, ref=23)
    tkt_blue_0_v = beam.histo1(6, nbins=101, xrange=xrange, ref=23)

    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 3], 1e6 * rays[:, 5], title='(X,Z) in microns', show=0)


    plot_image(tkt_zero_0["histogram"], 1e6 * tkt_zero_0["bin_h_center"], 1e6 * tkt_zero_0["bin_v_center"],
               figsize=(9, 8), xtitle=r"X [$\mu$m]", ytitle=r"Y [$\mu$m]", title='', show=0)

    plot_image(tkt_red_0["histogram"], 1e6 * tkt_red_0["bin_h_center"], 1e6 * tkt_red_0["bin_v_center"],
               figsize=(9, 8), xtitle=r"X [$\mu$m]", ytitle=r"Y [$\mu$m]", title='', show=0)

    plot_image(tkt_blue_0["histogram"], 1e6 * tkt_blue_0["bin_h_center"], 1e6 * tkt_blue_0["bin_v_center"],
               figsize=(9, 8), xtitle=r"X [$\mu$m]", ytitle=r"Y [$\mu$m]", title='', show=0)


    plot(1e6 * tkt_zero_0_v["bin_path"], tkt_zero_0_v["histogram_path"],
         1e6 * tkt_red_0_v["bin_path"], tkt_red_0_v["histogram_path"],
         1e6 * tkt_blue_0_v["bin_path"], tkt_blue_0_v["histogram_path"],
         legend=[r"$E_0$",r"$E_0 (1-1/(Nn))$",r"$E_0 (1+0.4/(Nn))$"],
         figsize=(10,8), xtitle=r"Z [$\mu$m]", ytitle=r"Intensity [arbitrary units]", color=['k','r','b'])

    plot_show()

