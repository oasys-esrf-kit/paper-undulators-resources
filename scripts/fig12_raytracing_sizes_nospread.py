def run_redshifted_internal():
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
        emin=9910.0,  # Photon energy scan from energy (in eV)
        emax=9910.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=1.63e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='internal',  # internal, pysru, srw
        flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.2,  # for SRW backpropagation, the range factor
        srw_resolution=100.0,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=1,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.05,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=1,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=1,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.2,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    return beam

def run_redshifted_pysru():
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
        emin=9910.0,  # Photon energy scan from energy (in eV)
        emax=9910.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=1.63e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='pysru',  # internal, pysru, srw
        flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.2,  # for SRW backpropagation, the range factor
        srw_resolution=100.0,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=1,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.01,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=1,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=0,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.2,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    return beam

def run_redshifted_srw():
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
        emin=9910.0,  # Photon energy scan from energy (in eV)
        emax=9910.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=5e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='srw',  # internal, pysru, srw
        flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.05,  # for SRW backpropagation, the range factor
        srw_resolution=100.0,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=0,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.01,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=1,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=0,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.2,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    return beam


def run_resonance_internal():
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
        emin=10000.0,  # Photon energy scan from energy (in eV)
        emax=10000.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=1.63e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='internal',  # internal, pysru, srw
        flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.2,  # for SRW backpropagation, the range factor
        srw_resolution=100.0,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=1,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.05,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=1,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=1,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.2,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    return beam


def run_resonance_pysru():
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
        emin=10000.0,  # Photon energy scan from energy (in eV)
        emax=10000.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=1.63e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='pysru',  # internal, pysru, srw
        flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.2,  # for SRW backpropagation, the range factor
        srw_resolution=100.0,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=1,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.01,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=1,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=0,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.2,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    return beam


def run_resonance_srw():
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
        emin=10000.0,  # Photon energy scan from energy (in eV)
        emax=10000.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=5e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='srw',  # internal, pysru, srw
        flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.05,  # for SRW backpropagation, the range factor
        srw_resolution=100.0,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=0,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.01,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=1,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=0,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.2,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    return beam


def run_blueshifted_internal():
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
        emin=10036.0,  # Photon energy scan from energy (in eV)
        emax=10036.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=1.63e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='internal',  # internal, pysru, srw
        flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.2,  # for SRW backpropagation, the range factor
        srw_resolution=100.0,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=1,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.05,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=1,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=1,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.2,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    return beam


def run_blueshifted_pysru():
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
        emin=10036.0,  # Photon energy scan from energy (in eV)
        emax=10036.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=2.5e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='pysru',  # internal, pysru, srw
        flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.2,  # for SRW backpropagation, the range factor
        srw_resolution=100.0,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=1,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.02,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=1,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=0,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.2,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    return beam


def run_blueshifted_srw():
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
        emin=10036.0,  # Photon energy scan from energy (in eV)
        emax=10036.0,  # Photon energy scan to energy (in eV)
        ng_e=1,  # Photon energy scan number of points
        maxangle=5e-05,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot='srw',  # internal, pysru, srw
        flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance=100.0,  # distance to far field plane
        srw_range=0.05,  # for SRW backpropagation, the range factor
        srw_resolution=100.0,  # for SRW backpropagation, the resolution factor
        srw_semianalytical=0,  # for SRW backpropagation, use semianalytical treatement of phase
        magnification=0.01,  # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source=1,
        # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight=0,  # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio=0.2,  # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread=0,  # for monochromatod sources, apply (1) or not (0) electron energy spread correction
    )

    # light source
    from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam, magnetic_structure=source,
                                          nrays=50000, seed=5676561)
    beam = light_source.get_beam()

    # test plot
    # from srxraylib.plot.gol import plot_scatter
    # rays = beam.get_rays()
    # plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')
    return beam



if __name__ == "__main__":
    from srxraylib.plot.gol import plot, plot_scatter, plot_image, plot_show
    import matplotlib.pylab as plt
    import matplotlib
    matplotlib.rcParams.update({'font.size': 18})

    xrange = [-75e-6, 75e-6]
    yrange = [-75e-6, 75e-6]
    hrange = [-50e-6, 50e-6]
    factor_nrays = 1.0

    flag_emittance = 0
    figure_number = 13 - flag_emittance


    beam = run_redshifted_internal()
    tkt_red_0 = beam.histo2(1,3, nbins=101, xrange=xrange, yrange=yrange, ref=23)
    tkt_red_0_v = beam.histo1(3, nbins=101, xrange=hrange, ref=23)
    beam = run_redshifted_pysru()
    tkt_red_1_v = beam.histo1(3, nbins=101, xrange=hrange, ref=23)
    beam =run_redshifted_srw()
    tkt_red_2_v = beam.histo1(3, nbins=101, xrange=hrange, ref=23)
    #
    beam = run_resonance_internal()
    tkt_zero_0 = beam.histo2(1,3, nbins=101, xrange=xrange, yrange=yrange, ref=23)
    tkt_zero_0_v = beam.histo1(3, nbins=101, xrange=hrange, ref=23)
    beam = run_resonance_pysru()
    tkt_zero_1_v = beam.histo1(3, nbins=101, xrange=hrange, ref=23)
    beam = run_resonance_srw()
    tkt_zero_2_v = beam.histo1(3, nbins=101, xrange=hrange, ref=23)
    #
    beam = run_blueshifted_internal()
    tkt_blue_0 = beam.histo2(1,3, nbins=101, xrange=xrange, yrange=yrange, ref=23)
    tkt_blue_0_v = beam.histo1(3, nbins=101, xrange=hrange, ref=23)
    beam = run_blueshifted_pysru()
    tkt_blue_1_v = beam.histo1(3, nbins=101, xrange=hrange, ref=23)
    beam = run_blueshifted_srw()
    tkt_blue_2_v = beam.histo1(3, nbins=101, xrange=hrange, ref=23)
    #
    # # rays = beam.get_rays()
    # # plot_scatter(1e6 * rays[:, 3], 1e6 * rays[:, 5], title='(X,Z) in microns', show=0)
    #
    #
    #
    #
    plot_image(tkt_red_0["histogram"], 1e6 * tkt_red_0["bin_h_center"], 1e6 * tkt_red_0["bin_v_center"],
               figsize=(9, 8), xtitle=r"X [$\mu$m]", ytitle=r"Y [$\mu$m]", title=r"$E=E_0 (1-(Nn)^{-1})$", show=0)
    plt.savefig("fig%da.png" % figure_number)
    #
    plot_image(tkt_zero_0["histogram"], 1e6 * tkt_zero_0["bin_h_center"], 1e6 * tkt_zero_0["bin_v_center"],
               figsize=(9, 8), xtitle=r"X [$\mu$m]", ytitle=r"Y [$\mu$m]", title=r"$E=E_0$", show=0)
    plt.savefig("fig%db.png" % figure_number)
    #
    plot_image(tkt_blue_0["histogram"], 1e6 * tkt_blue_0["bin_h_center"], 1e6 * tkt_blue_0["bin_v_center"],
               figsize=(9, 8), xtitle=r"X [$\mu$m]", ytitle=r"Y [$\mu$m]", title=r"$E=E_0 (1+0.4 (Nn)^{-1})$", show=0)
    plt.savefig("fig%dc.png" % figure_number)
    #
    #
    plot(1e6 * tkt_red_0_v["bin_path"], tkt_red_0_v["histogram_path"] + 50*factor_nrays,
         1e6 * tkt_red_1_v["bin_path"], tkt_red_1_v["histogram_path"] + 200*factor_nrays,
         1e6 * tkt_red_2_v["bin_path"], tkt_red_2_v["histogram_path"] + 350*factor_nrays,
         legend=[r"internal",r"pySRU",r"SRW"],
         # yrange=[0,850*factor_nrays],
         show=0, #linestyle=['solid','dotted','dashed'],
         figsize=(10,8), xtitle=r"Z [$\mu$m]", ytitle=r"Intensity [arbitrary units]", ) #color=['r','r','r'])
    plt.savefig("fig%dd.pdf" % figure_number)
    #
    plot(1e6 * tkt_zero_0_v["bin_path"], tkt_zero_0_v["histogram_path"] + 50*factor_nrays,
         1e6 * tkt_zero_1_v["bin_path"], tkt_zero_1_v["histogram_path"] + 200*factor_nrays,
         1e6 * tkt_zero_2_v["bin_path"], tkt_zero_2_v["histogram_path"] + 350*factor_nrays,
         legend=[r"internal",r"pySRU",r"SRW"],
         # yrange=[0,1250*factor_nrays],
         show=0, #linestyle=['solid','dotted','dashed'],
         figsize=(10,8), xtitle=r"Z [$\mu$m]", ytitle=r"Intensity [arbitrary units]", ) # color=['k','k','k'])
    plt.savefig("fig%de.pdf" % figure_number)
    #
    plot(1e6 * tkt_blue_0_v["bin_path"], tkt_blue_0_v["histogram_path"] + 50*factor_nrays,
         1e6 * tkt_blue_1_v["bin_path"], tkt_blue_1_v["histogram_path"] + 200*factor_nrays,
         1e6 * tkt_blue_2_v["bin_path"], tkt_blue_2_v["histogram_path"] + 350*factor_nrays,
         legend=[r"internal",r"pySRU",r"SRW"],
         # yrange=[0,1450*factor_nrays],
         show=0, #linestyle=['solid','dotted','dashed'],
         figsize=(10,8), xtitle=r"Z [$\mu$m]", ytitle=r"Intensity [arbitrary units]", ) # color=['b','b','b'])
    plt.savefig("fig%df.pdf" % figure_number)

    plot_show()

