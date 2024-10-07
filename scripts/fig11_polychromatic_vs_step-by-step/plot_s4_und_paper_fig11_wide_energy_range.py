import numpy
from syned.util.json_tools import load_from_json_file
from shadow4.sources.s4_electron_beam import S4ElectronBeam
from shadow4.sources.undulator.s4_undulator import S4Undulator
from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
import matplotlib.pyplot as plt
import h5py
from srxraylib.util.h5_simple_writer import H5SimpleWriter

""" Here we have a set of functions to calculate the source size and divergence for:
      
    - A polychromatic source for range energy min to energy max and save it in a HDF5

    - Calculates a monochromatic source for a each step in a range of energies
      and the performs a weighted sum and save it in a HDF5

    And then there is a function that reads the HDF5 and plots a comparison

In our S4 undulator paper (ver 04-Oct-24) it corresponds to figure 11.

"""


def run_source_full_harmonic(syned_json_file, res_photon_energy=10e3, energy_min= 9600,
              energy_max=10200, energy_points=11, harmonic=1, code_undul_phot = 'internal',
              emittance = True, nrays=100000, seed=5676561, flag_energy_spread=0,
              energy_spread=0.001, nbins=201, save_file=False):
    
    """ Using the polychromatic option in the S4UndulatorLightSource, this
    function calculates, at the source position, the spatial and divergence ray beam
    profiles and save them in a HDF5.

    Calcualtions are done for the photon energy range: energy_min to energy_max 
    and a given energy_points.

    """
    # electron beam
        
    syned_obj = load_from_json_file(syned_json_file)
    e_beam = syned_obj.get_electron_beam()
    undulator = syned_obj.get_magnetic_structure()    

    k_value = undulator.get_K_from_photon_energy(res_photon_energy,
                                                e_beam.gamma(),
                                                harmonic=harmonic)
    
    maxangle = 3 * 0.69 * (1/e_beam.gamma())*numpy.sqrt((1.0/(2.0*harmonic * undulator.number_of_periods())) * (1.0 + k_value**2/2.0))
    print(f'For an energy of {res_photon_energy} and harmonic {harmonic} we have a K of {k_value}')

    electron_beam = S4ElectronBeam(energy_in_GeV=e_beam.energy(), energy_spread=energy_spread, current=e_beam.current())

    if code_undul_phot == 'internal':
        magnification     = 0.0001
        flag_backprop_recalculate_source = 1
        flag_backprop_weight = 1
        weight_ratio         = 0.2
    elif code_undul_phot == 'pysru':
        magnification     = 0.01
        flag_backprop_recalculate_source = 1
        flag_backprop_weight = 1
        weight_ratio         = 0.5
    elif code_undul_phot == 'srw':
        magnification     = 0.01
        flag_backprop_recalculate_source = 0
        flag_backprop_weight = 1
        weight_ratio         = 0.5
    else:
        raise RuntimeError("ERROR: Unidentified calculation approach - plotting failed") 

    # electron beam   
    
    if emittance:
        flag_emittance = 1             
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = e_beam.get_sigmas_all()
    else:
        flag_emittance = 0
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = (0.0, 0.0, 0.0, 0.0) 
    
    electron_beam.set_sigmas_all(sigma_x=sigma_x,sigma_y=sigma_y,sigma_xp=sig_div_x,sigma_yp=sig_div_y)
    
    # magnetic structure
    
    source = S4Undulator(
        K_vertical        = k_value, # syned Undulator parameter
        period_length     = undulator.period_length(), # syned Undulator parameter
        number_of_periods = undulator.number_of_periods(), # syned Undulator parameter
        emin              = energy_min, # Photon energy scan from energy (in eV)
        emax              = energy_max, # Photon energy scan to energy (in eV)
        ng_e              = energy_points, # Photon energy scan number of points
        maxangle          = maxangle, # Maximum radiation semiaperture in RADIANS
        ng_t              = 100, # Number of points in angle theta
        ng_p              = 21, # Number of points in angle phi
        ng_j              = 20, # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot   = code_undul_phot, # internal, pysru, srw
        flag_emittance    = flag_emittance, # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size         = 2, # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        distance          = 100.0, # distance to far field plane
        srw_range         = 0.05, # for SRW backpropagation, the range factor
        srw_resolution    = 50, # for SRW backpropagation, the resolution factor
        srw_semianalytical= 0, # for SRW backpropagation, use semianalytical treatement of phase
        magnification     = magnification, # for internal/wofry backpropagation, the magnification factor
        flag_backprop_recalculate_source      = flag_backprop_recalculate_source, # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
        flag_backprop_weight = flag_backprop_weight, # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
        weight_ratio         = weight_ratio, # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
        flag_energy_spread   = flag_energy_spread, # for monochromatod sources, apply (1) or not (0) electron energy spread correction
        )
    
    
    # light source    
    light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam,
                                         magnetic_structure=source, nrays=nrays, seed=seed)
    beam = light_source.get_beam()    

    #we get the values 
    #spatial
    tkt_spatial = beam.histo2(1, 3, ref=23, xrange=[-120e-6, 120e-6], yrange=[-12e-6, 12e-6], nbins=nbins, nolost=1)
    x = tkt_spatial['bin_h_center']
    y = tkt_spatial['bin_v_center']
    histo_spatial = tkt_spatial['histogram']
    #divergence
    tkt_diver = beam.histo2(4, 6, ref=23, xrange=[-30e-6, 30e-6], yrange=[-20e-6, 20e-6], nbins=nbins, nolost=1)
    x_div = tkt_diver['bin_h_center']
    y_div = tkt_diver['bin_v_center']
    histo_div = tkt_diver['histogram']
    
    if save_file:
        
        h5w = H5SimpleWriter.initialize_file(f"s4_full_harmonic_profiles_{code_undul_phot}.h5", creator="h5_basic_writer.py")

        h5w.add_image(histo_spatial, x * 1e6, y * 1e6,
                      image_name="s4_full_harm_profile_spatial",
                      title_x="h [um]",title_y="v [um]")
        
        h5w.add_image(histo_div, x_div * 1e6, y_div * 1e6,
                      image_name="s4_full_harm_profile_div",
                      title_x="h_div [urad]",title_y="v_div [urad]")       
        
        print(f's4_full_harmonic_profiles_{code_undul_phot}.h5 has been save to disk')
        
def run_source_by_energy_step(syned_json_file, res_photon_energy=10e3, energy_min= 9600,
              energy_max=10200, energy_points=11, harmonic=1, code_undul_phot = 'internal',
              emittance = True, nrays=100000, seed=5676561, flag_energy_spread=0,
              energy_spread=0.001, nbins=201, save_file = True):
    
    """ Using the monochromatic option in the S4UndulatorLightSource, this
    function calculates, at the source position, the spatial and divergence ray beam
    profiles for a given energy_point between energy_min and energy_max. 

    Using the flux, it performs a weighted sum over the different energies and
    save the sum in a HDF5.      
    """
    
        # electron beam
        
    syned_obj = load_from_json_file(syned_json_file)
    e_beam = syned_obj.get_electron_beam()
    undulator = syned_obj.get_magnetic_structure()    

    k_value = undulator.get_K_from_photon_energy(res_photon_energy, e_beam.gamma(),
                                                harmonic=harmonic)
    
    maxangle = 3 * 0.69 * (1/e_beam.gamma())*numpy.sqrt((1.0/(2.0*harmonic * undulator.number_of_periods())) * (1.0 + k_value**2/2.0))
    print(f'For an energy of {res_photon_energy} and harmonic {harmonic} we have a K of {k_value}')

    electron_beam = S4ElectronBeam(energy_in_GeV=e_beam.energy(), energy_spread=energy_spread, current=e_beam.current())

    if code_undul_phot == 'internal':
        magnification     = 0.0001
        flag_backprop_recalculate_source = 1
        flag_backprop_weight = 1
        weight_ratio         = 0.2
    elif code_undul_phot == 'pysru':
        magnification     = 0.01
        flag_backprop_recalculate_source = 1
        flag_backprop_weight = 1
        weight_ratio         = 0.5
    elif code_undul_phot == 'srw':
        magnification     = 0.01
        flag_backprop_recalculate_source = 0
        flag_backprop_weight = 1
        weight_ratio         = 0.5
    else:
        raise RuntimeError("ERROR: Unidentified calculation approach - plotting failed") 

    # electron beam   
    
    if emittance:
        flag_emittance = 1             
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = e_beam.get_sigmas_all()
    else:
        flag_emittance = 0
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = (0.0, 0.0, 0.0, 0.0) 
    
    electron_beam.set_sigmas_all(sigma_x=sigma_x,sigma_y=sigma_y,sigma_xp=sig_div_x,sigma_yp=sig_div_y)
    
    # we creat the energy steps
    
    photon_energies = numpy.linspace(energy_min, energy_max, energy_points)
    
    # empty array to be filled with the calculations
    flux = numpy.zeros([energy_points])
    profile_spatial = numpy.zeros([nbins, nbins, energy_points])
    profile_div = numpy.zeros([nbins, nbins, energy_points])
    
    for i, photon_energy in enumerate(photon_energies):

        print(f'Calculating for energy {photon_energy} step {i+1} of {len(photon_energies)}')

        source = S4Undulator(
            K_vertical        = k_value, # syned Undulator parameter
            period_length     = undulator.period_length(), # syned Undulator parameter
            number_of_periods = undulator.number_of_periods(), # syned Undulator parameter
            emin              = photon_energy, # Photon energy scan from energy (in eV)
            emax              = photon_energy, # Photon energy scan to energy (in eV)
            ng_e              = energy_points, # Photon energy scan number of points
            maxangle          = maxangle, # Maximum radiation semiaperture in RADIANS
            ng_t              = 100, # Number of points in angle theta
            ng_p              = 21, # Number of points in angle phi
            ng_j              = 20, # Number of points in electron trajectory (per period) for internal calculation only
            code_undul_phot   = code_undul_phot, # internal, pysru, srw
            flag_emittance    = flag_emittance, # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_size         = 2, # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
            distance          = 100.0, # distance to far field plane
            srw_range         = 0.05, # for SRW backpropagation, the range factor
            srw_resolution    = 50, # for SRW backpropagation, the resolution factor
            srw_semianalytical= 0, # for SRW backpropagation, use semianalytical treatement of phase
            magnification     = magnification, # for internal/wofry backpropagation, the magnification factor
            flag_backprop_recalculate_source      = flag_backprop_recalculate_source, # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
            flag_backprop_weight = flag_backprop_weight, # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
            weight_ratio         = weight_ratio, # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
            flag_energy_spread   = flag_energy_spread, # for monochromatod sources, apply (1) or not (0) electron energy spread correction
            )        
        
        # light source    
        light_source = S4UndulatorLightSource(name='undulator', electron_beam=electron_beam,
                                             magnetic_structure=source, nrays=nrays, seed=seed)
        beam = light_source.get_beam()       
        # we get the all values
        flux[i] += light_source.get_flux_and_spectral_power()[0]   
        # spatial
        tkt_spatial = beam.histo2(1, 3, ref=23, xrange=[-120e-6, 120e-6], yrange=[-12e-6, 12e-6], nbins=201, nolost=1)
        x = tkt_spatial['bin_h_center']
        y = tkt_spatial['bin_v_center']
        histo_spatial = tkt_spatial['histogram']

        # divergence
        tkt_div = beam.histo2(4, 6, ref=23, xrange=[-30e-6, 30e-6], yrange=[-20e-6, 20e-6], nbins=201, nolost=1)                   
        x_div = tkt_div['bin_h_center']
        y_div = tkt_div['bin_v_center']
        histo_div = tkt_div['histogram']
        
        # we weigth each one by its flux
        profile_spatial[:, :, i] += histo_spatial * flux[i]
        profile_div[:, :, i] += histo_div * flux[i] 

    #we devide by the maxium flux, normalizing
    profile_spatial /= max(flux)
    profile_div /= max(flux)

    # to save the files we sum each energy step profile over the energy axis
    # and in microns and micro-rads:

    if save_file:

        h5w = H5SimpleWriter.initialize_file(f"s4_e_step_full_harmonic_profiles_{code_undul_phot}.h5", creator="h5_basic_writer.py")

        h5w.add_image(numpy.sum(profile_spatial, axis=2), x * 1e6, y * 1e6,
                      image_name="s4_full_harm_profile_spatial",
                      title_x="h [um]",title_y="v [um]")
        
        h5w.add_image(numpy.sum(profile_div, axis=2), x_div * 1e6, y_div * 1e6,
                      image_name="s4_full_harm_profile_div",
                      title_x="h_div [urad]",title_y="v_div [urad]")
        
        h5w.add_dataset(photon_energies, flux, dataset_name='Photon flux',
                        title_x='Photon energy [eV]', title_y='Photon flux [Phot/sec/0.1%bw]')
        
        print(f's4_e_step_full_harmonic_profiles_{code_undul_phot}.h5 has been save to disk')

def compare_results(h5_file_s4_full, h5_file_s4_steps, f_size=12, line_width=4):

    """ Here, by reading each created HDF5 and then we sum rows (for hortizontal)
    and columns (for vertical) for each spatial and divergence, and then it
    plots the comparison.
    """
    
    ### Reading both files ###

    h5_full = h5py.File(h5_file_s4_full, 'r')
    h5_steps = h5py.File(h5_file_s4_steps, 'r')    
    
    full_spatial_x = numpy.copy(h5_full['s4_full_harm_profile_spatial/axis_x'])
    full_spatial_y = numpy.copy(h5_full['s4_full_harm_profile_spatial/axis_y'])
    full_spatial_data = numpy.copy(h5_full['s4_full_harm_profile_spatial/image_data'])

    full_div_x = numpy.copy(h5_full['s4_full_harm_profile_div/axis_x']) 
    full_div_y = numpy.copy(h5_full['s4_full_harm_profile_div/axis_y'])
    full_div_data = numpy.copy(h5_full['s4_full_harm_profile_div/image_data'])

    steps_spatial_x = numpy.copy(h5_steps['s4_full_harm_profile_spatial/axis_x'])
    steps_spatial_y = numpy.copy(h5_steps['s4_full_harm_profile_spatial/axis_y'])
    steps_spatial_data = numpy.copy(h5_steps['s4_full_harm_profile_spatial/image_data'])

    steps_div_x = numpy.copy(h5_steps['s4_full_harm_profile_div/axis_x'])
    steps_div_y = numpy.copy(h5_steps['s4_full_harm_profile_div/axis_y'])
    steps_div_data = numpy.copy(h5_steps['s4_full_harm_profile_div/image_data'])

    h5_full.close()
    h5_steps.close()    

    ### HORIZONTAL SIZE ###

    #sum all the rows and normalize to the peak
    full_spatial_horizontal = numpy.sum(full_spatial_data, axis=0) / max(numpy.sum(full_spatial_data, axis=0))
    steps_spatial_horizontal = numpy.sum(steps_spatial_data, axis=0) / max(numpy.sum(steps_spatial_data, axis=0))

    plt.plot(full_spatial_x, full_spatial_horizontal, label='Polychromatic', lw=line_width)
    plt.plot(steps_spatial_x, steps_spatial_horizontal, label='Energy steps', lw=line_width)

    plt.xlabel("Horizontal size ($\mu$m)", fontsize= f_size)    
    plt.ylabel("Normalized intensity (a.u.)", fontsize= f_size)

    ##plt.xlim(4e3, 21e4)
    ##plt.ylim(1e12, 5e15)
    plt.xticks(fontsize= f_size)
    plt.yticks(fontsize= f_size)    
    ##plt.grid(which='both', axis='y')   
    ##plt.title(f"Method: {label}, value: {plot_value}", fontsize=f_size)
    plt.legend(fontsize = f_size)
    plt.show()

    ### VERTICAL SIZE ###

    #sum all the columns and normalize to the peak
    full_spatial_vertical = numpy.sum(full_spatial_data, axis=1) / max(numpy.sum(full_spatial_data, axis=1))
    steps_spatial_vertical = numpy.sum(steps_spatial_data, axis=1) / max(numpy.sum(steps_spatial_data, axis=1))

    plt.plot(full_spatial_y, full_spatial_vertical, label='Polychromatic', lw=line_width)
    plt.plot(steps_spatial_y, steps_spatial_vertical, label='Energy steps', lw=line_width)

    plt.xlabel("Vertical size ($\mu$m)", fontsize= f_size)    
    plt.ylabel("Normalized intensity (a.u.)", fontsize= f_size)

    ##plt.xlim(4e3, 21e4)
    ##plt.ylim(1e12, 5e15)
    plt.xticks(fontsize= f_size)
    plt.yticks(fontsize= f_size)    
    ##plt.grid(which='both', axis='y')   
    ##plt.title(f"Method: {label}, value: {plot_value}", fontsize=f_size)
    plt.legend(fontsize = f_size)
    plt.show()
   
    ##### HORIZONTAL DIVERGENCE #####

    #sum all the rows and normalize to the peak
    full_div_hor = numpy.sum(full_div_data, axis=0) / max(numpy.sum(full_div_data, axis=0))
    steps_div_hor = numpy.sum(steps_div_data, axis=0) / max(numpy.sum(steps_div_data, axis=0))

    plt.plot(full_div_x, full_div_hor, label='Polychromatic', lw=line_width)
    plt.plot(steps_div_x, steps_div_hor, label='Energy steps', lw=line_width)

    plt.xlabel("Horizontal divergence ($\mu$rad)", fontsize= f_size)    
    plt.ylabel("Normalized intensity (a.u.)", fontsize= f_size)

    ##plt.xlim(4e3, 21e4)
    ##plt.ylim(1e12, 5e15)
    plt.xticks(fontsize= f_size)
    plt.yticks(fontsize= f_size)    
    ##plt.grid(which='both', axis='y')   
    ##plt.title(f"Method: {label}, value: {plot_value}", fontsize=f_size)
    plt.legend(fontsize = f_size)
    plt.show()

    ### VERTICAL DIVERGENCE ###

    #sum all the columns and normalize to the peak
    full_div_ver = numpy.sum(full_div_data, axis=1) / max(numpy.sum(full_div_data, axis=1))
    steps_div_ver = numpy.sum(steps_div_data, axis=1) / max(numpy.sum(steps_div_data, axis=1))

    plt.plot(full_div_y, full_div_ver, label='Polychromatic', lw=line_width)
    plt.plot(steps_div_y, steps_div_ver, label='Energy steps', lw=line_width)

    plt.xlabel("Vertical divergence ($\mu$rad)", fontsize= f_size)    
    plt.ylabel("Normalized intensity (a.u.)", fontsize= f_size)

    ##plt.xlim(4e3, 21e4)
    ##plt.ylim(1e12, 5e15)
    plt.xticks(fontsize= f_size)
    plt.yticks(fontsize= f_size)    
    ##plt.grid(which='both', axis='y')   
    ##plt.title(f"Method: {label}, value: {plot_value}", fontsize=f_size)
    plt.legend(fontsize = f_size - 6)
    plt.show()    
    
if __name__=='__main__': 
    #pass
    
    code_undul_phot = 'internal'
    
    #run_source_full_harmonic('ESRF_ID06_EBS_CPMU18_1.json', res_photon_energy=10e3, energy_min= 9600,
    #          energy_max=10200, energy_points=101, harmonic=1, code_undul_phot = code_undul_phot,
    #          emittance = True, nrays=11000000, seed=5676561, flag_energy_spread=0,
    #          energy_spread=0.001, nbins=201, save_file=True)
    
    #run_source_by_energy_step('ESRF_ID06_EBS_CPMU18_1.json', res_photon_energy=10e3, energy_min= 9600,
    #          energy_max=10200, energy_points=101, harmonic=1, code_undul_phot = code_undul_phot,
    #          emittance = True, nrays=100000, seed=0, flag_energy_spread=0,
    #          energy_spread=0.001, nbins=201, save_file = True)

    compare_results('s4_full_harmonic_profiles.h5', 's4_e_step_full_harmonic_profiles.h5', f_size=24)