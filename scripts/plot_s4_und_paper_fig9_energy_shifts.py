import numpy
from syned.util.json_tools import load_from_json_file
import matplotlib.pyplot as plt

""" In this script we have a set of functions to calculate the data and then use
    it to plot the results. In our S4 undulator paper (ver 04-Oct-24) correspond to figures:

    - Figure 9 (Red and Blue shifts) 
    
    This script does not save any data and plot directly the results.
    """

def run_source(syned_json_file, photon_energy=10e3, harmonic=1,
                   code_undul_phot = 'internal', emittance = True,
                   nrays=100000, seed=5676561, flag_energy_spread=0,
                   energy_spread=0.001, font_size=12):
    
    """ Using SHADOW4 Undulator, this function calculates the vertical
    divergence source for 3 energies and a given harmonic (using Elleaume style)
    photon_energy = On-resonance:

        - 'E$_{n}$(1-1/$Nn$)'       (Red shift)
        - 'E$_{n}$'                 (On-Resonance)
        - 'E$_{n}$(1+0.4/$Nn$)'     (Blue shift)

    for code_undul_phot = 'internal', 'pysru' and 'srw'   
    and plot them together.             
         """

    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam
    
    syned_obj = load_from_json_file(syned_json_file)
    e_beam = syned_obj.get_electron_beam()
    undulator = syned_obj.get_magnetic_structure()    

    k_value = undulator.get_K_from_photon_energy(photon_energy, e_beam.gamma(),
                                                harmonic=harmonic)
    
    maxangle = 3 * 0.69 * (1/e_beam.gamma())*numpy.sqrt((1.0/(2.0*harmonic * undulator.number_of_periods())) * (1.0 + k_value**2/2.0))
    print(f'For an energy of {photon_energy} and harmonic {harmonic} we have a K of {k_value}')

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

    if emittance:
        flag_emittance = 1             
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = e_beam.get_sigmas_all()
    else:
        flag_emittance = 0
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = (0.0, 0.0, 0.0, 0.0) 
    
    electron_beam.set_sigmas_all(sigma_x=sigma_x,sigma_y=sigma_y,sigma_xp=sig_div_x,sigma_yp=sig_div_y)
    
    energies = [photon_energy * (1 - 1/(undulator.number_of_periods() * harmonic)) ,\
               photon_energy, photon_energy * (1 + 0.4/(undulator.number_of_periods() * harmonic))]
    
    e_labels = ['E$_{n}$(1-1/$Nn$)', 'E$_{n}$', 'E$_{n}$(1+0.4/$Nn$)']
    l_styles = [':', '-', '--']
    colors = ['r', 'g', 'b']
    

    for i, energy_step in enumerate(energies):
    # magnetic structure
        from shadow4.sources.undulator.s4_undulator import S4Undulator
        source = S4Undulator(
            K_vertical        = k_value, # syned Undulator parameter
            period_length     = undulator.period_length(), # syned Undulator parameter
            number_of_periods = undulator.number_of_periods(), # syned Undulator parameter
            emin              = energy_step, # Photon energy scan from energy (in eV)
            emax              = energy_step, # Photon energy scan to energy (in eV)
            ng_e              = 1, # Photon energy scan number of points
            maxangle          = maxangle, # Maximum radiation semiaperture in RADIANS
            ng_t              = 100, # Number of points in angle theta
            ng_p              = 11, # Number of points in angle phi
            ng_j              = 20, # Number of points in electron trajectory (per period) for internal calculation only
            code_undul_phot   = code_undul_phot, # internal, pysru, srw
            flag_emittance    = flag_emittance, # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_size         = 2, # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
            distance          = 100.0, # distance to far field plane
            srw_range         = 0.05, # for SRW backpropagation, the range factor
            srw_resolution    = 50, # for SRW backpropagation, the resolution factor
            srw_semianalytical= 0, # for SRW backpropagation, use semianalytical treatement of phase            
            magnification     = magnification, # for internal/wofry backpropagation, the magnification factor
            flag_backprop_recalculate_source = flag_backprop_recalculate_source, # for internal or pysru/wofry backpropagation: source reused (0) or recalculated (1)
            flag_backprop_weight = flag_backprop_weight, # for internal or pysru/wofry backpropagation: apply Gaussian weight to amplitudes
            weight_ratio         = weight_ratio, # for flag_backprop_recalculate_source=1, the ratio value sigma/window halfwidth
            flag_energy_spread   = flag_energy_spread, # for monochromatod sources, apply (1) or not (0) electron energy spread correction
            )
        
        
        # light source
        from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource
        light_source = S4UndulatorLightSource(name='undulator', 
                                            electron_beam=electron_beam,
                                            magnetic_structure=source,
                                            nrays=nrays,
                                            seed=seed)
        beam = light_source.get_beam()
        
        # test plot
        tkt = beam.histo1(6, nbins=101, nolost=0, ref=23)
        z = tkt['bin_center'] * 1e6  #microns
        histogram = tkt['histogram'] / max(tkt['histogram']) #normalized to the peak

        plt.plot(z, histogram, linestyle=l_styles[i], label=e_labels[i], color=colors[i])

    plt.xlabel("Vertical divergence [$\mu$rad]", fontsize=font_size)
    plt.ylabel("Normalized intensity [a. u.]", fontsize=font_size)    
    plt.legend(fontsize=font_size)
    plt.show()

if __name__=='__main__': 
    #pass  

    #run_pysru_wofry('ESRF_ID06_EBS_CPMU18_1.json', 10000, 1,
    #               code_undul_phot = 'srw', emittance = True,
    #               nrays=500000, seed=5676561, flag_energy_spread=0)

    code_undul_phot = 'srw'

    run_source('ESRF_ID06_EBS_CPMU18_1.json', 50000, 5,
               code_undul_phot = code_undul_phot, emittance = True,
               nrays=500000, seed=5676561, flag_energy_spread=0, font_size=12)
    
   # run_pysru_wofry('ESRF_ID06_EBS_CPMU18_1.json', 50000, 5,
   #            code_undul_phot = code_undul_phot, emittance = True,
   #            nrays=500000, seed=5676561, flag_energy_spread=1, energy_spread=0.001, font_size=12)