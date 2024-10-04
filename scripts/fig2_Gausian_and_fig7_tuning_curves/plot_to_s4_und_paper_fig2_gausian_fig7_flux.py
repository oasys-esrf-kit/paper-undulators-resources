import numpy
import matplotlib.pyplot as plt
from syned.util.json_tools import load_from_json_file
from srxraylib.plot.gol import plot as manolo_plot
import pandas as pd
from matplotlib import colors as mcolors
from shadow4.sources.s4_electron_beam import S4ElectronBeam
from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
from oasys.util.oasys_util import get_sigma

from collections import OrderedDict
from xoppylib.sources.srundplug import tuning_curves_on_slit


""" In this script we have a set of function to calculate the data and then use
    it to plot the results. In our S4 undulator paper (ver 04-Oct-24) correspond to figures:

    - Figure 2 (Gaussian sizes and divergences)

    - Figure 7 (Flux comparison between Gaussian source flux and XOPPY)

""" 

def get_gauss_lightsource(syned_json_file, photon_energy, energy_spread,
                          flag_energy_spread=1, harmonic_number=1, emittance=True,
                          nrays=50000, seed=5676561):
    
    """ This function gets the SHADOW4 Gaussian light_source from a SYNED JSON
    File, also the maximum K value (which is in the JSON file)"""

    syned_obj = load_from_json_file(syned_json_file)
    e_beam = syned_obj.get_electron_beam()
    undulator = syned_obj.get_magnetic_structure()  

    k_max = undulator.K_vertical()  

    electron_beam = S4ElectronBeam(energy_in_GeV=e_beam.energy(),
                                  energy_spread=energy_spread, current=e_beam.current())

    if emittance:
        flag_emittance = 1             
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = e_beam.get_sigmas_all()
    else:
        flag_emittance = 0
        (sigma_x, sig_div_x, sigma_y, sig_div_y) = (0.0, 0.0, 0.0, 0.0) 
    
    electron_beam.set_sigmas_all(sigma_x=sigma_x,sigma_y=sigma_y,sigma_xp=sig_div_x,sigma_yp=sig_div_y)
    #electron_beam.set_sigmas_all(sigma_x=3.01836e-05,sigma_y=4.36821e-06,sigma_xp=3.63641e-06,sigma_yp=1.37498e-06)
    
    # magnetic structure
    
    source = S4UndulatorGaussian(
        period_length     = undulator.period_length(),     # syned Undulator parameter (length in m)
        number_of_periods = undulator.number_of_periods(), # syned Undulator parameter
        photon_energy     = photon_energy, # Photon energy (in eV)
        delta_e           = 0.0, # Photon energy width (in eV)
        ng_e              = 50, # Photon energy scan number of points
        flag_emittance    = flag_emittance, # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_energy_spread = flag_energy_spread, # when sampling rays: Use e- energy spread (0=No, 1=Yes)
        harmonic_number    = harmonic_number, # harmonic number
        flag_autoset_flux_central_cone  = 0, # value to set the flux peak
        flux_central_cone  = 10000000000.0, # value to set the flux peak
        )    
    
    # light source
    
    light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam,
                                                  magnetic_structure=source,nrays=nrays,seed=seed)

    return light_source, k_max

def get_analytic_gauss_source(syned_json_file, photon_energy, energy_spread,
                    dimension='hor_size', k_min=0.2, max_harmonic_number=11,
                    flag_energy_spread=1, emittance=True, save_file=True):
    
    """ This function gets the photon energies (from k_max to k_min) and the standard deviation of a
        given Gauss source dimension: hor_size, ver_size, hor_div, ver_div, for 
        odd harmonics up to the max_harmonic_number, 
        it can save a CSV file with the photon energy in one column and
        the dimention in another """    
    
    var_dic = {'hor_size': 1, 'ver_size': 2, 'hor_div': 3, 'ver_div': 4}
     
    light_source, k_max = get_gauss_lightsource(syned_json_file, photon_energy,
                                               energy_spread, flag_energy_spread=flag_energy_spread,
                                               emittance=emittance)
    
    gauss_variables = light_source.get_size_and_divergence_vs_photon_energy(Kmin=k_min,
                                                                            Kmax=k_max,
                                                                            max_harmonic_number=max_harmonic_number)    
    
    photon_energies = gauss_variables[0]
    source_dim = gauss_variables[var_dic[dimension]]

    if save_file:
        for i in range(len(photon_energies)):
            if i == 0:
                df_analytic_gauss = pd.DataFrame({f'Photon Energy {i}':photon_energies[i], f'{dimension} {i}':source_dim[i]})
            else:
                df_analytic_gauss[f'Photon Energy {1+((2*i)-2)}'] = photon_energies[i]
                df_analytic_gauss[f'{dimension} {1+((2*i)-2)}'] = source_dim[i]
        df_analytic_gauss.to_csv(f'Analytic_Gauss_{dimension}_spread_{energy_spread}.csv', index=False)
        print(f'file Analytic_Gauss_{dimension}_spread_{energy_spread}.csv has been saved to disk')
    
    return photon_energies, source_dim


def run_s4_und_gauss(syned_json_file, photon_energy, energy_spread, dimension='hor_size',
                    num_points=3, num_rep=50, k_min=0.2, max_harmonic_number=11,
                    flag_energy_spread=1, emittance=True, nrays=10000, seed=5676561,
                    save_file=True):
        
    """ By generating a SHADOW4 source, this function calculates the standard 
        deviation of a given dimension, for a set of energies (num points) which
        are equally distributed for each harmonic, also by repeating a num_rep
        times the source generation can calculates a statistical STD for each value,
        these results are all odd harmonics up to max_harmonic_number.
        It can save a CSV file with the 
        photon energy in one column and the dimention in another """ 

    col_dic = {'hor_size':1, 'ver_size':3, 'hor_div':4, 'ver_div':6}

    light_source, k_max = get_gauss_lightsource(syned_json_file, photon_energy,
                                               energy_spread, flag_energy_spread=flag_energy_spread,
                                               emittance=emittance, nrays=nrays, seed=seed)

    beam = light_source.get_beam()

    #first it gets the photon_energies
    photon_energies,_,_,_,_,_ = light_source.get_size_and_divergence_vs_photon_energy(Kmin=k_min,
                                                                                     Kmax=k_max,
                                                                                     max_harmonic_number=max_harmonic_number)
    
    s4_energy = []
    s4_dim = []
    s4_dim_std = []

    for i in range(1, len(photon_energies)):
        harmonic = 1 + (2*(i-1))
        interval = int(len(photon_energies[i])/(num_points+1))

        tmp_energy = []
        tmp_sigma = []
        tmp_dim_stad = []
        
        for j in range(1, num_points+1):
            
            step_energy = photon_energies[i][j*interval]
            
            tmp_energy.append(step_energy)
            stat_dim = []
            #for repeat the calculation for the num_rep:
            for k in range(num_rep):
                light_source, k_max = get_gauss_lightsource(syned_json_file,
                                                           step_energy , energy_spread,
                                                           flag_energy_spread=flag_energy_spread,
                                                           harmonic_number=harmonic,
                                                           emittance=emittance,
                                                           nrays=nrays,
                                                           seed=seed)
                beam = light_source.get_beam()
                #fwhm = beam.histo1(col_dic[dimension], nbins=301, ref=23)['fwhm']
                #sigma_k = get_sigma(beam.histo1(col_dic[dimension], nbins=301, ref=23)['histogram'],
                #                   beam.histo1(col_dic[dimension], nbins=301, ref=23)['bin_center'])
                
                sigma_k = beam.get_standard_deviation(col_dic[dimension], nolost=1, ref=1)
                print(f"For {step_energy} and harmonic {harmonic} we have a sigma of {sigma_k} m")
                stat_dim.append(sigma_k)
            tmp_sigma.append(numpy.average(stat_dim))
            tmp_dim_stad.append(numpy.std(stat_dim))
        
        s4_energy.append(tmp_energy)
        s4_dim.append(tmp_sigma)
        s4_dim_std.append(tmp_dim_stad)

    if save_file:
        for i in range(len(s4_energy)):
            if i == 0:
                df_s4_gauss = pd.DataFrame({f'Photon Energy {(2*i)+1}':s4_energy[i],
                                           f's4_{dimension} {(2*i)+1}':s4_dim[i],
                                           f's4_std_{dimension} {(2*i)+1}':s4_dim_std[i]})
            else:
                df_s4_gauss[f'Photon Energy {(2*i)+1}'] = s4_energy[i]
                df_s4_gauss[f's4_{dimension} {(2*i)+1}'] = s4_dim[i]
                df_s4_gauss[f's4_std_{dimension} {(2*i)+1}'] = s4_dim_std[i]
        df_s4_gauss.to_csv(f'shadow4_Gauss_{dimension}_spread_{energy_spread}.csv', index=False)
        print(f'file shadow4_Gauss_{dimension}_spread_{energy_spread}.csv has been saved to disk')

    return s4_energy, s4_dim, s4_dim_std

def plot_gauss_shadow4(analytic_gauss_csv, shadow4_gauss_csv, dimension = 'hor_size',
                      var_name= 'Horizontal size', units='$\mu$m', error_bar=True,
                      save_fig=True, f_size=11):

    """ This function plots for a given dimension the analityc Gauss results:
     - Results from function: get_analytic_gauss_source
     - Rsult from ray generated source run_s4_und_gauss (fun run_s4_und_gauss)
           
     dimension options: hor_size, ver_size, hor_div, ver_div    
          """

    df_analytic = pd.read_csv(analytic_gauss_csv, sep=',', comment = '#', engine='python')

    df_s4_gauss = pd.read_csv(shadow4_gauss_csv, sep=',', comment = '#', engine='python')

    colors = list((dict(mcolors.TABLEAU_COLORS)).values())

    for i in range(int(len(df_analytic.columns)/2)):
        if i==0:
            plt.plot(df_analytic.iloc[:,i]*1e-3, df_analytic.iloc[:,i+1] * 1e6, label='Zero energy spread')
        else:
            plt.plot(df_analytic.iloc[:, 2*i]*1e-3, df_analytic.iloc[:,(2*i)+1] * 1e6, label=f'Harmonic {(2*i)-1}')
            if error_bar:
                plt.errorbar(df_s4_gauss[f'Photon Energy {(2*i)-1}']*1e-3, df_s4_gauss[f's4_{dimension} {(2*i)-1}']* 1e6,
                            yerr=df_s4_gauss[f's4_std_{dimension} {(2*i)-1}'] * 1e6, marker='o', color=colors[i],
                            linestyle='none', capsize=5)
            else:
                plt.plot(df_s4_gauss[f'Photon Energy {(2*i)-1}'] * 1e-3,
                        df_s4_gauss[f's4_{dimension} {(2*i)-1}']* 1e6,
                        marker='o', color=colors[i], linestyle='none')

    plt.xlabel("Photon energy [keV]", fontsize= f_size)    
    plt.ylabel(f"{var_name} ({units})", fontsize= f_size)

    #plt.xlim(4e3, 21e4)
    #plt.ylim(1e12, 5e15)
    plt.xticks(fontsize= f_size)
    plt.yticks(fontsize= f_size)    
    #plt.grid(which='both', axis='y')   
    #plt.title("EBS ID06 CMPU18, On Resonance, $\epsilon$=0, $\delta$=0.001", fontsize=f_size)
    plt.legend(fontsize = f_size - 6)

    if save_fig:
        plt.savefig(f"Source_Gauss_{dimension}.png", dpi=600)
        print(f"Figure Source_Gauss_{dimension}.png has been saved to disk")
          
    plt.show()    

def get_flux_xoppy(syned_json_file, Kpoints=3, harms="1, 3, 5, 7, 9",
                  energy_spread=0.001, kmin=0.002, flag_energy_spread=1,
                  emittance=True, nrays=10000, seed=5676561, save_file=True):
    
    syned_obj = load_from_json_file(syned_json_file)
    e_beam = syned_obj.get_electron_beam()
    undulator = syned_obj.get_magnetic_structure()

    kmax = undulator.K_vertical() 

    (sigma_x, sig_div_x, sigma_y, sig_div_y) = e_beam.get_sigmas_all() 


    bl = OrderedDict()
    bl['ElectronBeamDivergenceH'] = sig_div_x
    bl['ElectronBeamDivergenceV'] = sig_div_y
    bl['ElectronBeamSizeH']       = sigma_x
    bl['ElectronBeamSizeV']       = sigma_y
    bl['ElectronCurrent']         = e_beam.current()
    bl['ElectronEnergy']          = e_beam.energy()
    bl['ElectronEnergySpread']    = energy_spread
    bl['NPeriods']                = undulator.number_of_periods()
    bl['PeriodID']                = undulator.period_length()
    bl['distance']                = 10.0
    bl['gapH']                    = 0.008
    bl['gapV']                    = 0.008
    bl['gapHcenter']              = 0.0
    bl['gapVcenter']              = 0.0
    
    harmonics = harms.split(",")
    
    K_scan,harmonics,power_array, energy_values_at_flux_peak,flux_values = tuning_curves_on_slit(bl,
        Kmin=kmin,
        Kmax=kmax,
        Kpoints=Kpoints,
        harmonics=harmonics,
        zero_emittance=False,
        do_plot_peaks=False,
        code="srw")
    
    print(energy_values_at_flux_peak)
    
    for i in range(len(energy_values_at_flux_peak[0])):
        harm = (2*i) + 1
        photon_energy_harm = energy_values_at_flux_peak[:, i]
        print(photon_energy_harm)
        s4_flux = []

        for photon_energy in photon_energy_harm:
            light_source, k_max = get_gauss_lightsource(syned_json_file,
                                       photon_energy , energy_spread,
                                       flag_energy_spread=flag_energy_spread,
                                       harmonic_number=harm,
                                       emittance=emittance,
                                       nrays=nrays,
                                       seed=seed)
            
            s4_flux.append(light_source.get_flux_central_cone())

        

        if i == 0:
            df_flux = pd.DataFrame({f'Photon Energy {harm}':photon_energy_harm,
                                    f'xoppy_flux_{harm}':flux_values[:, i],
                                    f's4_flux_{harm}':s4_flux})
            print(df_flux)
                
        else:            
            df_flux[f'Photon Energy {harm}'] = photon_energy_harm
            df_flux[f'xoppy_flux_{harm}'] = flux_values[:, i]
            df_flux[f's4_flux_{harm}'] = s4_flux

    if save_file:
        df_flux.to_csv(f'Flux_comparison_xoppy_s4.csv', index=False)
        print(f'Flux_comparison_xoppy_s4.csv.csv has been saved to disk')
      

    return df_flux

def plot_flux(flux_csv, save_fig=True):

    f_size = 11

    df_flux = pd.read_csv(flux_csv, sep=',', comment = '#', engine='python')

    for i in range(int(df_flux.shape[1]/3)):
        harmonic = (2*i)+1
        plt.plot(df_flux[f'Photon Energy {harmonic}']*1e-3, df_flux[f'xoppy_flux_{harmonic}'], marker='o', label=f'XOPPY_{harmonic}')
        plt.plot(df_flux[f'Photon Energy {harmonic}']*1e-3, df_flux[f's4_flux_{harmonic}'], marker='s', label=f'S4_{harmonic}')   

    plt.xlabel("Photon energy [keV]", fontsize= f_size)    
    plt.ylabel("Flux [Photons/s/0.1%BW]", fontsize= f_size)

    #plt.xlim(4e3, 21e4)
    #plt.ylim(10, 5e17)
    plt.xticks(fontsize= f_size)
    plt.yticks(fontsize= f_size)
    plt.yscale('log')    
    #plt.grid(which='both', axis='y')       
    plt.legend(fontsize = f_size)

    if save_fig:

        plt.savefig(f"Flux_comparison_xoppy_s4.png", dpi=600)
        print(f"Flux_comparison_xoppy_s4.png has been saved to disk")

    plt.show()


if __name__=='__main__':
    #pass
    #examples of use:    
   
    plotting = True
    calculating = False

    if calculating:

        dimension = 'hor_size'

        get_analytic_gauss_source('ESRF_ID06_EBS_CPMU18_1.json', 10000, 0.001,
                    dimension=dimension, k_min=0.2, max_harmonic_number=11,
                    flag_energy_spread=1, emittance=True, save_file=True)
        s4_energy, s4_dim, s4_dim_std = run_s4_und_gauss('ESRF_ID06_EBS_CPMU18_1.json',
                                                        10000, 0.001, dimension=dimension,
                                                        num_points=3, num_rep=50, k_min=0.2,
                                                        max_harmonic_number=11, flag_energy_spread=1,
                                                        emittance=True, nrays=100000, seed=0,
                                                        save_file=True)
    else:
        pass

    #get_analytic_gauss_source('ESRF_ID06_EBS_CPMU18_1.json', 10000, 0.001,
    #                dimension='hor_div', k_min=0.2, max_harmonic_number=11,
    #                flag_energy_spread=1, emittance=True, save_file=True)
    #s4_energy, s4_dim, s4_dim_std = run_s4_und_gauss('ESRF_ID06_EBS_CPMU18_1.json',
    #                                                10000, 0.001, dimension='hor_size',
    #                                                num_points=3, num_rep=50, k_min=0.2,
    #                                                max_harmonic_number=11, flag_energy_spread=1,
    #                                                emittance=True, nrays=100000, seed=0,
    #                                                save_file=True)
    
    # Figure 2

    if plotting:
        dimension = 'ver_size'
        var_name = 'Vertical size' #Vertical divergence'
        units = '$\mu$m'#'$\mu$rad'# 

        
        plot_gauss_shadow4(f'Analytic_Gauss_{dimension}_spread_0.001.csv',
                          f'shadow4_Gauss_{dimension}_spread_0.001.csv',
                           dimension = dimension, var_name = var_name,
                           units=units, error_bar=True, save_fig=True, f_size=24) 

    # for figure 7
        
    # plot_flux('Flux_comparison_xoppy_s4.csv')        