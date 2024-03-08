import spectra.python.spectra as spectra
import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy
import scipy.constants as codata
from syned.util.json_tools import load_from_json_file
import h5py

def gamma(energy_in_GeV):
    return 1e9 * energy_in_GeV / (codata.m_e *  codata.c**2 / codata.e)
    
def resonance_photon_energy(energy_in_GeV, B_target, id_period, harmonic=1):
    #id_period in meters

    k = B_target * (id_period) * codata.e / (2 * numpy.pi * codata.m_e * codata.c)   

    resonance_wavelength =  ((id_period / (2.0 * gamma(energy_in_GeV) **2)) * (1 + k**2 / 2.0))/harmonic   
    
    resonance_energy =  round(codata.c * codata.h / codata.e / resonance_wavelength, 3)
    
    print(f'Resonance_energy {resonance_energy} eV')
    
    return resonance_energy


def run_spectra_div(energy_in_GeV, syned_json_file, target_photon_energy, harmonic=1, zero_emittance=False, zero_spread=False, plot=False):

    #for now it works only for 
    f_size = 12
    # open an input file
    f = open("C://Work_JR//spectra_win//config_set_undulator_ang_dist_100m.json")
    #f = open("C://Work_JR//spectra_win//test_config_set.json")
    prm = json.load(f)
    
    #load the syned file
    syned_obj = load_from_json_file(syned_json_file)
    #get the electron beam and undulator parameters
    electron_beam = syned_obj.get_electron_beam()
    undulator = syned_obj.get_magnetic_structure()
    undulator_name = syned_obj.get_name()  
    #get the K value  
    k_value = undulator.get_K_from_photon_energy(target_photon_energy, electron_beam.gamma(), harmonic=harmonic)
        
    #Give all the values for the calcualtion
    #change the electron beam energy
    prm["Accelerator"]["Energy (GeV)"] = energy_in_GeV
    
    #undulator settings
    prm["Light Source"]["&lambda;<sub>u</sub> (mm)"] = undulator.period_length() * 1e3 #mm
    prm["Light Source"]["Device Length (m)"] =  undulator.length()    
    prm["Light Source"]["K value"] = k_value
    
    
    #Define the target energy    
    prm["Configurations"]["Target Energy (eV)"] = target_photon_energy
    prm["Configurations"]["X Range (mm)"] = [-4, 4]
    prm["Configurations"]["Points (x)"] = 101
    prm["Configurations"]["Y Range (mm)"] = [-4, 4]
    prm["Configurations"]["Points (y)"] = 101
    
    #distance on propagation
    distance = prm["Configurations"]["Distance from the Source (m)"] 
    
    prm["Accelerator"]["Options"]["Zero Emittance"] = zero_emittance
    
    prm["Accelerator"]["Options"]["Zero Energy Spread"] = zero_spread
    
    prmstr = json.dumps(prm)
    
    # call solver with the input string (JSON format)
    solver = spectra.Solver(prmstr)
    
    # check if the paramete load is OK
    isready = solver.IsReady()
    if isready == False:
        print("Parameter load failed.")
        sys.exit()
    
    # start calculation
    solver.Run()
    
    data = solver.GetData()
    xyarray = data["variables"]
    dat = data["data"]
    """ for this type of calculation SPECTRA data has a structure of nested list (:,(:,:,:,:))
    in which the flux density is the following """     
    flux_dens = dat[0][:]     
        
    x = numpy.array(xyarray[0])
    y = numpy.array(xyarray[1])
        
    nx = len(x)
    ny = len(y)
    zdata = numpy.array(flux_dens).reshape([ny, nx])
    # to mili rad
    x_div = x/distance
    y_div = y/distance
    
    if plot:
        # get results
        # get titles and units
        captions = solver.GetCaptions()
        
        titles = captions["titles"]
        
        #units = captions["units"]
               
        #plt.ion()
        #plt.figure()
        plt.pcolormesh(x_div, y_div, zdata, cmap=plt.cm.viridis, shading='auto')	
        plt.colorbar().ax.tick_params(axis='y', labelsize=f_size)
        plt.title("{} Photon Energy {} eV, harmonic {}, {} GeV".format(\
                    undulator_name, target_photon_energy, harmonic, energy_in_GeV), fontsize=f_size)
        
        plt.xlabel('Horizontal divergence [mrad]', fontsize=f_size)
        plt.ylabel('Vertical divegence [mrad]', fontsize=f_size)
        plt.xticks(fontsize=f_size)
        plt.yticks(fontsize=f_size)
        
        plt.show()
        
    return x_div, y_div, zdata
    
def get_flux_density(electron_energies, syned_json_file, target_photon_energy, harmonic=1, spectra_data= 'div_100m', save_file=True, zero_emittance=False, zero_spread=False):

        
    syned_obj = load_from_json_file(syned_json_file)
    
    flux_densities = []
    
    for electron_energy in electron_energies:
        
        x_div, y_div, zdata = run_spectra_div(electron_energy, syned_json_file, target_photon_energy=target_photon_energy, harmonic=harmonic, zero_emittance=zero_emittance, zero_spread=zero_spread)
        flux_densities.append(zdata)        
    
    if save_file:
    
        if len(electron_energies) == 1:
            #single energy
            HDF5_FILE = 'Spectra_{}_{}_{}_eV'.format(spectra_data, syned_obj.get_name(), electron_energies[0])
        else:
            HDF5_FILE = 'Spectra_{}_{}'.format(spectra_data, syned_obj.get_name())
            
        if zero_emittance :
            HDF5_FILE += '_0_emitt'
        else:
           pass
        if zero_spread:
            HDF5_FILE += '_0_spread'
        else:
            pass
        
        f = h5py.File(HDF5_FILE+'.hdf5', 'a')
        
        nxentry = f.create_group('Scan')
        nxentry.attrs['NX_class'] = 'NXentry'
        
        nxdata = nxentry.create_group('harmonic {}'.format("%02d"%harmonic))
        nxdata.attrs['NX_class'] = 'NXdata'
        nxdata.attrs['signal'] = '%s'%("flux_dens")
        nxdata.attrs['axes'] = ['electron_energy', 'y_div', 'x_div']
            
        # Image data
        fs = nxdata.create_dataset("flux_dens", data=flux_densities)
        fs.attrs['long_name'] = 'Flux density'
    
        es = nxdata.create_dataset('electron_energy', data=electron_energies)
        es.attrs['units'] = 'GeV'
        es.attrs['long_name'] = 'Electron Energy (GeV)'    
    
        # X axis data
        xs = nxdata.create_dataset('x_div', data=x_div)
        xs.attrs['units'] = 'urad'
        xs.attrs['long_name'] = 'horizontal divegence (urad)'    
    
        # Y axis data
        ys = nxdata.create_dataset('y_div', data=y_div)
        ys.attrs['units'] = 'urad'
        ys.attrs['long_name'] = 'vertical divergence (urad)'    
    
        f.close()
        
        print ('File {} save to the disk'.format(HDF5_FILE)) 
    else:
        pass  
  
    
if __name__=='__main__':

    #test = resonance_photon_energy(6.01, 0.797933, 0.018, harmonic=1)

    
    #run_spectra_div(6.1, "ESRF_ID06_EBS_CPMU18_1.json", target_photon_energy=10000, harmonic=1, zero_emittance=True, zero_spread=True, plot=True )
    
    #test_save_h5([9500, 10000,], "ESRF_ID06_EBS_CPMU18_1.json", harmonic=1, save_file=True)
    
    electron_energies = numpy.linspace(5.9, 6.1, 201).tolist()
    get_flux_density(electron_energies, "ESRF_ID06_EBS_CPMU18_1.json", target_photon_energy=30000, harmonic=3, spectra_data= 'div_harm3', save_file=True, zero_emittance=True, zero_spread=True)
    
    
    #save_full_h5('spectra_results/SPECTRA_EBS_cpmu18_hor_size.txt', 'ESRF_ID06_EBS_CPMU18_1.json', spectra_data= 'charac_at_source_point')   
    
    
    #test
    #f = open("C://Work_JR//spectra_win//config_set_undulator_at_point_source_spatial_profile.json")
    #prm = json.load(f)
    
    #print(prm["Configurations"]["X Range (mm)"])
   # print(type(prm["Configurations"]["X Range (mm)"]))
    
