import spectra.python.spectra as spectra
import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy
import scipy.constants as codata
from syned.util.json_tools import load_from_json_file
import h5py

def run_spectra_div(res_energy, shift, syned_json_file, harmonic=1, plot=False, zero_emittance=False, zero_spread=False):

    #for now it works only for 
    f_size = 12
    # open an input file
    f = open("config_set_undulator_ang_dist_100m.json")
    #load the json file
    prm = json.load(f)
    
    #load the syned file
    syned_obj = load_from_json_file(syned_json_file)
    #get the electron beam and undulator parameters
    electron_beam = syned_obj.get_electron_beam()
    undulator = syned_obj.get_magnetic_structure()
    undulator_name = syned_obj.get_name()
    #change the undulator magnetic period lenght
    prm["Light Source"]["&lambda;<sub>u</sub> (mm)"] = undulator.period_length() * 1e3 #mm
    prm["Light Source"]["Device Length (m)"] =  undulator.length()
    
    #get the K paramaeter value for the given photon energy resonance
    #we can use the undulator  method 
    K = undulator.get_K_from_photon_energy(res_energy, electron_beam.gamma(), harmonic=harmonic)
    
    prm["Light Source"]["K value"] = K    
    prm["Configurations"]["Target Energy (eV)"] = res_energy + shift #Off resonance
    
    # This has to be tunned depending of each undulator
    prm["Configurations"]["X Range (mm)"] = [-4, 4]
    prm["Configurations"]["Points (x)"] = 101
    prm["Configurations"]["Y Range (mm)"] = [-4, 4]
    prm["Configurations"]["Points (y)"] = 101
    
    #we read the distance for calculation
    distance = prm["Configurations"]["Distance from the Source (m)"]
    
    # change flags for zero emittance and zero energy spread 
    
    prm["Accelerator"]["Options"]["Zero Emittance"]=zero_emittance            
    prm["Accelerator"]["Options"]["Zero Energy Spread"]=zero_spread   
    
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
        plt.title("{} Photon Energy {} eV".format(\
                    und_name, photon_energy), fontsize=f_size)
        
        plt.xlabel('Horizontal divergence [mrad]', fontsize=f_size)
        plt.ylabel('Vertical divegence [mrad]', fontsize=f_size)
        plt.xticks(fontsize=f_size)
        plt.yticks(fontsize=f_size)
        
        plt.show()
        
    return x_div, y_div, zdata
    
def test_save_h5(res_energy, shifts, syned_json_file, harmonic=1, spectra_data= 'off_res_div_100m_case2', save_file=True, zero_emittance=False, zero_spread=False):

        
    syned_obj = load_from_json_file(syned_json_file)
    
    flux_densities = []
    
    for shift in shifts:
        
        x_div, y_div, zdata = run_spectra_div(res_energy, shift, syned_json_file, harmonic=harmonic, zero_emittance=zero_emittance, zero_spread=zero_spread)
        flux_densities.append(zdata)
        #tz_data = numpy.array(zdata).T
        #flux_densities.append(tz_data)
    
    #flux_densities = numpy.array(flux_densities)
    #array_energies = numpy.array(photon_energies)
    #array_energies = array_energies.astype(numpy.float)   
    
    if save_file:
    
        if len(shifts) == 1:
            #single energy
            HDF5_FILE = 'Spectra_{}_{}_{}_eV'.format(spectra_data, syned_obj.get_name(),shifts[0])
        else:
            HDF5_FILE = 'test_spectra_{}_{}'.format(spectra_data, syned_obj.get_name())
        if zero_emittance:
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
        nxdata.attrs['axes'] = ['shift', 'y_div', 'x_div']
            
        # Image data
        fs = nxdata.create_dataset("flux_dens", data=flux_densities)
        fs.attrs['long_name'] = 'Flux density'
    
        es = nxdata.create_dataset('shift', data=shifts)
        es.attrs['units'] = 'eV'
        es.attrs['long_name'] = 'Photon Energy Shift (eV)'    
    
        # X axis data
        xs = nxdata.create_dataset('x_div', data=x_div)
        xs.attrs['units'] = 'mrad'
        xs.attrs['long_name'] = 'horizontal divegence (mrad)'    
    
        # Y axis data
        ys = nxdata.create_dataset('y_div', data=y_div)
        ys.attrs['units'] = 'mrad'
        ys.attrs['long_name'] = 'vertical divergence (mrad)'    
    
        f.close()
        
        print ('File {} save to the disk'.format(HDF5_FILE)) 
    else:
        pass
    
    #TODO change the variable to div in the next function
    
def save_full_h5(ref_photon_energy_file, syned_json_file, spectra_data= 'charac_at_source_point', zero_emittance=False, zero_spread=False):

    df_e_ref = pd.read_csv(ref_photon_energy_file, sep='\t', comment = '#', engine='python')
    
    syned_obj = load_from_json_file(syned_json_file)
    
    #get the energies
    list_energy = []
    for i in range(int(df_e_ref.shape[1]/2)):
        list_energy.append(df_e_ref.iloc[:, i*2])    
    
    HDF5_FILE = 'Spectra_{}_{}_full'.format(spectra_data, syned_obj.get_name())
    
    if zero_emittance:
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
    
    
    for i in range(len(list_energy)):
        flux_densities = []
        for photon_energy in list_energy[i]:
            print('Calculating SPECTRA for photon energy {} and harmonic {}'.format(photon_energy, "%02d"%(i*2+1)))
            div_x, div_y, zdata = run_spectra_div(photon_energy, syned_json_file, harmonic=((i*2)+1), zero_emittance=zero_emittance, zero_spread=zero_spread)
            flux_densities.append(zdata)
        #dump in the h5file
        print('Done for the harmonic {} dumping in the h5 file...'.format("%02d"%(i*2+1)))
        
        nxdata = nxentry.create_group('harmonic {}'.format("%02d"%(i*2+1)))
        nxdata.attrs['NX_class'] = 'NXdata'
        nxdata.attrs['signal'] = '%s'%("flux_dens")
        nxdata.attrs['axes'] = ['energy', 'y_div', 'x_div']
            
        # Image data
        fs = nxdata.create_dataset("flux_dens", data=flux_densities)
        fs.attrs['long_name'] = 'Flux density'
    
        es = nxdata.create_dataset('energy', data=list_energy[i])
        es.attrs['units'] = 'eV'
        es.attrs['long_name'] = 'Photon Energy (eV)'    
    
        # X axis data
        xs = nxdata.create_dataset('x', data=x)
        xs.attrs['units'] = 'mm'
        xs.attrs['long_name'] = 'horizontal (mm)'    
    
        # Y axis data
        ys = nxdata.create_dataset('y', data=y)
        ys.attrs['units'] = 'mm'
        ys.attrs['long_name'] = 'vertical (mm)'
    
    f.close()
    
    print ('File {} save to the disk'.format(HDF5_FILE))   
    
if __name__=='__main__':
    #run_spectra_div(10000, "ESRF_ID06_EBS_CPMU18_1.json", harmonic=1, plot=True, zero_emittance=False, zero_spread=False)
    
    #test_save_h5([9500, 10000,], "ESRF_ID06_EBS_CPMU18_1.json", harmonic=1, save_file=True)
    #test_save_h5([ 9990,  9992,  9994,  9996,  9998, 10000, 10002, 10004, 10006,
    #   10008], "ESRF_ID06_EBS_CPMU18_1.json", harmonic=1, save_file=True, zero_emittance=True, zero_spread=True)
    
    #save_full_h5('spectra_results/SPECTRA_EBS_cpmu18_hor_size.txt', 'ESRF_ID06_EBS_CPMU18_1.json', spectra_data= 'charac_at_source_point')   
    #shifts = numpy.arange(-400, 100, 1).tolist()
    shifts = numpy.arange(-1000, 1000, 2).tolist()
    #print(len(shifts))    
    test_save_h5(10000,shifts, "ESRF_ID06_EBS_CPMU18_1.json", harmonic=1, spectra_data= 'off_res_div_case2', save_file=True, zero_emittance=True, zero_spread=True)
    
    #test
    #f = open("C://Work_JR//spectra_win//config_set_undulator_at_point_source_spatial_profile.json")
    #prm = json.load(f)
    
    #print(prm["Configurations"]["X Range (mm)"])
   # print(type(prm["Configurations"]["X Range (mm)"]))
    
