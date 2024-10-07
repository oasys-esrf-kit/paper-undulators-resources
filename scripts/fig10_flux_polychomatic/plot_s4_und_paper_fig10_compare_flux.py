import pandas as pd
import matplotlib.pyplot as plt
import numpy

""" This simple script plots the comparison of polychromatic source S4:undulator
    theoretical flux and the actual energies of the generated rays, it reads two
    files created directly by S4 GUI: Undulator Light Source
    
    You can check the OASYS workflow (in this same folder): s4_und_flux_polychromatic_source.ows
    
    - The theoretical flux file was created by saving the espectrum in a CSV file from the tab:
       Undulators Plots/Flux spectrum/save_curve
       
    - The histogram ray spectrum was saved by using the tab: 
       Plots/Energy/save_curve    
    
    """


def plot_flux_spectrum(s4_get_flux, s4_histo_flux, font_size= 12):
    
    """ This function reads the files saved directly from the S4 GUI: Undulator Light Source,
    then, it weights the flux from the histogram to convert it to Ph/s/eV
    """         

    df_s4_get_f = pd.read_csv(s4_get_flux, sep=',', comment = '#', engine='python')
    e_get = df_s4_get_f['Photon energy [eV]']
    f_get = df_s4_get_f['Photons/s/0.1%bw']
    ff_get  = f_get  / (e_get * 1e-3) # in ph/eV

    df_s4_histo_f = pd.read_csv(s4_histo_flux, sep=',', comment = '#', engine='python', skiprows=1)
    e_histo = df_s4_histo_f.iloc[: , 0]
    f_histo = df_s4_histo_f.iloc[: , 1]    

    # Manolo fashion calibration histogram
    n_total = numpy.trapz(ff_get, e_get)
    h_total = numpy.trapz(f_histo, e_histo)
    
    plt.plot(e_get, ff_get, label='Analytical')
    plt.plot(e_histo, f_histo * n_total / h_total, label='Histogram')    
       
    plt.xlabel("Photon energy [eV]", fontsize=font_size)
    plt.ylabel("Flux [Photons/s/eV]", fontsize=font_size)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size) 
    plt.legend(fontsize=font_size)
    #plt.title()
    plt.show()   
    
if __name__=='__main__': 
    
    #pass    

    plot_flux_spectrum('s4_spectrum_1st_harm_slit_0.001_poly.csv', 's4_spectrum_1st_harm_slit_0.001_poly_histo.csv', font_size=12)
    plot_flux_spectrum('s4_spectrum_3rd_harm_slit_0.0005_poly.csv', 's4_spectrum_3rd_harm_slit_0.0005_poly_histo.csv', font_size=12)
