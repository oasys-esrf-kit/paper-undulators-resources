
#
# Import section
#
import numpy

from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters

from wofry.propagator.wavefront1D.generic_wavefront import GenericWavefront1D

from wofryimpl.propagator.propagators1D.fresnel import Fresnel1D
from wofryimpl.propagator.propagators1D.fresnel_convolution import FresnelConvolution1D
from wofryimpl.propagator.propagators1D.fraunhofer import Fraunhofer1D
from wofryimpl.propagator.propagators1D.integral import Integral1D
from wofryimpl.propagator.propagators1D.fresnel_zoom import FresnelZoom1D
from wofryimpl.propagator.propagators1D.fresnel_zoom_scaling_theorem import FresnelZoomScaling1D

from srxraylib.plot.gol import plot, plot_image
plot_from_oe = 1000 # set to a large number to avoid plots


def run_beamline(shift=0.0):

    photon_energy = 10000.0 + shift
    ##########  SOURCE ##########


    #
    # create output_wavefront
    #
    #
    from wofryimpl.propagator.light_source_pysru import WOPySRULightSource
    light_source = WOPySRULightSource.initialize_from_keywords(
        energy_in_GeV=6,
        current=0.2,
        K_vertical=1.34109,
        period_length=0.018,
        number_of_periods=111.111,
        distance=100,
        gapH=0.005,
        gapV=0.005,
        photon_energy=photon_energy,
        h_slit_points=3,
        v_slit_points=250,
        number_of_trajectory_points=1666,
        traj_method=0, # 0=TRAJECTORY_METHOD_ANALYTIC, 1=TRAJECTORY_METHOD_ODE
        rad_method=2, # 0=RADIATION_METHOD_NEAR_FIELD, 1= RADIATION_METHOD_APPROX, 2=RADIATION_METHOD_APPROX_FARFIELD
        )

    output_wavefront = light_source.get_wavefront().get_Wavefront1D_from_profile(1, 0.0)


    if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='SOURCE')

    if shift < 0:
        subgroupname = "wfr_%04d" % shift
    else:
        subgroupname = "wfr_%03d" % shift
    output_wavefront.save_h5_file("wfr1D_farfield.h5", subgroupname=subgroupname,
                                  intensity=True, phase=False, overwrite=False, verbose=True)




if __name__ == "__main__":
    # delete  wfr1D_farfield.h5
    Shift = numpy.arange(-1000, 1000, 2)
    print(Shift)
    for shift in Shift:
        print(">>>> shift: ", shift)
        run_beamline(shift=shift)
