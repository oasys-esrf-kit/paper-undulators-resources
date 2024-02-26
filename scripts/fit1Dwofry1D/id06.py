
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
plot_from_oe = 0 # set to a large number to avoid plots


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
    photon_energy=10000,
    h_slit_points=20,
    v_slit_points=250,
    number_of_trajectory_points=1666,
    traj_method=1, # 0=TRAJECTORY_METHOD_ANALYTIC, 1=TRAJECTORY_METHOD_ODE
    rad_method=2, # 0=RADIATION_METHOD_NEAR_FIELD, 1= RADIATION_METHOD_APPROX, 2=RADIATION_METHOD_APPROX_FARFIELD
    )

output_wavefront = light_source.get_wavefront().get_Wavefront1D_from_profile(1, 0.0)


if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='SOURCE')


filename = "id06_farfield_y.dat"
f = open(filename,'w')
y = output_wavefront.get_abscissas()
z = output_wavefront.get_intensity()
for i in range(y.size):
    f.write("%g  %g\n" % (y[i], z[i]))
f.close()
print("File written to disk: %s" % filename)


##########  OPTICAL SYSTEM ##########





##########  OPTICAL ELEMENT NUMBER 1 ##########



input_wavefront = output_wavefront.duplicate()
from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

optical_element = WOScreen1D()

# no drift in this element
output_wavefront = optical_element.applyOpticalElement(input_wavefront)


#
#---- plots -----
#
if plot_from_oe <= 1: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 1')


##########  OPTICAL ELEMENT NUMBER 2 ##########



input_wavefront = output_wavefront.duplicate()
from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen1D

optical_element = WOScreen1D()

# drift_before -100 m
#
# propagating
#
#
propagation_elements = PropagationElements()
beamline_element = BeamlineElement(optical_element=optical_element,    coordinates=ElementCoordinates(p=-100.000000,    q=0.000000,    angle_radial=numpy.radians(0.000000),    angle_azimuthal=numpy.radians(0.000000)))
propagation_elements.add_beamline_element(beamline_element)
propagation_parameters = PropagationParameters(wavefront=input_wavefront,    propagation_elements = propagation_elements)
#self.set_additional_parameters(propagation_parameters)
#
propagation_parameters.set_additional_parameters('magnification_x', 0.0125)
#
propagator = PropagationManager.Instance()
try:
    propagator.add_propagator(FresnelZoom1D())
except:
    pass
output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_1D')


#
#---- plots -----
#
if plot_from_oe <= 2: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='OPTICAL ELEMENT NR 2')

filename = "id06_backpropagated_y.dat"
f = open(filename,'w')
y = output_wavefront.get_abscissas()
z = output_wavefront.get_intensity()
for i in range(y.size):
    f.write("%g  %g\n" % (y[i], z[i]))
f.close()
print("File written to disk: %s" % filename)