################## PYTHON CODE #########################


#
# Import section
#
import numpy

from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters

from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D

from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
from wofryimpl.propagator.propagators2D.fresnel import Fresnel2D
from wofryimpl.propagator.propagators2D.fresnel_convolution import FresnelConvolution2D
from wofryimpl.propagator.propagators2D.fraunhofer import Fraunhofer2D
from wofryimpl.propagator.propagators2D.integral import Integral2D
from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D

from srxraylib.plot.gol import plot, plot_image
plot_from_oe = 0 # set to a large number to avoid plots


##########  SOURCE ##########


#
# create output_wavefront
#
#
from orangecontrib.esrf.wofry.util.light_source import WOPySRULightSource # TODO: from wofryimpl...
light_source = WOPySRULightSource.initialize_from_keywords(
    energy_in_GeV=6,
    current=0.2,
    K_vertical=1.34109,
    period_length=0.018,
    number_of_periods=111.111,
    distance=100,
    gapH=0.005, # 0.00325963,
    gapV=0.005, # 0.00325963,
    photon_energy=10000,
    h_slit_points=250,
    v_slit_points=250,
    number_of_trajectory_points=2222,
    traj_method=1,
    rad_method=2,)

output_wavefront = light_source.get_wavefront()


if plot_from_oe <= 0: plot_image(output_wavefront.get_intensity(),output_wavefront.get_coordinate_x(),output_wavefront.get_coordinate_y(),aspect='auto',title='SOURCE')


z = output_wavefront.get_intensity()
x = output_wavefront.get_coordinate_x()
y = output_wavefront.get_coordinate_y()


filename = "id06_farfield_x.dat"
f = open(filename,'w')
for i in range(x.size):
    f.write("%f  %f\n" % (x[i], z[i, y.size // 2]))
f.close()
print("File written to disk: %s" % filename)

filename = "id06_farfield_y.dat"
f = open(filename,'w')
for i in range(y.size):
    f.write("%f  %f\n" % (y[i], z[x.size // 2, i]))
f.close()
print("File written to disk: %s" % filename)


##########  OPTICAL SYSTEM ##########





##########  OPTICAL ELEMENT NUMBER 1 ##########



input_wavefront = output_wavefront.duplicate()
from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen

optical_element = WOScreen()

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
propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
propagation_parameters.set_additional_parameters('magnification_x', 0.0125)
propagation_parameters.set_additional_parameters('magnification_y', 0.0125)
#
propagator = PropagationManager.Instance()
try:
    propagator.add_propagator(FresnelZoomXY2D())
except:
    pass
output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,    handler_name='FRESNEL_ZOOM_XY_2D')


#
#---- plots -----
#
if plot_from_oe <= 1: plot_image(output_wavefront.get_intensity(),output_wavefront.get_coordinate_x(),output_wavefront.get_coordinate_y(),aspect='auto',title='OPTICAL ELEMENT NR 1')

z = output_wavefront.get_intensity()
x = output_wavefront.get_coordinate_x()
y = output_wavefront.get_coordinate_y()

# filename = "id06_backpropagated.dat"
# f = open(filename,'w')
# for i in range(y.size):
#     f.write("%f  %f\n" % (y[i], z[x.size // 2, i]))
# f.close()
# print("File written to disk: %s" % filename)

filename = "id06_backpropagated_x.dat"
f = open(filename,'w')
for i in range(x.size):
    f.write("%g  %g\n" % (x[i], z[i, y.size // 2]))
f.close()
print("File written to disk: %s" % filename)

filename = "id06_backpropagated_y.dat"
f = open(filename,'w')
for i in range(y.size):
    f.write("%g  %g\n" % (y[i], z[x.size // 2, i]))
f.close()
print("File written to disk: %s" % filename)
################## END PYTHON CODE #########################