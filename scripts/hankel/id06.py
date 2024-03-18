
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
    h_slit_points=5,
    v_slit_points=10001,
    number_of_trajectory_points=1666,
    traj_method=1, # 0=TRAJECTORY_METHOD_ANALYTIC, 1=TRAJECTORY_METHOD_ODE
    rad_method=2, # 0=RADIATION_METHOD_NEAR_FIELD, 1= RADIATION_METHOD_APPROX, 2=RADIATION_METHOD_APPROX_FARFIELD
    )

output_wavefront = light_source.get_wavefront().get_Wavefront1D_from_profile(1, 0.0)


# if plot_from_oe <= 0: plot(output_wavefront.get_abscissas(),output_wavefront.get_intensity(),title='SOURCE')


# filename = "id06_farfield_y.dat"
# f = open(filename,'w')
y = output_wavefront.get_abscissas()
# z = output_wavefront.get_intensity()
c = output_wavefront.get_complex_amplitude()
# for i in range(y.size):
#     f.write("%g  %g\n" % (y[i], z[i]))
# f.close()
# print("File written to disk: %s" % filename)


##########  OPTICAL SYSTEM ##########

# ##########  OPTICAL ELEMENT NUMBER 1 ##########
from hankel_propagator import hankel_propagate

igood = numpy.argwhere(y >= 0)

yy = y[igood]
cc = c[igood]
yy.shape = -1
cc.shape = -1

lambda1 = output_wavefront.get_wavelength()
z_max = 100.0
# if plot_from_oe <= 0: plot(yy, numpy.abs(cc)**2, title='RADIAL SOURCE')

print(yy.shape, cc.shape, lambda1)
Erz = hankel_propagate(cc, yy, lambda1, z_max=-z_max)
Isource = numpy.abs(cc)**2
Ibackpropagated = numpy.abs(Erz)**2
if plot_from_oe <= 0: plot(yy, Isource / Isource.max(),
                           yy, Ibackpropagated / Ibackpropagated.max(),
                           title='RADIAL BACKPROPAGATED', legend=['source','backpropagated'])


filename = "id06_backpropagated_y.dat"
f = open(filename,'w')
yy0 = numpy.flip(yy)
Ibackpropagated0 = numpy.flip(Ibackpropagated)
for i in range(yy0.size):
    if yy0[i] < 30e-6: f.write("%g  %g\n" % (-yy0[i], Ibackpropagated0[i]))
for i in range(1, yy.size):
    if yy[i] < 30e-6: f.write("%g  %g\n" % (yy[i], Ibackpropagated[i]))
f.close()
print("File written to disk: %s" % filename)