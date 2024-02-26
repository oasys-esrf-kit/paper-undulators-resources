import numpy

from silx.math.fit.functions import sum_gauss, sum_agauss
from silx.math.fit import fittheories
from silx.math.fit.fitmanager import FitManager


x = numpy.linspace(-4, 4, 1000)
Gamma = numpy.pi / 2 * x**2
y = (numpy.sin(Gamma) / Gamma)**2


# Fitting
p = [y.max(), 0, 2]
fit = FitManager()
fit.setdata(x=x, y=y)
fit.loadtheories(fittheories)
fit.settheory('Gaussians')
# fit.settheory('Area Gaussians')
# fit.settheory('Slit')
fit.estimate()
fit.runfit()

print("Searched parameters = %s" % p)
print("Obtained parameters : ")
dummy_list = []
for param in fit.fit_results:
    # print(param)
    print(param['name'], ' = ', param['fitresult'])
    if param['name'] == "FWHM1":
        print("sigma1 = ", param['fitresult'] / (2 * numpy.sqrt(2 * numpy.log(2))))
    dummy_list.append(param['fitresult'])
print("chisq = ", fit.chisq)

y1 = sum_gauss(x, *dummy_list)
# y1 = sum_agauss(x, *dummy_list)

print("Area original: ", y.sum() * (x[1] - x[0]))
print("Area fit: ", y1.sum() * (x[1] - x[0]))


# plot and save file
from srxraylib.plot.gol import plot
p = plot(x, y, x, y1,
         xtitle=r'$\theta_r=\theta\sqrt{L/\lambda}$', ytitle=r'$(sin(\Gamma)/\Gamma)^2$', show=0)

import matplotlib.pylab as plt
plt.savefig('undulator_divergence.eps')
plt.show()


import scipy.constants as codata
u_nperiods = 111.111
u_period = 0.018
wavelength = codata.h * codata.c / codata.e / 10000.

cte = numpy.sqrt(u_nperiods * u_period / wavelength)

print("wavelength: ", wavelength)
filename = "id06_farfield.dat"
f = open(filename,'w')
for i in range(y.size):
    f.write("%g  %g\n" % (x[i] * 100 / cte, y[i]))
f.close()
print("File written to disk: %s" % filename)




