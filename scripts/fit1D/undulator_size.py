import numpy

from silx.math.fit.functions import sum_gauss, sum_agauss
from silx.math.fit import fittheories
from silx.math.fit.fitmanager import FitManager

from srxraylib.plot.gol import plot

a = numpy.loadtxt("undulator_size.dat")
# x = a[:, 0].copy()
# y = a[:, 1].copy()
# y /= y.max()
# p = [y.max(), 0.01, 2]


x = a[:, 0].copy() * 2 * numpy.pi * numpy.sqrt(2) + 1e-8 # tryinng to reproduce Fig 3.4 in Elleaume
y = a[:, 1].copy()
y /= y.max()
# plot(x, y)

# Fitting

fit = FitManager()
fit.setdata(x=x, y=y)
fit.loadtheories(fittheories)
fit.settheory('Gaussians')
# fit.settheory('Area Gaussians')
# fit.settheory('Slit')
fit.estimate()
fit.runfit()

# print("Searched parameters = %s" % p)
print("Obtained parameters : ")
dummy_list = []
for param in fit.fit_results:
    # print(param)
    print(param['name'], ' = ', param['fitresult'])
    if param['name'] == "FWHM1":
        print("sigma1 = ", param['fitresult'] / (2 * numpy.sqrt(2 * numpy.log(2))))
        # print("sigma1 * 4pi= ", 4 * numpy.pi * param['fitresult'] / (2 * numpy.sqrt(2 * numpy.log(2))))
        # print("sigma1 * 4pi * sqrt(2)= ", 4 * numpy.pi * numpy.sqrt(2) * param['fitresult'] / (2 * numpy.sqrt(2 * numpy.log(2))))
        print("sigma1 * sqrt(2) = ", numpy.sqrt(2) * param['fitresult'] / (2 * numpy.sqrt(2 * numpy.log(2))))
        print("emittance = ", 0.69 * numpy.sqrt(2) * param['fitresult'] / (2 * numpy.sqrt(2 * numpy.log(2))))


    dummy_list.append(param['fitresult'])
print("chisq = ", fit.chisq)

y1 = sum_gauss(x, *dummy_list)
# y1 = sum_agauss(x, *dummy_list)

print("Area original: ", y.sum() * (x[1] - x[0]))
print("Area fit: ", y1.sum() * (x[1] - x[0]))


# plot and save file

p = plot(x, y, x, y1,
         xtitle=r'$2 \pi r / \sqrt{ 2 \lambda L }$', ytitle=r'Normalized intensity', show=0)

import matplotlib.pylab as plt
plt.savefig('undulator_size.eps')
plt.show()

####################
import scipy.constants as codata
u_nperiods = 111.111
u_period = 0.018
wavelength = codata.h * codata.c / codata.e / 10000.

cte = 4 * numpy.pi / numpy.sqrt(2 * u_nperiods * u_period * wavelength) # todo: check 4!!!

print("wavelength: ", wavelength)
filename = "id06_backpropagated.dat"
f = open(filename,'w')
for i in range(y.size):
    f.write("%g  %g\n" % (x[i] / cte, y[i]))
f.close()
print("File written to disk: %s" % filename)





