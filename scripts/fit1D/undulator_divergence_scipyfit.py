import numpy
from scipy.optimize import curve_fit
from scipy.stats import kurtosis, kurtosistest
import matplotlib.pylab as plt



def gauss(x, height, center, sigma):
    return height * numpy.exp(-(x - center)**2 / (2 * sigma**2))

def fitgauss(x, y, p0 = [1.0, 0, 0.5]):
    popt, pcov = curve_fit(gauss, x, y, p0=p0, method='trf')  # ‘lm’, ‘trf’, ‘dogbox’
    y1 = gauss(x, *popt)
    print("FWHM = ", 2.355 * popt[2])
    return popt[2], popt, pcov


def doublegauss(x, height, center, sigma, distance): #, correction):
    out = numpy.exp(-(x - center - distance / 2)**2 / (2 * sigma**2)) + \
          numpy.exp(-(x - center + distance / 2)**2 / (2 * sigma**2)) #+ \
          # correction * numpy.exp(-(x - center) ** 2 / (2 * sigma**2))
    # icen = numpy.argwhere(numpy.abs(x) <= 0.5 * distance)
    # out[icen] = out[out.argmax()]
    return out * height


if __name__ == "__main__":

    ngaussians = 1


    x = numpy.linspace(-4, 4, 1000)
    Gamma = numpy.pi / 2 * x**2
    y = (numpy.sin(Gamma) / Gamma)**2

    import scipy.constants as codata


    # Fitting
    if ngaussians == 1:
        # popt, pcov = curve_fit(gauss, x, y, p0=p0, method='trf') # ‘lm’, ‘trf’, ‘dogbox’
        sd1, popt, pcov = fitgauss(x, y, p0 = [1.0, 0, 0.5])
        y1 = gauss(x, *popt)
        print("FWHM = ", 2.355 * sd1)
    elif ngaussians == 2:
        p0 = [1.0, 0, 0.5, 0.2] #, 0.0]
        popt, pcov = curve_fit(doublegauss, x, y, p0=p0, method='trf') # ‘lm’, ‘trf’, ‘dogbox’
        y1 = doublegauss(x, *popt)
        print("FWHM = ", 2.355 * popt[2] + popt[3])

    print(popt, numpy.diag(pcov))

    # print("Searched parameters = %s" % popt)
    # print("Obtained parameters : ")
    # print("Area original: ", y.sum() * (x[1] - x[0]))
    # print("Area fit: ", y1.sum() * (x[1] - x[0]))


    # plot and save file
    from srxraylib.plot.gol import plot
    p = plot(x, y,
             x, y1,
             xtitle=r'$\theta_r=\theta\sqrt{L/\lambda}$', ytitle=r'$(sin(\Gamma)/\Gamma)^2$', show=0)
    # plt.savefig('undulator_divergence.eps')
    plt.show()


    import scipy.constants as codata
    u_nperiods = 111.111
    u_period = 0.018
    wavelength = codata.h * codata.c / codata.e / 10000.

    cte = numpy.sqrt(u_nperiods * u_period / wavelength)

    print("wavelength: ", wavelength)
    # filename = "id06_farfield.dat"
    # f = open(filename,'w')
    # for i in range(y.size):
    #     f.write("%g  %g\n" % (x[i] * 100 / cte, y[i]))
    # f.close()
    # print("File written to disk: %s" % filename)


    # print("Stdv: ", y.std(), y1.std())
    # print("Kurtosis: ", kurtosis(y, fisher=True), kurtosis(y1))
    # print("KurtosisTest: ", kurtosistest(y), kurtosistest(y1))