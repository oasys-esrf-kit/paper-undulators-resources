
from srxraylib.plot.gol import plot
import numpy

def get_stdev(x, y):
    delta = (x[1] - x[0])
    Y = y.copy()
    Y /= y.sum() * delta
    m1 = (x ** 1 * Y).sum() * delta
    m2 = (x ** 2 * Y).sum() * delta
    return numpy.sqrt(m2 - m1**2)

x = numpy.linspace(-100.0, 100, 201)

sigma= 20.0
# y = numpy.exp(-x**2 / 2 / sigma**2) / (sigma * numpy.sqrt(2 * numpy.pi))

y0 = numpy.exp(-x**2 / 2 / sigma**2)
y = y0 / (y0.sum() * (x[1] - x[0]))

# plot(x, y)

m1 = (x**1 * y).sum() * (x[1] - x[0])
m2 = (x**2 * y).sum() * (x[1] - x[0])
m3 = (x**3 * y).sum() * (x[1] - x[0])
m4 = (x**4 * y).sum() * (x[1] - x[0])
print( "1st moment mean: ", m1, 0)
print( "2nd moment variance: ", m2 - m1**2, get_stdev(x, y)**2, sigma**2)
print( "3rd moment Skewness: ", m3)
print( "4th moment Kurtosis: ", m4)