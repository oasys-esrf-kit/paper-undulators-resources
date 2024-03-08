import numpy
import os

def check_missing(i):
    dir = "/home/esrf/srio/paper-undulators-resources/scripts/enerSpreadScanWofry1D/"
    filename = "wfr1D_elec_energy_scan_farfield_n1_shift%i.h5" % i
    print(i, os.path.isfile(dir + filename))
    return os.path.isfile(dir + filename)

if __name__ == "__main__":

    Shift = numpy.arange(-300, 601, 1)
    missing = []
    for i in range(Shift.size):
        if check_missing(i):
            missing.append(i)

    print("missing: ", missing)