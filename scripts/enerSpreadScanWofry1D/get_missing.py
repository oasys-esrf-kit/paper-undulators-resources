import numpy
import os

def check_missing(i):
    dir = "/scisoft/data/srio/paper-undulator/enerSpreadScanWofry1D/"
    filename = "wfr1D_elec_energy_scan_farfield_n1_shift%i.h5" % i
    return os.path.isfile(dir + filename)


def get_missing():
    Shift = numpy.arange(-300, 501, 1)
    missing = []
    for i in range(Shift.size):
        if not check_missing(Shift[i]):
            missing.append(Shift[i])
    return missing

if __name__ == "__main__":

    print("missing: ", get_missing())
