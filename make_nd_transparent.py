#NDFilter_transmittance

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import csv
from specReader import specReader

if __name__ == '__main__' :
    od_list = ['OD0','OD1','OD2']
    home_directory = os.environ['HOME']
    dataFile = home_directory + '/spectrometer/filter_measurement/20210115/STL-BF1-AT-BF2-SS-COSsp2-L1-spec'

    cps_set,wavelength = specReader.get_cps_set(dataFile)
    transparent = specReader.get_transparent(cps_set)
    dirname = './lib/filter'
    os.makedirs(dirname,exist_ok=True)
    for od in od_list[1:]:
        specReader.csvWrite(wavelength, transparent[od], dirname, od)
