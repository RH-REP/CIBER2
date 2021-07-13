import numpy as np
import csv
from scipy import interpolate
import os
import glob
from . import spec_reader
import datetime
import pathlib

p_sub = pathlib.Path(__file__)
parent_dir = str(p_sub.parent)

class CalibrationClass:
    def __init__(self, measuredFiles,bgFiles,) :
        self.measuredFiles = measuredFiles
        self.bgFiles = bgFiles
        w,y = getHLirradiance()
        self.theoretical_irradiance = spec_reader.Wave_and_Y(w,y)
    
    def convertToCF(self):
        integratedTime,waves,counts = spec_reader.read_data_set(self.measuredFiles)
        integratedTime,waves,countsOfBg = spec_reader.read_data_set(self.bgFiles)
        measuredCounts = counts - countsOfBg
        dx = spec_reader.makedx(waves)
        theoretical_irradiance_ip = interpolate.interp1d(self.theoretical_irradiance.wave,self.theoretical_irradiance.y,fill_value="extrapolate")(waves)
        calibrationFactor = theoretical_irradiance_ip * ((0.39/2)**2*np.pi) * dx / measuredCounts * integratedTime 
        return waves, calibrationFactor
    
    def saveCF(self,spec_name,config_name):
        waves, calibrationFactor = self.convertToCF()
        np_to_string = np.array2string(np.array([waves,calibrationFactor]).T,separator="\t",threshold=np.inf)
        # write_string = np_to_string
        write_string = " " + np_to_string.replace("[","").replace("]","").replace("\n  ","\n")[:]
        with open(parent_dir + '/spec_lib/cal/'+spec_name+"/" +str(datetime.date.today()) +"-"+ config_name + '.cal', 'w') as f:
            for _ in range(9):
                f.write("\n")
            f.writelines(write_string)
        return waves, calibrationFactor


def getHLirradiance():
    file = parent_dir + "/spec_lib/HL-Spectrum.csv"
    with open(file) as csvfile:
        reader = csv.reader(csvfile)
        rawData = [ row for row in reader]
    rawData = np.array(rawData,dtype=float)
    waves = rawData.T[0]
    irradiance = rawData.T[1]
    return waves, irradiance

"""
do not deleate
"""
# def getCalRawData(file = "/Users/hashimotoryo/Desktop/CIBER2/Experiment/Spectrometer/SL/19Oct_2020/SL-COSv-L1_2ms_1/2ms_cal_FLMS143431_15-17-51-293.txt"):
#     with open(file) as csvfile:
#         reader = csv.reader(csvfile)
#         rawData = [ row for row in reader]
#         # rawData = [ row.split('\t') for row in reader]
#     rawData = rawData[14:]
#     for i,r in enumerate(rawData) :
#         rawData[i] = r[0].split("\t")
#     rawData = np.array(rawData[14:],dtype=np.float)
#     waves = rawData.T[0]
#     counts = rawData.T[1]
#     return waves, counts


"""
do not deleate
"""
# def readCalFacFile():
#     file = '/Users/hashimotoryo/Desktop/CIBER2/Experiment/Spectrometer/SL/19Oct_2020/SL-COSv-L1_2ms_1/20201019_2ms_1_OOIIrrad.cal'
#     with open(file) as csvfile:
#         reader = csv.reader(csvfile)
#         rawData = [ row for row in reader]
#         # rawData = [ row.split('\t') for row in reader]
#     rawData = rawData[9:]
#     # print(rawData)
#     for i,r in enumerate(rawData) :
#         rawData[i] = r[0].split("\t")
#     rawData = np.array(rawData,dtype=np.float)
#     waves = rawData.T[0]
#     calibrationFactorBySoft = rawData.T[1]
#     return waves, calibrationFactorBySoft

if __name__ == "__main__":
    home_directory = os.environ['HOME']
    datadirectory = home_directory + '/spectrometer/SL/20210405/L1-glass/'
    measuredFiles = glob.glob(datadirectory +"2ms*")
    print(measuredFiles)
    bgFiles = glob.glob(datadirectory +"BG_2ms*")
    cal_class = CalibrationClass(measuredFiles,bgFiles)
    cal_class.saveCF()
    