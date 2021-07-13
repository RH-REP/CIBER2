import numpy as np
import csv
from scipy import interpolate
import os
import glob
import sys
import pathlib

p_sub = pathlib.Path(__file__)
parent_dir = str(p_sub.parent)

class Wave_and_Y:
    def __init__(self, wave,y) :
        self.wave = wave
        self.y = y

od_list = ['OD0','OD1','OD2']

def read_data_set(files):
    if files == []:
        print("no files")
        # sys.exit()
    y_np = np.full((len(files),2048),np.nan)
    for i,f in enumerate(files):
        t,x,y = read_data(f)
        y_np[i] = y
    
    return t,x,np.average(y_np,axis=0)

def csv_to_np(file,isHeaderExist=False):
    with open(file) as csvfile:
        reader = csv.reader(csvfile)
        if isHeaderExist:
            next(reader)
        rawData = [np.array(list(map(float, row))) for row in reader]

    return np.array(rawData)

def read_data(fname):
    with open(fname, "r") as load_data:
        lines = load_data.readlines()
        data = lines[14:]
        data_np = np.full((len(data),2),np.nan)
        for i,d in enumerate(data):
            data_np[i] = np.array(d.split('\t'),dtype=np.float32)
        x = data_np.T[0]
        y = data_np.T[1]
        times = np.array(lines[6][24:],dtype=np.float32)
    times = np.array(lines[6][24:],dtype=np.float32)
    return times, x, y

def correctNonlinearity(counts):
    """
    correct Counts from 2020 11 data
    cosvis - L1 -spec
    """
    newFitParam2 = np.array([ 
    1.91304720e-33,-4.33364502e-28,
    3.93044161e-23,-1.86963569e-18,
    5.13784287e-14,-8.96437198e-10,
    1.18028149e-05,8.97667951e-01
    ])
    correctFunction = np.poly1d(newFitParam2)(counts)
    correctedData = counts/correctFunction
    return correctedData

def csvWrite(x, y, dirname, name):
    list_data = [x, y]
    list_data = np.array(list_data).T
    with open(dirname + '/' + name + '.csv', 'w') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(['wavelength[nm]', 'transmittance']) 
        writer.writerows(list_data)

def get_cps(dataFile,od_name):
    files = sorted(glob.glob(dataFile + '/' + od_name + '*.txt', recursive=True))
    bg_files = sorted(glob.glob(dataFile + '/BG_*.txt'))
    wavelength = read_data(files[0])[1]
    count_np = np.full((len(files),len(wavelength)),np.nan)
    for i in range(len(files)):
        times,_,data = read_data(files[i])
        _,_,data_bg = read_data(bg_files[i])
        count_np[i] = np.array(data) - np.array(data_bg)
    count_ave = np.mean(count_np, axis=0)
    cps= count_ave/times
    return cps,wavelength

def get_cps_set(dataFile):
    cps_set={}
    for od in od_list:
        cps_set[od],wavelenght = get_cps(dataFile,od)
    return cps_set,wavelenght

def makedx(waves):
    dx = []
    for i in range(len(waves)-1):
        dx.append(waves[i+1] - waves[i])
    dx.append(dx[-1])
    dx = np.array(dx,dtype=float)
    return dx

def calTransparent(cps_set,main_od):
    return cps_set[main_od]/cps_set['OD0']

def get_transparent(cps_set):
    transparent = {}
    for main_od in od_list[1:]:
        transparent[main_od] = calTransparent(cps_set,main_od)
    return transparent

def readCalFacFile(datafile):
    # datafile = parent_dir + "/spec_lib/cal/2021-06-29.cal"
    with open(datafile) as csvfile:
        reader = csv.reader(csvfile)
        rawData = [ row for row in reader]
    rawData = rawData[9:]
    for i,r in enumerate(rawData) :
        rawData[i] = np.array(r[0].split("\t"),np.float32)
    rawData = np.array(rawData,dtype=np.float32)
    waves = rawData.T[0]
    calibrationFactorBySoft = rawData.T[1]
    return waves, calibrationFactorBySoft


    
