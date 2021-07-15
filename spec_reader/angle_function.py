import csv
import numpy as np
from scipy import interpolate
import pathlib
import matplotlib.pyplot as plt
import math 
p_sub = pathlib.Path(__file__)
parent_dir = str(p_sub.parent)
parent_dir = "."

class BaffleClass:
    def __init__(self, angle_file = "./cosM.csv",angle_base_file = "./cosM.csv",) :
        self.angle_base_file = angle_base_file
        self.angle_file = angle_file
        self.angle, self.sensitivity = self._get_baffle_function(angle_file)
        self.angle, self.sensitivity_base = self._get_baffle_function(angle_base_file)
        self.deleta_sense,self.deleta_sense_base = self._calc_delta_sense()
        self.ratio = np.sum(self.deleta_sense)/np.sum(self.deleta_sense_base)

    def _get_baffle_function(self,file):

        lib = parent_dir + "/spec_lib/angle_vs_sensitivity/"
        full_angle, sensitivity = get_rad_vs_sensitivity_from_csv(lib + file)
        centerIndex = int(len(full_angle)/2)
        left_sensitivity = sensitivity[:centerIndex+1][::-1]
        right_sensitivity = sensitivity[centerIndex:]
        sensitivity_average = (left_sensitivity + right_sensitivity)/2
        half_angle = full_angle[centerIndex:]
        half_angle,sensitivity_average = extend_to_90deg(half_angle,sensitivity_average)
        return half_angle, sensitivity_average

    def _calc_delta_sense(self):
        sphere_angle, spherea_area_np =culc_angle_vs_sphere_area()
        sphere_area_np_ip= interpolate.interp1d(sphere_angle,spherea_area_np)(self.angle)  
        delta_sense =  sphere_area_np_ip * self.sensitivity
        delta_sense_base =  sphere_area_np_ip * self.sensitivity_base
        return delta_sense,delta_sense_base

    def show_ratio(self):
        fig,ax = plt.subplots()
        ax.scatter(self.angle,self.sensitivity,label="w/  baffle")
        ax.scatter(self.angle,self.sensitivity_base,label="w/o baffle")
        ax.set_xlabel("angle [deg]")
        ax.set_ylabel("angle sensitivity R [a.u.]")
        ax.legend()
        return fig,ax

    def show_delta_sense(self):
        fig,ax = plt.subplots()
        x,y = self.angle,self.deleta_sense
        ax.scatter(x,y,label="baffle")
        x_b,y_b = self.angle,self.deleta_sense_base
        ax.scatter(x_b,y_b,label="base")
        ax.legend()
        ax.set_xlabel("angle [deg]")
        ax.set_ylabel("delta sensitivity /angle")
        return fig,ax

    def show_sphere_area(self):
        fig,ax = plt.subplots()
        x,y = culc_angle_vs_sphere_area()
        print("sum is ",np.sum(y))
        ax.scatter(x,y)
        ax.set_xlabel("angle [deg]")
        ax.set_ylabel("delta s")
        return fig,ax

def culc_angle_vs_sphere_area():
    size = 9000
    sphere_area_np = np.full(size + 1,0,dtype=np.float32)
    angle = np.linspace(0,90,size + 1,)
    for i,ang in enumerate(angle):
        ds = 2*math.pi*math.sin(math.radians(ang)) * math.radians(90/size)/(2*math.pi)
        sphere_area_np[i] = ds
    return angle,sphere_area_np

def get_rad_vs_sensitivity_from_csv(file):
    with open(file) as csvfile:
        reader = csv.reader(csvfile)
        rawData = [np.array(list(map(float, row))) for row in reader]
    rawData_np = np.array(rawData)
    full_angle = rawData_np.T[0]
    sensitivity = rawData_np.T[1]
    return full_angle, sensitivity

def extend_to_90deg(half_angle,sensitivity_average):
    maxhalf_Angle = np.max(half_angle) 
    while maxhalf_Angle < 90:
        space = half_angle[1] - half_angle[0]
        half_angle = np.append(half_angle,maxhalf_Angle + space)
        sensitivity_average = np.append(sensitivity_average,0)
        maxhalf_Angle = np.max(half_angle[:])
    return half_angle,sensitivity_average
