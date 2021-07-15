import numpy as np
from scipy import interpolate
from . import spec_reader
import matplotlib.pyplot as plt
from . import angle_function
import pathlib
import math

p_sub = pathlib.Path(__file__)
parent_dir = str(p_sub.parent)
parent_dir = '.'

class IrradianceClass:
    def __init__(self, measuredFiles,bgFiles,baffle,cal_file="2021-06-29.cal") :
        self.measuredFiles = measuredFiles
        self.bgFiles = bgFiles
        self.baffle = baffle
        self._set_cal(cal_file)
        self.waves,self.irradiance = self._calc_wave_and_irradiance()

    def _calc_wave_and_irradiance(self):
        integratedTime,waves,counts = spec_reader.read_data_set(self.measuredFiles)
        integratedTime,waves,countsOfBg = spec_reader.read_data_set(self.bgFiles)
        measuredCounts = counts - countsOfBg
        angle_ratio = self.baffle.ratio
        irradiance = get_irradeiance(waves,measuredCounts,integratedTime,self.cal) / angle_ratio
        return waves,irradiance

    def show_wave_and_irradiance(self):
        x,y = self.waves,self.irradiance
        fig,ax = plt.subplots()
        ax.scatter(x[100:],y[100:])
        ax.set_xlabel("wavelength[nm]")
        ax.set_ylabel("Absolute Irradiance [W/m^2/nm]")
        plt.plot()
        return fig, ax
    def get_wave_and_irradiance(self):
        return self.waves,self.irradiance 


    def _set_cal(self,cal_file):
        _,cal = spec_reader.readCalFacFile(parent_dir + '/spec_lib/cal/' + cal_file)
        self.cal = cal
        # self.waves,self.irradiance = self._calc_wave_and_irradiance()


class CF_Class:
    def __init__(self, eps_shelf,ls_irradiance_class,OD,arm_name) :
        self.OD = OD
        self.arm_name = arm_name
        self.eps_shelf = eps_shelf
        self.waves_lvf, self.eps_lvf = self.get_wv_vs_eps_at_lvf()
        self.waves_spec = ls_irradiance_class.waves
        self.ls_irradiance = ls_irradiance_class.irradiance
        _,self.ls_irradiance_reduced = self._reduce_by_ND_fiter()
        self.ls_irradiance_reduced_ip_to_lvf = interpolate.interp1d(self.waves_spec,self.ls_irradiance_reduced,fill_value="extrapolate")(self.waves_lvf)    
        self.cf = self.ls_irradiance_reduced_ip_to_lvf/self.eps_lvf/math.pi * 0.76*1.2
    
    def _reduce_by_ND_fiter(self):
        _,transparent = self.get_ND()
        irradiance_spec_reduced = self.ls_irradiance * transparent
        return _,irradiance_spec_reduced

    def get_ND(self):
        file = './lib/filter/'+self.OD+'.csv'
        rawData_np = spec_reader.csv_to_np(file,isHeaderExist=True)
        w,transparent = rawData_np.T[0],rawData_np.T[1]
        return w,transparent
        
    def get_wv_vs_eps_at_lvf(self):
        return self.eps_shelf.getLVFxy(self.OD,self.arm_name)

    def get_eps_img(self):
        return self.eps_shelf.get_eps_img(self.OD,self.arm_name)

    def get_LVF_img(self):
        return self.eps_shelf.get_LVF_img(self.OD,self.arm_name)
    # def show_wv_vs_irradiance(self):
    #     fig,ax = plt.subplots()

    #     ax.scatter(self.waves_spec,self.ls_irradinance_reduced,label= self.OD+" reduced")
    #     ax.legend()
    #     ax.set_xlabel("wavelength [nm]")
    #     ax.set_ylabel("abs irradiance [W/m^2/nm]")
    #     ax.set_ylim(1*10**-7,0.01)
    #     ax.set_yscale("log")

    # def show_wv_vs_cf(self):
    #     fig,ax = plt.subplots()
    #     ax.scatter(self.waves_lvf,self.cf,label="measured")
    #     ax.legend()
    #     ax.set_xlim(500,1000)
    #     ax.set_ylim(10**-10,10**-8)
    #     ax.set_yscale("log")
    #     ax.set_xlabel("wavelength [nm]")
    #     ax.set_ylabel("CF[radiance/eps]")
    #     return 


def get_cf_from_theory():
    f = parent_dir + '/spec_lib/cal/20210629_cf.csv'
    raw = spec_reader.csv_to_np(f,isHeaderExist=True)
    w,ti = raw.T[0],raw.T[1]    
    return w,ti
# def ip_easy(x1,y1,x2):
#     y2 = interpolate.interp1d(x1,y1,fill_value="extrapolate")(x2)
#     return y2

def get_irradeiance(waves,measured_counts,integratedTime,cf):
    dx = spec_reader.makedx(waves)
    uwcmtowm = 10**(-6+4)
    return cf /(((0.39/2)**2*np.pi) * dx * integratedTime) * measured_counts*uwcmtowm
    # calibrationFactor = theoretical_irradiance_ip * ((0.39/2)**2*np.pi) * dx / measuredCounts * integratedTime 


if __name__ == "__main__":
    pass
    # home_directory = os.environ['HOME']
    # datadirectory = home_directory + '/spectrometer/SL/20200919/*/'
    # measuredFiles = glob.glob(datadirectory +"3ms*")
    # bgFiles = glob.glob(datadirectory +"BG_3ms*")
    # cal_class = CalibrationClass(measuredFiles,bgFiles)
    # cal_class.saveCF()
    