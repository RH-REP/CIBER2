import numpy as np
from scipy import interpolate
from . import spec_reader
import matplotlib.pyplot as plt
from . import angle_function
import pathlib
import math
import csv
import pandas as pd

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
        self.waves_lvf, self.eps_lvf,self.eps_lvf_sterr = self.get_wv_vs_eps_at_lvf()
        self.waves_spec = ls_irradiance_class.waves
        self.ls_irradiance = ls_irradiance_class.irradiance

        _,self.ls_irradiance_reduced_by_ND = self._reduce_by_ND_fiter()
        self.ls_irradiance_reduced_by_ND_ip_to_lvf = interpolate.interp1d(self.waves_spec,self.ls_irradiance_reduced_by_ND,fill_value="extrapolate")(self.waves_lvf)    

        loading = 1.2
        window_transaprent = 0.99**6
        self.cf = self.ls_irradiance_reduced_by_ND_ip_to_lvf * window_transaprent * loading /self.eps_lvf/math.pi 
        self.sterr = self.cf /self.eps_lvf*self.eps_lvf_sterr
    
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
        return self.eps_shelf.getLVF_with_sterr(self.OD,self.arm_name)

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

class NIR_CF_Class():
    
    def __init__(self, reduction_shelf,kelvin) :
        self.reduction_shelf = reduction_shelf
        self.kelvin =kelvin

    def getIrradiance(self,book_name):
        w,irr =calc_ls_radiance(book_name,self.kelvin)
        return  w,irr
    
    def get_wv_vs_eps_at_lvf(self,book_name,arm_name):
        waves_lvf, eps_lvf = self.reduction_shelf.getLVFxy(book_name,arm_name,)
        return waves_lvf, eps_lvf 


    def get_lvf_img(self,book_name,arm_name,):
        img = self.reduction_shelf.get_eps_img(book_name,arm_name,)
        img_at_lvf = img[10:150,]
        return img_at_lvf 


    def calc_CF(self,book_name,arm_name):
        w,ls_radiance = calc_ls_radiance(book_name,self.kelvin)
        waves_lvf, eps_lvf, eps_sterr = self.reduction_shelf.getLVF_with_sterr(book_name,arm_name,)

        ls_radaince_ip = interpolate.interp1d(w,ls_radiance,fill_value="extrapolate")(waves_lvf)
        cf = ls_radaince_ip/eps_lvf
        cf_sterr = ls_radaince_ip/eps_lvf/eps_lvf*eps_sterr
        return waves_lvf,cf, cf_sterr

    # def get_LVF_img(self,book_name,arm_name,):
    #     img = self.get_eps_img(book_name,arm_name,)
    #     img_at_lvf = img[10:150,]
    #     return img_at_lvf 
        # w,irr = calc_ls_radiance(book_name=book_name,kelvin=kelvin)

        # self.waves_for_ir = 
        # self.irradiance = 
        # self.OD = OD
        # self.arm_name = arm_name
        # self.eps_shelf = eps_shelf
        # self.waves_lvf, self.eps_lvf = self.get_wv_vs_eps_at_lvf()
        # self.waves_spec = ls_irradiance_class.waves
        # self.ls_irradiance = ls_irradiance_class.irradiance
        # _,ls_irradiance_reduced_by_ND = self._reduce_by_ND_fiter()
        # ls_irradiance_reduced_by_ND_ip_to_lvf = interpolate.interp1d(self.waves_spec,ls_irradiance_reduced_by_ND,fill_value="extrapolate")(self.waves_lvf)    
        # self.cf = ls_irradiance_reduced_by_ND_ip_to_lvf/self.eps_lvf/math.pi * 0.76*1.2
        


def calc_ls_radiance(book_name,kelvin):
    """
    W/m^2/sr/nm
    """
    waves = np.arange(500,2200,)
    irradiance_from_plank=get_BB_from_planck(waves,kelvin)
    pinholesize = (1.732/2)**2*np.pi
    solid_angle = np.pi*(25.4/2/67.1)**2
    filter_transmittance = get_filter_conbination(book_name)
    ls_radiance_in = irradiance_from_plank*pinholesize*solid_angle*filter_transmittance
    ls_aperture_size = (330/2)**2*np.pi
    # ls_transparent = 0.4
    ls_transparent = get_ls_transparent() * 0.4
    window_transparent = (0.95)**6
    loading = 1.2
    ls_radiance_out = ls_radiance_in/ls_aperture_size/np.pi*ls_transparent*window_transparent*loading

    return waves,ls_radiance_out

def get_filter_conbination(book_name):
    filt1_name = book_name.split("①")[1].split("②")[0]
    filt2_name = book_name.split("②")[1].split("③")[0]
    filt3_name = book_name.split("③")[1]
    w,filt1_trans = get_filter_transmittance(filt1_name)
    w,filt2_trans = get_filter_transmittance(filt2_name)
    w,filt3_trans = get_filter_transmittance(filt3_name)
    
    total_trans_parent = filt1_trans * filt2_trans * filt3_trans
    waves = np.arange(500,2200,)
    total_trans_parent_ip = interpolate.interp1d(w,total_trans_parent,fill_value="extrapolate")(waves)
    return total_trans_parent_ip

def get_filter_transmittance(filtername, filter_path = './lib/NIR_filter/20210716_filter_transmittance.csv'):
    with open(filter_path) as csvfile:
        reader = csv.reader(csvfile)
    #     header = reader
        rawData = [row for row in reader]
    # print(header/)
    header = rawData[0]
    body_T = np.array(rawData[1:]).T
    filter_np = np.array([1])
    for i ,v in enumerate(header):
        if v in filtername:
            filter_np = np.array(body_T[i],dtype=np.float32)
    if len(filter_np) ==1:
        print(filtername ,"doesn't match")
    # print(filtename)
    # print(filter_np)
    w =  np.array(body_T[1],dtype=np.float32)
    return w,filter_np

def get_ls_transparent(file_path = "./lib/sphere_transmittance/Lsphere_relativeT.csv"):
    # file_path = "./lib/sphere_transmittance/Lsphere_relativeT.csv"
    with open(file_path) as csvfile:
            reader = csv.reader(csvfile)
            rawData = [row for row in reader]
    # print(rawData)
    rawData_np = np.array(rawData,dtype=np.float32)
    w = rawData_np.T[0]
    transparent = rawData_np.T[1]
    waves = np.arange(500,2200,)
    transparent_ip = interpolate.interp1d(w,transparent,fill_value="extrapolate")(waves)
    return transparent_ip


# class NIR_CF_Class:
#     def __init__(self, eps_shelf,kelvin,book_name) :
#         self.eps_shelf = eps_shelf
        # NIR_irradiance(kelvin)
        # self.waves_lvf, self.eps_lvf = self.get_wv_vs_eps_at_lvf()
        # self.waves_spec = ls_irradiance_class.waves
        # self.ls_irradiance = ls_irradiance_class.irradiance
        # _,ls_irradiance_reduced_by_ND = self._reduce_by_ND_fiter()
        # ls_irradiance_reduced_by_ND_ip_to_lvf = interpolate.interp1d(self.waves_spec,ls_irradiance_reduced_by_ND,fill_value="extrapolate")(self.waves_lvf)    
        # self.cf = ls_irradiance_reduced_by_ND_ip_to_lvf/self.eps_lvf/math.pi * 0.76*1.2
    
# class NIR_CF:
#     def __init__(self, measuredFiles,bgFiles,baffle,cal_file="2021-06-29.cal") :
#         self.measuredFiles = measuredFiles
#         self.bgFiles = bgFiles
#         self.baffle = baffle
#         self._set_cal(cal_file)
#         self.waves,self.irradiance = self._calc_wave_and_irradiance()

#     def _calc_wave_and_irradiance(self):
#         integratedTime,waves,counts = spec_reader.read_data_set(self.measuredFiles)
#         integratedTime,waves,countsOfBg = spec_reader.read_data_set(self.bgFiles)
#         measuredCounts = counts - countsOfBg
#         angle_ratio = self.baffle.ratio
#         irradiance = get_irradeiance(waves,measuredCounts,integratedTime,self.cal) / angle_ratio
#         return waves,irradiance

#     def show_wave_and_irradiance(self):
#         x,y = self.waves,self.irradiance
#         fig,ax = plt.subplots()
#         ax.scatter(x[100:],y[100:])
#         ax.set_xlabel("wavelength[nm]")
#         ax.set_ylabel("Absolute Irradiance [W/m^2/nm]")
#         plt.plot()
#         return fig, ax
#     def get_wave_and_irradiance(self):
#         return self.waves,self.irradiance 


#     def _set_cal(self,cal_file):
#         _,cal = spec_reader.readCalFacFile(parent_dir + '/spec_lib/cal/' + cal_file)
#         self.cal = cal
#         # self.waves,self.irradiance = self._calc_wave_and_irradiance()



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




def get_BB_from_planck(waves, T):
    """
    input waves as nano
    return irradiance(W/m^2/nm)
    """
    waves = waves *10**-9
    h = 6.6260693 * 10 ** -34
    k = 1.38* 10 ** -23
    c = 2.99792458 * 10 ** 8
    mmtonm = 10**-9
    y = 2*h*c**2/(waves**5)/(np.exp(h*c/k/waves/T)-1)*mmtonm
    
    return y 

if __name__ == "__main__":
    pass
    # home_directory = os.environ['HOME']
    # datadirectory = home_directory + '/spectrometer/SL/20200919/*/'
    # measuredFiles = glob.glob(datadirectory +"3ms*")
    # bgFiles = glob.glob(datadirectory +"BG_3ms*")
    # cal_class = CalibrationClass(measuredFiles,bgFiles)
    # cal_class.saveCF()
    