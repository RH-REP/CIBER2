import spec_calibration
import os
import glob
import matplotlib.pyplot as plt
import angle_function
# print(angle_function.calcAngleRatio())
# half_angel,sensitivity = angle_function.getBaffleF()
# plt.ion()
# plt.scatter(half_angel,sensitivity)
# plt.show()
# input()

home_directory = os.environ['HOME']
datadirectory = home_directory + '/spectrometer/calibration_data/20210118/cosSp2-F1/0deg/'
measuredFiles = glob.glob(datadirectory +"3ms*")
bgFiles = glob.glob(datadirectory +"BG_3ms*")
test_data = spec_calibration.IrradianceClass(measuredFiles,bgFiles,angle_csv='matsumi40deg.csv')
plt.ion()
test_data.show_wave_and_irradiance()
input()

