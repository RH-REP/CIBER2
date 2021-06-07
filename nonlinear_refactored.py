import matplotlib.pyplot as plt
import glob
from fitsReader import fitsReader
import os

home_directory = os.environ['HOME']
target_directory = "/WSMR2021May/abs_calib/vis/frames/ch3/RUN05222021_115228"
ch3 = sorted(glob.glob(home_directory + target_directory + '/*_PIX.fts'))

times, data_image_minus_reset = fitsReader.getTimeAndImage(ch3, ch3[0])
cutArea = fitsReader.AreaSizeClass(1100,1200,1400,1500)
average_of_left_area,error_right_area_average =fitsReader.getAverage(data_image_minus_reset,cutArea)

plt.plot(times, average_of_left_area, marker='o', lw=0)
plt.ylabel('electrons')
plt.xlabel('time[s]')
plt.title('RUN05222021_115228, Arm-S')
plt.show()
plt.savefig('./fig/Test_3.png')