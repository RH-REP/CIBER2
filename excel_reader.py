import os
from fitsReader import fitsReader


home_directory = os.environ['HOME']
# target_file = "/WSMR2021May/abs_calib/worksheets/worksheet_20210524_NIRAbsCal.xlsx"
season_name = "/WSMR2021May"
experiment_name = "/abs_calib"
worksheet_name = "/worksheets/worksheet_20210521_VisAbsCal.xlsx"
filename = home_directory + season_name + worksheet_name 
print(filename)
book = fitsReader.Get_experiment_book(filename)


for d,v in book._fitDataDict.items():
    print(d)
    print(v._dataSet)
    print(v._bg_dataSet)