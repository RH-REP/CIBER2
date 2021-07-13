import matplotlib
import astropy.io.fits as fits
import numpy as np
import glob
import os
from memory_profiler import profile

class AreaSizeClass:
    def __init__(self,xs,xe,ys,ye):
        self.xs = xs
        self.xe = xe
        self.ys = ys
        self.ye = ye

def loadFits(fileName):
    if os.path.getsize(fileName) < 16781760:
        print(fileName," is broken")
        print(np.full((2048,2048,), None))
        return np.full((2048,2048,), None),None
    hdulist=fits.open(fileName)
    # hdu=hdulist[0]
    data_np = np.full((len(hdulist),2048,2048),0,dtype=np.float32)
    for i,hdu in enumerate(hdulist):
        data=np.array(hdu.data / (4*10**(-6)),dtype=np.float32)
        data_np[i] = data
    header=hdu.header
    return data_np, header

def saveFits(data,fileName,*data2):
    hdu = fits.PrimaryHDU(np.array(data*(4*10**(-6)),dtype=np.float32))
    hdu2_list = [fits.ImageHDU(d) for d in data2]
    strage_list = [hdu]+ hdu2_list
    hdulist = fits.HDUList(strage_list)
    hdulist.writeto(fileName,overwrite=True)

def getFrameNumber(filename):
    return int(filename.split('/FRM')[-1].split('_RFR')[0].split('_PIX')[0])

def convertNumberToSecond(number):
    frameRate = 1.5
    return number * frameRate

def convertWorkingdirectoryToSecond(working_directory):
    return convertNumberToSecond(getFrameNumber(getSortedPixFiles(working_directory[0])[0]))


def zero_padding_sort(unsorted_files):
    unsorted_dict ={}

    for filename in unsorted_files:
        unsorted_dict[filename] = getFrameNumber(filename)
    sorted_list = sorted(unsorted_dict.items(), key=lambda x:x[1])
    sorted_files = []
    for v,_ in sorted_list:
        sorted_files.append(v)
    return sorted_files

def getSortedPixFiles(folder_name):
    data_before_sort = glob.glob(folder_name + '/*_PIX.fts')
    files = zero_padding_sort(data_before_sort)
    return files

# @profile
def getTimeAndImage(folder_name):
    files = getSortedPixFiles(folder_name)
    if len(files) == 0:
        print("no fits data. please check ",folder_name)
        return [],[]
    reset_fits_file = files[0]
    reset_number = getFrameNumber(reset_fits_file)
    file_length = len(files)
    times = np.full(file_length,np.nan)
    data_images = np.full((file_length,2048,2048,),np.nan,dtype=np.float32)
    ##do not deleate##
    # print('file length is ',file_length)
    # print('memory size is ','{:.2g}'.format(data_images.__sizeof__()/2**30),"GB")

    for i,f in enumerate(files):
        frame_number = getFrameNumber(f)
        times[i]= convertNumberToSecond(frame_number-reset_number)
        data_images[i],_ = loadFits(f)
        if i+1 == file_length:
            break
    return times[:-1], data_images[:-1]

###     Do NOT deleate    ###
# def getAverageImage(data_set,):
#     data_length = len(data_set)
#     area_average_set = np.arange(data_length)
#     area_average_set_error = np.arange(data_length)
#     for i,data in enumerate(data_set):
#         # area_average_set[i],area_average_set_error[i] = getAreaAverage(data)
#         # area_average_set[i],area_average_set_error[i] = np.average(data,axis=0),np.average(data,axis=0)
#         area_average_set[i],area_average_set_error[i] = np.nanmean(data[cutArea.ys:cutArea.ye,cutArea.xs:cutArea.xe]),np.nanstd(data[cutArea.ys:cutArea.ye,cutArea.xs:cutArea.xe])
#         # area_average_set[i],area_average_set_error[i] 
    
#     return area_average_set,area_average_set_error

def getFitsAfterReset(folders):
    frames_after_reset = []
    for f in folders:
        rfrs = glob.glob(f + "/*RFR.fts")
        if len(rfrs) == 0:continue
            
        reset_number = rfrs[0].split("/FRM")[-1].split("_RFR")[0]
        frame_after_reset = glob.glob(f + "/*" + str(int(reset_number)+1) + "*PIX.fts")[0]
        frames_after_reset.append(frame_after_reset)
    return frames_after_reset
    # data_set = np.full((len(frames_after_reset),2048,2048,),np.nan)

def getResetFits(folders):
    frames = []
    for f in folders:
        rfrs = glob.glob(f + "/*RFR.fts")
        if len(rfrs) == 0:continue
            
        frames.append(rfrs[0])
    return frames
    # data_set = np.full((len(frames_after_reset),2048,2048,),np.nan)
def getResetFitsFolder(folders):
    folders_contain_reset = []
    for f in folders:
        rfrs = glob.glob(f + "/*RFR.fts")
        if len(rfrs) == 0:continue
        folders_contain_reset.append(f)
    return folders_contain_reset
    # data_set = np.full((len(frames_after_reset),2048,2048,),np.nan)
    return

def getFitsDataSet(fits_files):
    data_set = np.full((len(fits_files),2048,2048,),np.nan)
    for i, f in enumerate(fits_files):
        data,_ = loadFits(f)
        data_set[i] = data
    return data_set


def loadMask(ch,):
    f = "./lib/mask/" + ch + "mask.fts"
    masks, header = loadFits(f)
    mask = np.nanmean(masks, axis=0)
    return mask,header


def loadRawMask(ch,):
    f = "./lib/mask/" + ch + "mask.fts"
    raw_mask, header = loadFits(f)
    return raw_mask,header



def save_mask(mask_data,index,ch):
    masks, header = loadRawMask(ch)
    if len(masks) < index:
        masks = np.append(masks,mask_data).reshape((index,2048,2048))
    else:
        masks[index] = mask_data

    f = "./lib/mask/" + ch +"mask.fts"
    saveFits(masks,f)
    return masks, header
