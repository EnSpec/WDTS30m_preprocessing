#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Author: Zhiwei Ye
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

#Step0:
#python aws_down_topo.py --inputfiles /home/ye6/hys_test/gdrive_test/new_obs_ort/f130522t01p00r07_s2sub_obs_v1a_img_geocorr_norot --outputdir '/home/ye6/hys_test/gdrive_test/aws_dsm'

# Step1:
# python /home/ye6/hys_test/gdrive_test/northup_2_rot_v4.py -img /home/ye6/hys_test/gdrive_test/aws_dsm/topo_f130522t01p00r07_tmp_terrain.tif   -ref /home/ye6/hys_test/gdrive_test/new_obs_ort/f130522t01p00r07_s2sub_obs_v1a_img_geocorr_norot -od /home/ye6/hys_test/gdrive_test/aws_dsm/

# Step2:
# python merge_aws_obs_ort.py -img /home/ye6/hys_test/gdrive_test/aws_dsm/topo_f130522t01p00r07_tmp_terrain_rot.tif   -ref /home/ye6/hys_test/gdrive_test/new_obs_ort/f130522t01p00r07_s2sub_obs_v1a_img_geocorr_norot  -od /home/ye6/hys_test/gdrive_test/aws_dsm/

import argparse, os, sys
import numpy as np
from shutil import copyfile

import logging
logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.DEBUG)

try:    
    from osgeo import gdal
except:
    import gdal
   

# copy from HyTools
dtypeDict = {1:np.uint8,
                   2:np.int16,
                   3:np.int32,
                   4:np.float32,
                   5:np.float64,
                   12:np.uint16,
                   13:np.uint32,
                   14:np.int64,
                   15:np.uint64}
                   
# copy from HyTools
def parse_ENVI_header(hdrFile):
    """Parse ENVI header into dictionary
    """

    # Dictionary of all types
    fieldDict = {"acquisition time": "str",
                 "band names":"list_str", 
                 "bands": "int", 
                 "bbl": "list_float",
                 "byte order": "int",
                 "class lookup": "str",
                 "class names": "str",
                 "classes": "int",
                 "cloud cover": "float",
                 "complex function": "str",
                 "coordinate system string": "str",
                 "correction factors": "list_float",
                 "data gain values": "list_float",
                 "data ignore value": "float",
                 "data offset values": "list_float",
                 "data reflectance gain values": "list_float",
                 "data reflectance offset values": "list_float",
                 "data type": "int",
                 "default bands": "list_float",
                 "default stretch": "str",
                 "dem band": "int",
                 "dem file": "str",
                 "description": "str",
                 "envi description":"str",
                 "file type": "str",
                 "fwhm": "list_float",
                 "geo points": "list_float",
                 "header offset": "int",
                 "interleave": "str",
                 "lines": "int",
                 "map info": "list_str",
                 "pixel size": "list_str",
                 "projection info": "str",
                 "read procedures": "str",
                 "reflectance scale factor": "float",
                 "rpc info": "str",
                 "samples":"int",
                 "security tag": "str",
                 "sensor type": "str",
                 "smoothing factors": "list_float",
                 "solar irradiance": "float",
                 "spectra names": "list_str",
                 "sun azimuth": "float",
                 "sun elevation": "float",
                 "wavelength": "list_float",
                 "wavelength units": "str",
                 "x start": "float",
                 "y start": "float",
                 "z plot average": "str",
                 "z plot range": "str",
                 "z plot titles": "str"}

    headerDict = {}

    headerFile = open(hdrFile,'r')
    line = headerFile.readline()
      
    while line :
        if "=" in line:
            key,value = line.rstrip().split("=",1)
            # Add field not in ENVI default list
            if key.strip() not in fieldDict.keys():
                fieldDict[key.strip()] = "str"
            
            valType = fieldDict[key.strip()]
            
            if "{" in value and not "}" in value: 
                while "}" not in line:
                    line = headerFile.readline()
                    value+=line

            if '{}' in value: 
                value = np.nan
            elif valType == "list_float":
                value= np.array([float(x) for x in value.translate(str.maketrans("\n{}","   ")).split(",")])
            elif valType == "list_int":
                value= np.array([int(x) for x in value.translate(str.maketrans("\n{}","   ")).split(",")])
            elif valType == "list_str":
                value= [x.strip() for x in value.translate(str.maketrans("\n{}","   ")).split(",")]
            elif valType == "int":
                value = int(value.translate(str.maketrans("\n{}","   ")))
            elif valType == "float":
                value = float(value.translate(str.maketrans("\n{}","   ")))
            elif valType == "str":
                value = value.translate(str.maketrans("\n{}","   ")).strip().lower()

            headerDict[key.strip()] = value
                            
        line = headerFile.readline()
    
    # Fill unused fields with nans
    for key in fieldDict.keys():
        if key not in headerDict.keys():
            headerDict[key] = np.nan
    
    headerFile.close()
    return headerDict


def copy_hdr(inimg, outimg):
    copyfile(inimg+'.hdr', outimg+'.hdr')

def copy_binary(inimg, outimg):
    copyfile(inimg, outimg)
    
def update_slope_aspect(slope, aspect, mask):
        #mask = slope>no_data_value
        slope[mask] = 90 -  slope[mask]
        aspect[mask] = ((aspect[mask] - 90) % 360)
        #aspect[aspect > 180] = aspect[aspect > 180] - 360
    
# input angles are in degrees    
def cal_cosine_i(slope, aspect, to_sun_zenith, to_sun_azimuth, mask, no_data_value):
    out_cos_i =  np.cos(np.radians(slope))*np.cos(np.radians(to_sun_zenith))+np.sin(np.radians(slope))*np.sin(np.radians(to_sun_zenith))*np.cos(np.radians(aspect - to_sun_azimuth))
    out_cos_i[~mask] = no_data_value
    return out_cos_i


def merge_replace_slope_aspect(Slope,Aspect, mask, new_slope, new_aspect, new_mask):
    Slope[mask & new_mask] = new_slope[mask & new_mask]
    Aspect[mask & new_mask] = new_aspect[mask & new_mask]
    
def main():
    
    parser = argparse.ArgumentParser(description = "Merge tool.")
    parser.add_argument("-img", help="Input new topo image pathname",required=True, type = str)
    parser.add_argument("-od", help="Output directory for all resulting products", required=True, type = str)
    parser.add_argument("-ref", help="Reference old topo image pathname",required=True, type = str)
    parser.add_argument("-corr", help="Execute correction on old version of obs_ort",required=False, default=False, type = bool)    
    #parser.add_argument("-b", help="Band to match", nargs='+', required=False, default=[7,8],type = int)


    args = parser.parse_args()
    
    infile = args.img
    refimg = args.ref
    outpath = args.od
    do_corr = args.corr

    logging.info("Input {} Reference {}".format(infile, refimg))    

    outimg = outpath + '/' + os.path.basename(refimg)+'_corr'
    copy_hdr(refimg, outimg)
    copy_binary(refimg, outimg)    

    headerDict = parse_ENVI_header(outimg+'.hdr')    
    no_data_value = headerDict["data ignore value"]
    
    ds_new = gdal.Open(infile)
    
    #print(ds_new.RasterYSize ,headerDict['lines'], ds_new.RasterYSize != headerDict['lines'])
    #print(ds_new.RasterXSize , headerDict['samples'],ds_new.RasterXSize != headerDict['samples'])
    
    if (ds_new.RasterYSize != headerDict['lines']) or (ds_new.RasterXSize != headerDict['samples']):
        logging.warning('Dimensions do not match. Exit.')
        quit()
    #quit()
    
    new_slope = ds_new.GetRasterBand(1).ReadAsArray()
    #print(new_slope.min(),new_slope.max())


    new_aspect = ds_new.GetRasterBand(2).ReadAsArray()
    #print(new_aspect.min(),new_aspect.max())
    
    #print((new_slope>0).dtype)
    #new_mask = ~ ((new_slope==0) & (new_aspect==0))
    new_mask = new_slope>0.5*(-9999)  # ~ ((new_slope==-9999) & (new_aspect==-9999))
    #print(np.count_nonzero(new_mask))
    #print(np.count_nonzero(new_slope>0))
    #print(np.count_nonzero(new_aspect>0))
    ds_new=None
    #print(np.count_nonzero(new_mask))
    #print(np.count_nonzero(new_slope))
    #print(np.count_nonzero(new_aspect))
    #print(np.count_nonzero(new_aspect>-9999))    
    #new_mask.astype(np.uint8).tofile('/home/ye6/hys_test/new_mmmm.bin')
    #quit()
    
    if headerDict["interleave"]=='bip' or headerDict["interleave"]=='BIP':
        image_mammap = np.memmap(outimg, dtype = dtypeDict[headerDict["data type"]], mode='r+', shape = (headerDict['lines'],headerDict['samples'], int(headerDict['bands'])), offset=0)        
        to_sun_azimuth = image_mammap[:,:,3]
        to_sun_zenith = image_mammap[:,:,4]
        Slope = np.copy(image_mammap[:,:,6])
        Aspect = np.copy(image_mammap[:,:,7])
        
        mask = Slope>no_data_value

        #mask.astype(np.uint8).tofile('/home/ye6/hys_test/old_mmmm.bin')
        #new_mask.astype(np.uint8).tofile('/home/ye6/hys_test/new_mmmm.bin')
        #quit()
        
        if do_corr:
            update_slope_aspect(Slope, Aspect, mask)
            
        merge_replace_slope_aspect(Slope,Aspect, mask, new_slope, new_aspect, new_mask)
            
        cos_i = cal_cosine_i(Slope, Aspect,  to_sun_zenith, to_sun_azimuth, mask, no_data_value)        
        
        image_mammap[:,:,6] = Slope
        image_mammap[:,:,7] = Aspect
        image_mammap[:,:,8] = cos_i        

    elif headerDict["interleave"]=='bil' or headerDict["interleave"]=='BIL':
        image_mammap = np.memmap(outimg, dtype = dtypeDict[headerDict["data type"]], mode='r+', shape = (headerDict['lines'], int(headerDict['bands']),headerDict['samples']), offset=0)
        to_sun_azimuth = image_mammap[:,3,:]
        to_sun_zenith = image_mammap[:,4,:]
        Slope = np.copy(image_mammap[:,6,:])
        Aspect = np.copy(image_mammap[:,7,:])
        
        mask = Slope>no_data_value


        if do_corr:
            update_slope_aspect(Slope, Aspect, mask)
            
        merge_replace_slope_aspect(Slope,Aspect, mask, new_slope, new_aspect, new_mask)
        cos_i = cal_cosine_i(Slope, Aspect,  to_sun_zenith, to_sun_azimuth, mask, no_data_value)

        image_mammap[:,6,:] = Slope
        image_mammap[:,7,:] = Aspect
        image_mammap[:,8,:] = cos_i
        
    else:
    # 'BSQ'
        image_mammap = np.memmap(outimg, dtype = dtypeDict[headerDict["data type"]], mode='r+', shape = (int(headerDict['bands']),headerDict['lines'], headerDict['samples']), offset=0)    
        to_sun_azimuth = image_mammap[3,:,:]
        to_sun_zenith = image_mammap[4,:,:]
        Slope = np.copy(image_mammap[6,:,:])
        Aspect = np.copy(image_mammap[7,:,:])
        
        mask = Slope>no_data_value
        
        if do_corr:
            update_slope_aspect(Slope, Aspect, mask)
        
        merge_replace_slope_aspect(Slope,Aspect, mask, new_slope, new_aspect, new_mask)
        cos_i = cal_cosine_i(Slope, Aspect,  to_sun_zenith, to_sun_azimuth, mask, no_data_value)
        
        image_mammap[6,:,:] = Slope
        image_mammap[7,:,:] = Aspect
        image_mammap[8,:,:] = cos_i
        
    del image_mammap
        
if __name__== "__main__":
    main()    