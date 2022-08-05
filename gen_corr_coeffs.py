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

# python /home/ye6/hys_test/gdrive_test/gen_corr_coeffs.py /home/ye6/hys_test/gdrive_test/new_obs_ort/f130522t01p00r07_s2sub_rfl_v1a_img /home/ye6/hys_test/gdrive_test/temp/

import os, sys, json
import numpy as np
from collections import Counter

#import matplotlib.pyplot as plt

import logging
logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.DEBUG)


bad_band = [[300,400],[1320,1430],[1800,1960],[2450,2600]]

#from HyTools
# ENVI datatype conversion dictionary
dtype_dict = {1:np.uint8,
             2:np.int16,
             3:np.int32,
             4:np.float32,
             5:np.float64,
             12:np.uint16,
             13:np.uint32,
             14:np.int64,
             15:np.uint64}

# Dictionary of all ENVI header fields
field_dict = {"acquisition time": "str",
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
              "correction_factors": "list_float", 
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

def envi_header_dict():
    """
    Returns:
        dict: Empty ENVI header dictionary.
    """
    return {key:None for (key,value) in field_dict.items()}

def parse_envi_header(header_file):
    """
    Args:
        header_file (str): Header file pathname.
    Returns:
        dict: Populated header dictionary.
    """

    header_dict = envi_header_dict()
    header_file = open(header_file,'r')
    line = header_file.readline()

    while line :
        if "=" in line:
            key,value = line.rstrip().split("=",1)
            # Add fields not in ENVI default list
            if key.strip() not in field_dict.keys():
                field_dict[key.strip()] = "str"
            val_type = field_dict[key.strip()]

            if "{" in value and not "}" in value:
                while "}" not in line:
                    line = header_file.readline()
                    value+=line

            if '{}' in value:
                value = None
            elif val_type == "list_float":
                value= np.array([float(x) for x in value.translate(str.maketrans("\n{}","   ")).split(",")])
            elif val_type == "list_int":
                value= np.array([int(x) for x in value.translate(str.maketrans("\n{}","   ")).split(",")])
            elif val_type == "list_str":
                value= [x.strip() for x in value.translate(str.maketrans("\n{}","   ")).split(",")]
            elif val_type == "int":
                value = int(value.translate(str.maketrans("\n{}","   ")))
            elif val_type == "float":
                value = float(value.translate(str.maketrans("\n{}","   ")))
            elif val_type == "str":
                value = value.translate(str.maketrans("\n{}","   ")).strip().lower()

            header_dict[key.strip()] = value
        line = header_file.readline()

    # Fill unused fields with None
    for key in field_dict:
        if key not in header_dict.keys():
            header_dict[key] = None

    header_file.close()
    return header_dict

def make_bad_band_list(wave_list,bad_band):
    
    good_list = np.array([1]*len(wave_list),dtype=np.uint8)
    wave_list_array = np.array(wave_list)
    
    for i in range(len(bad_band)):
        #print(bad_band[i][0], bad_band[i][1])
        wave_ind = (wave_list_array>=bad_band[i][0]) & (wave_list_array<=bad_band[i][1])
        good_list[wave_ind]=0
      
    return good_list.tolist()

def open_envi_4_nodata(file_name, header_dict,use_band=10):
    # use_band : starts from 0

    lines =  header_dict["lines"]
    columns =  header_dict["samples"]
    bands =   header_dict["bands"]
    interleave =  header_dict["interleave"]
    dtype = dtype_dict[header_dict["data type"]]
    no_data = header_dict['data ignore value']
    byte_order = header_dict['byte order']

    if byte_order == 1:
        img_endianness = 'big'
    else:
        img_endianness = 'little'
        
    offset = header_dict["header offset"]
    if offset is None:
        offset = 0

    if interleave == 'bip':
        img_shape = (lines, columns, bands)
    elif interleave == 'bil':
        img_shape = (lines, bands, columns)
    elif interleave == 'bsq':
        img_shape = (bands, lines, columns)
    else:
        print("ERROR: Unrecognized interleave type.")
    
    # If no_data value is not specified guess using image corners.
    if no_data is None:
        
        img_data = np.memmap(file_name,dtype = dtype, mode='r',
                                  shape = img_shape, offset=offset )
                                  
        if interleave == 'bip':
            band_pick = img_data[:,:,use_band]
        elif interleave == 'bil':
            band_pick = img_data[:,use_band,:]
        elif interleave == 'bsq':
            band_pick = img_data[use_band,:,:]                          
        
        up_l = band_pick[0,0]
        up_r = band_pick[0,-1]
        low_l = band_pick[-1,0]
        low_r = band_pick[-1,-1]
        
        if img_endianness != sys.byteorder:
            up_l = up_l.byteswap()
            up_r = up_r.byteswap()
            low_l = low_l.byteswap()
            low_r = low_r.byteswap()

        counts = {v: k for k, v in Counter([up_l,up_r,low_l,low_r]).items()}
        no_data = counts[max(counts.keys())]
        band_upper = np.percentile(band_pick[band_pick>0.5*no_data],98)
        band_lower = np.percentile(band_pick[band_pick>0.5*no_data],2)
        del img_data
        img_data=None

    
    return no_data, band_upper, band_lower


def main(argv):
    
    img_file = argv[0]
    out_dir = argv[1]
    
    bname = os.path.basename(img_file)
    new_json = "{}/{}_corr_ceoffs.json".format(out_dir,bname)
    
    logging.info(img_file)
    logging.info(new_json)
    
    iband = 30
    
    ref_header_dict = parse_envi_header(img_file+'.hdr')

    out_header_dict = {}
    #out_header_dict["fwhm"] =  list(ref_header_dict["fwhm"])
    out_header_dict["wavelength"] = ref_header_dict["wavelength"].tolist()
    
    out_header_dict["good_band_list"] = make_bad_band_list(out_header_dict["wavelength"],bad_band)
    
    #logging.info(dtype_dict[ref_header_dict["data type"]])
    #print(ref_header_dict["data ignore value"])
    no_data_val, band_98_val,band_02_val = open_envi_4_nodata(img_file, ref_header_dict,use_band = iband)
    #print(no_data_val,no_data_val.dtype,no_data_val.shape)
    logging.info("Band {:03d}: no data value {:.4f} 2 percentile {:.4f} 98 percentile {:.4f}".format(iband+1, no_data_val, band_02_val, band_98_val))
    
    if ref_header_dict["correction factors"] is not None:
        corr_list = ref_header_dict["correction factors"]
        no_data_val /= corr_list[iband]
    elif ref_header_dict["correction_factors"] is not None:
        corr_list = ref_header_dict["correction_factors"]
        no_data_val /= corr_list[iband]
    else:
        corr_list = [1.0]*ref_header_dict["bands"]

        
    corr_coeffs = [ 1/x if x !=0 else 1.0 for x in corr_list]
    
    if ref_header_dict["data type"] in [2,3,12,13,14,15]:
        logging.info("Integer Image, range 0-10000")
        corr_coeffs = [ x/10000 for x in corr_coeffs]
        refl_upper_lim = 10000
        
    elif ref_header_dict["data type"] in [4,5]:        
        if band_98_val> 20:
            corr_coeffs = [ x/10000 for x in corr_coeffs]
            refl_upper_lim = 10000
        else:
            refl_upper_lim = 1.0
            
        logging.info("Floating Point Image, range 0-{:.1f}".format(refl_upper_lim))    

    out_header_dict["reflectance_upper_limit"] = refl_upper_lim
    out_header_dict["file_name"] = img_file
    out_header_dict["no_data_value"] = no_data_val.item()
    out_header_dict["corr_coeffs"] = corr_coeffs
    out_header_dict["total_bands"] = ref_header_dict["bands"]
    out_header_dict["interleave"] = ref_header_dict["interleave"]
    
    with open(new_json, "w") as outfile:
        json.dump(out_header_dict, outfile, indent = 4)
    
if __name__ == "__main__":
    main(sys.argv[1:])