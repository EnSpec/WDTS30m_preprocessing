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

import os, sys
import numpy as np
import json

import logging
logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.DEBUG)

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

def write_envi_header(output_name,header_dict):
    """Write ENVI header file to disk.
    Args:
        output_name (str): Header file pathname.
        header_dict (dict): Populated ENVI header dictionary..
    Returns:
        None.
    """

    header_file = open(output_name + ".hdr",'w+')
    header_file.write("ENVI\n")

    for key in header_dict.keys():
        value = header_dict[key]
        # Convert list to comma seperated strings
        if isinstance(value,(list,np.ndarray)):
            value = "{%s}" % ",".join(map(str, value))
        else:
            value = str(value)
        # Skip entires with nan as value
        if value != 'None':
            header_file.write("%s = %s\n" % (key,value))
    header_file.close()

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

def load_json_corr_dict(json_file):
    with open(json_file, "r") as outfile:
        corr_coeff_info = json.load(outfile)
        good_band_list = corr_coeff_info["good_band_list"]

        return good_band_list

def main(argv):
    
    ref_hdr_file = argv[0]
    target_hdr_file = argv[1]
    corr_json = argv[2]
    
    logging.info(ref_hdr_file)
    logging.info(target_hdr_file)
    
    ref_header_dict = parse_envi_header(ref_hdr_file)
    out_header_dict = parse_envi_header(target_hdr_file)
    out_header_dict["fwhm"] =  ref_header_dict["fwhm"]
    #out_header_dict["bbl"] =  ref_header_dict["bbl"]
    out_header_dict["bbl"] = load_json_corr_dict(corr_json)
    out_header_dict["wavelength"] = ref_header_dict["wavelength"]
    del out_header_dict["band names"]
    
    write_envi_header(target_hdr_file.split(".hdr")[0],out_header_dict)
                                        
if __name__ == "__main__":
    main(sys.argv[1:])
