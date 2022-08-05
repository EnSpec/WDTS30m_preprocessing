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

# env: new_geocorr
# Usage:
# python aws_down_topo.py -h
# Example:
# python aws_down_topo.py --inputfiles /home/ye6/hys_test/gdrive_test/new_obs_ort/f130522t01p00r07_s2sub_obs_v1a_img_geocorr_norot --outputdir '/home/ye6/hys_test/gdrive_test/aws_dsm'

import os, sys, subprocess, argparse
import numpy as np

try:
    from osgeo import osr, ogr, gdal
except:
    import gdal, osr , ogr

import time
import logging

from get_aws_terrain import dem_generate 

def convert_coord(epsg_code,x,y):
    
    point_srs = osr.SpatialReference()
    point_srs.ImportFromEPSG ( int(epsg_code) ) 

    latlon_wgs84 = osr.SpatialReference()
    latlon_wgs84.ImportFromEPSG ( 4326 )

    #pntsrs2imgsrs = osr.CoordinateTransformation ( latlon_wgs84,point_srs)
    imgsrs2latlon = osr.CoordinateTransformation (  point_srs,latlon_wgs84)
    
    lat, lon, z_coord = imgsrs2latlon.TransformPoint(x,y)
    
    return lon, lat


def get_bound(infile,buffer_in_pixel=5):
    
    ds = gdal.Open(infile)

    geotrans = ds.GetGeoTransform()
    trans_mat = np.array([[geotrans[1],geotrans[2],geotrans[0]],[geotrans[4],geotrans[5],geotrans[3]]])

    spatial_res = np.sqrt(geotrans[1]**2+geotrans[2]**2)
    
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    proj.AutoIdentifyEPSG() #Key to get EPSG
    epsg_code = proj.GetAttrValue('AUTHORITY',1)

    corner_mat = np.array([
        [-buffer_in_pixel,-buffer_in_pixel,1],                       
        [-buffer_in_pixel,ds.RasterYSize+buffer_in_pixel,1],
        [ds.RasterXSize+buffer_in_pixel,ds.RasterYSize+buffer_in_pixel,1],
        [ds.RasterXSize+buffer_in_pixel,-buffer_in_pixel,1]
    ])

    geo_mat = trans_mat @ corner_mat.T

    latlon_mat = np.array([ convert_coord(int(epsg_code),x,y) for (x,y) in geo_mat.T])
    
    ds = None
    
    return [np.min(latlon_mat[:,0]),np.max(latlon_mat[:,0]),np.min(latlon_mat[:,1]),np.max(latlon_mat[:,1])], epsg_code,spatial_res
    #return [np.min(geo_mat[:,0]),np.max(geo_mat[:,0]),np.min(geo_mat[:,1]),np.max(geo_mat[:,1])], epsg_code,spatial_res


def main():

    parser = argparse.ArgumentParser()  
    parser.add_argument('--inputfiles', type=str, required=True, help="Input image file, format can be ENVI or GeoTIFF")
    parser.add_argument('--outputdir', type=str, required=True, help="Local output directory for the reference image downloaded from AWS")

    args = parser.parse_args()

    #n_input = len(args.inputfiles)
    out_dir = args.outputdir
    #print("{} input file(s) found.".format(n_input)) 

    latlon_mat, epsg_code,spatial_res = get_bound(args.inputfiles,buffer_in_pixel=10)
    #utm_bound, epsg_code,spatial_res = get_bound(args.inputfiles)
    
    flight_basename = os.path.basename(args.inputfiles).split('_')[0]
    #latlon_bound = [9999, -9999, 9999, -9999]  # lon_min, lon_max, lat_min, lat_max
    #latlon_bound = [min(latlon_bound[0],latlon_mat[0]),max(latlon_bound[1],latlon_mat[1]),min(latlon_bound[2],latlon_mat[2]),max(latlon_bound[3],latlon_mat[3])]

    out_prefix = 'topo_' #'sentinel2_'
    #print(spatial_res)

    aws_site="https://copernicus-dem-30m.s3.amazonaws.com/"
    #print(out_prefix+flight_basename, epsg_code, latlon_mat,spatial_res, aws_site,out_dir)
    #dem_generate(out_prefix+flight_basename, epsg_code, utm_bound,spatial_res, aws_site,out_dir)
    dem_generate(out_prefix+flight_basename, epsg_code, latlon_mat,spatial_res, aws_site,out_dir)
    
    
    #out_name = out_prefix+ img_info['imagename']

if __name__ == '__main__':
    main()
