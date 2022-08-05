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
https://copernicus-dem-30m.s3.amazonaws.com/
https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com/
"""

import os
import glob
#import tarfile
import logging
import numpy as np
import subprocess
import pandas as pd

from rtree import index

try:
    from osgeo import gdal
except:
    import gdal

import requests


def Gdal_write_band(dst_raster, band_number, band_array, band_name, no_data=-9999):
    '''
    Write band to a gdal raster.
    band_number: int
    band number in the dst raster
    band_array: np array
    the array needs to be written
    band_name: str
    name of the band
    no_data: int or float
    no data value, default to -9999
    '''
    outband = dst_raster.GetRasterBand(band_number)
    outband.SetNoDataValue(no_data)
    outband.SetDescription(band_name)
    outband.WriteArray(band_array)
    outband.FlushCache()

# https://github.com/EnSpec/sister/blob/cf28213698da057f72358555e50befb8e1a9c2ff/sister/utils/misc.py
def download_file(file,url):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(file, 'wb') as f:
            for chunk in r.iter_content(chunk_size=int(1E8)):
                f.write(chunk)

# modified from https://github.com/EnSpec/sister/blob/master/sister/utils/terrain.py
def dem_generate(flightline_name, flightline_epsg_code,latlon_bound,spatial_res, elev_dir,temp_dir):
    '''
    Args:
        flightline_name (str): File basename of the flightline
        flightline_epsg_code (int): Projection's EPSG code of the flightline
        latlon_bound (float list): Longitude Min, Longitude Max, Latitude Min, Latitude Max
        elev_dir (str): Directory of zipped elevation tiles
        temp_dir (str): Temporary output directory
    Returns:
        dem gdal dataset : Elevation dataset.
    '''
    
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.DEBUG)
    
    # Get extents of image
    lon_min = latlon_bound[0]
    lon_max = latlon_bound[1]
    lat_min = latlon_bound[2]
    lat_max = latlon_bound[3]
    #x_min = utm_bound[0]
    #x_max = utm_bound[1]
    #y_min = utm_bound[2]
    #y_max = utm_bound[3]    

    if 'aws' in elev_dir:
        tiles = pd.read_csv(elev_dir + 'tileList.txt',header = None).values.flatten()
    else:
        tiles = glob.glob(elev_dir + '/Copernicus_DSM_COG*.tif')

    idx = index.Index(properties=index.Property())

    #Get list of intersecting tiles
    for i, tile in enumerate(tiles):
        lat,ign,lon = os.path.basename(tile).replace('_COG','').split('_')[3:6]
        if 'W' in lon:
            lon = -1*float(lon[1:])
        else:
            lon = float(lon[1:])
        if 'S' in lat:
            lat = -1*float(lat[1:])
        else:
            lat = float(lat[1:])
        idx.insert(i,(lon,lat,lon+1,lat+1))
    tiles_intersect = [tiles[n] for n in idx.intersection((lon_min, lat_min, lon_max, lat_max))]

    tile_file_list = []
    #print(len(tiles_intersect))
    if len(tiles_intersect) == 0:
        constant_elev = float(input("No overlapping tiles found, enter constant elevation for scene (m): "))
        elevation = np.ones(longitude.shape) * constant_elev
    else:
        tile_string = "Found %s intersecting elevation tiles:" % len(tiles_intersect)
        for tile in tiles_intersect:
            tile_string+= '\n\t%s' % tile
            #print(tile_string)
            if 'aws' in elev_dir:
                tile_url = "%s%s/%s.tif" % (elev_dir,tile,tile)
                tile_file = "%s%s.tif" % (temp_dir,tile)
                if not os.path.exists(tile_file):
                    #print("download...")
                    logging.info("Downloading {}...".format(tile))
                    download_file(tile_file,tile_url)
                tile_file_list+=[tile_file]
            else:
                tile_file_list+=["%s%s.tif" % (elev_dir,tile)]
            #    with tarfile.open(tile, 'r') as tar_ref:
            #        tar_ref.extractall(temp_dir)
        logging.info(tile_string)
    #print(tile_file_list)    
    #quit()
    logging.info('Merging DEM tiles (VRT)')
    dem_file  = '%stemp_dem.vrt' % temp_dir
    #vrt_opt = gdal.BuildVRTOptions(outputBounds=[lon_min, lat_max, lon_max, lat_min])
    gdal.BuildVRT(dem_file, tile_file_list) #,options=vrt_opt)
    
    warp_opt = gdal.WarpOptions(dstSRS="EPSG:{}".format(flightline_epsg_code),outputBounds=[lon_min, lat_min, lon_max, lat_max],outputBoundsSRS="EPSG:4326", xRes=spatial_res,yRes=spatial_res,resampleAlg="cubicspline", targetAlignedPixels=True) #,outputBounds=[lon_min, lat_min, lon_max, lat_max] #projWinSRS="EPSG:4326",projWin=[lon_min,lat_max,lon_max,lat_min],
    dsm_ds_clip = gdal.Warp("/vsimem/dsm_subset.tif", dem_file, options=warp_opt)
    
    elevation = dsm_ds_clip.GetRasterBand(1).ReadAsArray()
    
    if np.sum(elevation<0) > 0:
        logging.warning('Elevations below sea level found, setting to 0m')
        elevation[elevation<0] =0
        dsm_ds_clip.GetRasterBand(1).WriteArray(elevation)
        dsm_ds_clip.GetRasterBand(1).FlushCache()
    
    #dem_opt = gdal.DEMProcessingOptions(azimuth=file_dict["Solar_Azimuth_Angle"],altitude=90-file_dict["Solar_Zenith_Angle"])

    logging.info('Calculating slope')
    dsm_slope_ds = gdal.DEMProcessing("/vsimem/dsm_processing_slope.tif", dsm_ds_clip, "slope") #, options=dem_opt)
    logging.info('Calculating aspect')
    dsm_aspect_ds = gdal.DEMProcessing("/vsimem/dsm_processing_aspect.tif", dsm_ds_clip, "aspect") #, options=dem_opt)

    slope=dsm_slope_ds.GetRasterBand(1).ReadAsArray()
    aspect=dsm_aspect_ds.GetRasterBand(1).ReadAsArray()
    
    
    # write output file
    out_fn = '{}/{}_tmp_terrain.tif'.format(temp_dir,flightline_name)
    driver = gdal.GetDriverByName('GTiff')
    ds_out = driver.Create(out_fn, elevation.shape[1], elevation.shape[0], 3, gdal.GDT_Float32)
    ds_out.SetGeoTransform(dsm_ds_clip.GetGeoTransform())
    ds_out.SetProjection(dsm_ds_clip.GetProjection())

    Gdal_write_band(ds_out, 1, slope, "Slope", no_data=-9999)
    Gdal_write_band(ds_out, 2, aspect, "Aspect", no_data=-9999)
    Gdal_write_band(ds_out, 3, elevation, "DSM", no_data=-9999)
    
    dsm_ds_clip=None
    dsm_aspect_ds=None
    dsm_slope_ds=None
    ds_out=None
    


