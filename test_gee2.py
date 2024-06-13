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

# Usage:
# python test_gee2.py -h
# Example:
# python test_gee2.py --inputfiles f130612t01p00r15_rfl_v1a_img_topo_brdf_b29.tif --outputdir '/output/' --gdrivedir "tmp_gee_folder"

import os, sys, subprocess, argparse
import ee
import numpy as np

try:    
    from osgeo import gdal, ogr, osr
except:
    import gdal, osr , ogr

import time

ee.Initialize()

SENTINEL2_match_band = ['B4']
LANDSAT8_match_band = ['B8']


def maskL8TOA(image):
    
  CirrusBitMask =  (1 << 12)  #(1 << 11) +
  cloudShadowBitMask =  (1 << 8) #(1 << 7) +
  cloudsBitMask = 1 << 4

  # Get the pixel QA band.
  qa = image.select('BQA')  #('pixel_qa')

  # Both flags should be set to zero, indicating clear conditions.
  mask = qa.bitwiseAnd(cloudsBitMask).eq(0).And(qa.bitwiseAnd(cloudShadowBitMask).lt(1))\
                                                .And(qa.bitwiseAnd(CirrusBitMask).lt(1))
  return image.updateMask(mask)


def maskL8sr(image):
    
  # Bits 3 and 5 are cloud shadow and cloud, respectively.
  CirrusBitMask = 1 << 2
  cloudShadowBitMask = 1 << 3
  cloudsBitMask = 1 << 4
    
  # Get the pixel QA band.
  qa = image.select('QA_PIXEL')  #('pixel_qa')
    
  # Both flags should be set to zero, indicating clear conditions.
  mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(qa.bitwiseAnd(cloudsBitMask).eq(0))\
                                                .And(qa.bitwiseAnd(CirrusBitMask).eq(0))
  return image.updateMask(mask)


def maskS2clouds(image):

  qa = image.select('QA60')

  # Bits 10 and 11 are clouds and cirrus, respectively.
  cloudBitMask = 1 << 10
  cirrusBitMask = 1 << 11

  # Both flags should be set to zero, indicating clear conditions.
  mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))

  return image.updateMask(mask) #.divide(10000)

def addNDVI(image):
  ndvi = image.normalizedDifference(['B8', 'B4']).multiply(10000).toInt16().rename('NDVI')
  return image.addBands(ndvi)

def addNDVI_L8(image):
  ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).multiply(10000).toInt16().rename('NDVI')
  return image.addBands(ndvi)


def convert_coord(epsg_code,x,y):
    
    # in the case of from 4326 to 4326
    if int(epsg_code)==4326:
      return x, y
    
    point_srs = osr.SpatialReference()
    point_srs.ImportFromEPSG ( int(epsg_code) ) 

    latlon_wgs84 = osr.SpatialReference()
    latlon_wgs84.ImportFromEPSG ( 4326 )

    #pntsrs2imgsrs = osr.CoordinateTransformation ( latlon_wgs84,point_srs)
    imgsrs2latlon = osr.CoordinateTransformation (  point_srs,latlon_wgs84)
    
    lat, lon, z_coord = imgsrs2latlon.TransformPoint(x,y)
    
    return lon, lat


def get_bound(infile):
    
    ds = gdal.Open(infile)

    geotrans = ds.GetGeoTransform()
    trans_mat = np.array([[geotrans[1],geotrans[2],geotrans[0]],[geotrans[4],geotrans[5],geotrans[3]]])

    proj = osr.SpatialReference(wkt=ds.GetProjection())
    proj.AutoIdentifyEPSG() #Key to get EPSG
    epsg_code = proj.GetAttrValue('AUTHORITY',1)

    corner_mat = np.array([
        [-2,-2,1],                       
        [-2,ds.RasterYSize+2,1],
        [ds.RasterXSize+2,ds.RasterYSize+2,1],
        [ds.RasterXSize+2,-2,1]
    ])

    geo_mat = trans_mat @ corner_mat.T

    latlon_mat = np.array([ convert_coord(int(epsg_code),x,y) for (x,y) in geo_mat.T])
    #print(latlon_mat.shape, np.min(latlon_mat[:,0]),np.max(latlon_mat[:,0]),np.min(latlon_mat[:,1]),np.max(latlon_mat[:,1]))
    
    ds = None
    
    return [np.min(latlon_mat[:,0]),np.max(latlon_mat[:,0]),np.min(latlon_mat[:,1]),np.max(latlon_mat[:,1])], epsg_code

# if input is a list of images
# boundary should contain all images, but the date and EPSG is from the first image
def get_img_list_geo(infile_list):
    
    epsg_code_list = []
    latlon_bound = [9999, -9999, 9999, -9999]  # lon_min, lon_max, lat_min, lat_max
    
    for infile in infile_list:
        latlon_mat, epsg_code = get_bound(infile)
        epsg_code_list=epsg_code_list+[epsg_code]
        latlon_bound = [min(latlon_bound[0],latlon_mat[0]),max(latlon_bound[1],latlon_mat[1]),min(latlon_bound[2],latlon_mat[2]),max(latlon_bound[3],latlon_mat[3])]

    polygonLinearRing = ee.Geometry.Polygon(
      [
        ee.Geometry.LinearRing(
          [
            [latlon_bound[0],latlon_bound[2]],
            [latlon_bound[1],latlon_bound[2]],
            [latlon_bound[1],latlon_bound[3]],
            [latlon_bound[0],latlon_bound[3]]
          ]
        )
      ],
      proj="EPSG:4326"
    )

    filebase = os.path.basename(infile_list[0])

    if filebase.startswith('ang'):
        date_str1 = filebase[3:7]+'-'+ filebase[7:9] + '-'+ filebase[9:11]
    elif filebase.startswith('f'):
        date_str1 = '20'+filebase[1:3]+'-'+ filebase[3:5] + '-'+ filebase[5:7] 

    return {'imagename':filebase.split('.tif')[0], 'boundary':polygonLinearRing,'epsg':epsg_code_list[0],'date_str':date_str1 }

## if input is a single images
def get_img_geo(infile):
    
    ds = gdal.Open(infile)

    geotrans = ds.GetGeoTransform()

    trans_mat = np.array([[geotrans[1],geotrans[2],geotrans[0]],[geotrans[4],geotrans[5],geotrans[3]]])

    proj = osr.SpatialReference(wkt=ds.GetProjection())
    proj.AutoIdentifyEPSG() #Key to get EPSG
    epsg_code = proj.GetAttrValue('AUTHORITY',1)

    corner_mat = np.array([
        [-2,-2,1],                       
        [-2,ds.RasterYSize+2,1],
        [ds.RasterXSize+2,ds.RasterYSize+2,1],
        [ds.RasterXSize+2,-2,1]
    ])

    geo_mat = trans_mat @ corner_mat.T

    latlon_mat = [ convert_coord(int(epsg_code),x,y) for (x,y) in geo_mat.T]

    polygonLinearRing = ee.Geometry.Polygon(
      [
        ee.Geometry.LinearRing(
          [
            [latlon_mat[0][0],latlon_mat[0][1]],
            [latlon_mat[1][0],latlon_mat[1][1]],
            [latlon_mat[2][0],latlon_mat[2][1]],
            [latlon_mat[3][0],latlon_mat[3][1]]
          ]
        )
      ],
      proj="EPSG:4326"
    )

    filebase = os.path.basename(infile)

    if filebase.startswith('ang'):
        date_str1 = filebase[3:7]+'-'+ filebase[7:9] + '-'+ filebase[9:11]
    elif filebase.startswith('f'):
        date_str1 = '20'+filebase[1:3]+'-'+ filebase[3:5] + '-'+ filebase[5:7] 

    return {'imagename':filebase.split('.tif')[0], 'boundary':polygonLinearRing,'epsg':epsg_code,'date_str':date_str1 }



def main():

    parser = argparse.ArgumentParser()  
    parser.add_argument('--inputfiles', type=str, required=True, nargs='+', help="Input image file or files, format can be ENVI or GeoTIFF")
    parser.add_argument('--outputdir', type=str, required=True, help="Local output directory for the reference image downloaded from Google Earth Engine via Google Drive")
    parser.add_argument('--gdrivedir', type=str, required=False, default="tmp_gee_download", help="Output temporarily Google Drive folder for the reference image from Google Earth Engine")
    parser.add_argument('--daybuffer', type=int, required=False, default=24, help="A buffer of day before and after the input date extracted from the input image")
    args = parser.parse_args()

    n_input = len(args.inputfiles)
    out_dir = args.outputdir
    out_gdrive_dir = args.gdrivedir
    print("{} input file(s) found.".format(n_input))     

    if n_input==1:
        img_info = get_img_geo(args.inputfiles[0])
    elif n_input>1:
        img_info = get_img_list_geo(args.inputfiles)    
    else:
        print("No image found.")
        quit()

    target_date = ee.Date(img_info['date_str'])

    lag = args.daybuffer
    image_col = ee.ImageCollection("COPERNICUS/S2_SR") \
                    .filterDate(target_date.advance(-lag, 'days'),target_date.advance(lag, 'days')) \
                    .filterBounds(img_info['boundary']) \
                    .map(maskS2clouds) \
                    .select(SENTINEL2_match_band)
                    #.map(addNDVI) \
                    #.select(['B2','B3','B4','NDVI'])
                    #.select(['B2','B3','B4']) 

    image_count = image_col.size().getInfo()
    print('Sentinel-2 Count: ', image_count)
    spatial_res = 10
    out_prefix = 'ref_' #'sentinel2_'
    band_subset = SENTINEL2_match_band  #['B2','B3','B4','NDVI']
    scale_factor = 1
    
    # old version "LANDSAT/LC08/C01/T1_SR"   4
    # "LANDSAT/LC08/C02/T1_L2" 4
    # LANDSAT/LC08/C01/T2_TOA 0 pan B8
    # LANDSAT/LC08/C01/T1_TOA 4 pan B8
    # LANDSAT/LC08/C01/T1_RT_TOA 4 pan B8
    # LANDSAT/LC08/C01/T2_SR 0 no pan
    # LANDSAT/LC08/C01/T1_SR 4 no pan

    if image_count==0:
        spatial_res = 15  #30
        out_prefix = 'ref_' #'l8_pan_'
        band_subset = LANDSAT8_match_band
        scale_factor = 10000
        image_col = ee.ImageCollection("LANDSAT/LC08/C01/T1_TOA") \
          .filterDate(target_date.advance(-lag, 'days'),target_date.advance(lag, 'days')) \
          .filterBounds(img_info['boundary']) \
          .filterMetadata('CLOUD_COVER', 'less_than', 80) \
          .map(maskL8TOA) \
          .select(band_subset)
          #.map(maskL8sr) \
          #.map(addNDVI_L8) \
          #.select(['SR_B2','SR_B3','SR_B4','NDVI'],['B2','B3','B4','NDVI'])

        image_count = image_col.size().getInfo()
        print('Landsat 8 Count: ', image_count)
        if image_count==0:
            sys.exit(1)

    out_name = out_prefix+ img_info['imagename']

    #greenest = image_col.qualityMosaic('NDVI')
    median = image_col.median().multiply(scale_factor).toInt16()

    task = ee.batch.Export.image.toDrive(image=ee.Image(median.select(band_subset)),  # an ee.Image object.
                                         region=img_info['boundary'],  # an ee.Geometry object.
                                         description=out_name,
                                         folder=out_gdrive_dir,
                                         fileNamePrefix=out_name,
                                         maxPixels=1e9,
                                         scale=spatial_res,
                                         crs='EPSG:{}'.format(img_info['epsg']))

    task.start()

    str_list = ['.   ','..  ','... ','....']
    ccc = 0
    
    # wait until the task is completed
    while task.active():
        ccc+=1
        time.sleep(10)
        print('Task Status: {} {}'.format(str_list[ccc % 4], task.status()['state']),end='\r')

    print('')

    subprocess.call('python gdrive_download_file.py {} {}'.format(out_name+'.tif',out_dir),shell=True)

if __name__ == '__main__':
    main()
