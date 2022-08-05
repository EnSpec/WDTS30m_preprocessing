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

#python apply_geocorr_resize_rot_chunk_for_obs.py -t ./new_obs_ort/f130522t01p00r07_s2sub_obs_v1a_img -o ./new_obs_ort/ -g ./new_obs_ort/f130522t01p00r07_s2sub_rfl_v1a_img_b29_GCPs.json -p 30 -r ./new_obs_ort/f130522t01p00r07_s2sub_rfl_v1a_img_b29.tif


import os, sys, argparse, json
import numpy as np

import logging
logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.DEBUG)

try:
  from osgeo import gdal, gdalconst, osr
except:
  import gdal, gdalconst, osr

import warnings
warnings.filterwarnings("ignore")

# correct the averaged sensor azimuth near north (360 deg), make range to be (-180,180)
def update_sensor_azimuth_before_warp(in_ds):
    sensor_az=in_ds.GetRasterBand(2).ReadAsArray()
    sensor_az[sensor_az>180]=sensor_az[sensor_az>180]-360
    in_ds.GetRasterBand(2).WriteArray(sensor_az)
    in_ds.GetRasterBand(2).FlushCache()

def update_sensor_azimuth_after_warp(in_ds):  
    sensor_az=in_ds.GetRasterBand(2).ReadAsArray()
    data_mask = (sensor_az > -9999/2 ) & (sensor_az<0)
    sensor_az[data_mask]=sensor_az[data_mask]+360
    in_ds.GetRasterBand(2).WriteArray(sensor_az)
    in_ds.GetRasterBand(2).FlushCache()
    
def load_json_gcp_dict(json_file):
    with open(json_file, "r") as outfile:
        coreg_info = json.load(outfile)
        
        for item, gcp_info in enumerate(coreg_info['GCPList']):
            x,y,z,pixel,line = gcp_info
            coreg_info['GCPList'][item] = gdal.GCP(x,y,z,pixel,line)
        
        logging.info("{} GCP points loaded.".format(len(coreg_info['GCPList'])))
        return coreg_info

def load_json_corr_dict(json_file):
    with open(json_file, "r") as outfile:
        corr_coeff_info = json.load(outfile)
        
        corr_coeffs = corr_coeff_info['corr_coeffs'] 
        no_data_val = corr_coeff_info['no_data_value']
        reflectance_upper_limit = corr_coeff_info['reflectance_upper_limit']
        good_band_list = corr_coeff_info["good_band_list"]
        interleave = corr_coeff_info["interleave"]
        
        print("{} band of correction coeffs are loaded.".format(len(corr_coeffs)))
        return corr_coeffs, no_data_val, reflectance_upper_limit,good_band_list,interleave    
    
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

def geocorr_warp_bands(in_ds,GCPs,band2process,outputBounds,out_nodata,proj,out_pixel_size):
    destname = "/vsimem/translate.tif"
    
    translate_option = gdal.TranslateOptions(GCPs=GCPs, bandList=band2process)
    translate_ds = gdal.Translate(destname,in_ds,options=translate_option)
    
    if translate_ds.RasterCount>1: #obs_ort
        update_sensor_azimuth_before_warp(translate_ds)
    
    warp_option = gdal.WarpOptions(outputBounds=outputBounds,
                                   format="ENVI",
                                   dstNodata=out_nodata, 
                                   dstSRS=proj,
                                   xRes=out_pixel_size,
                                   yRes=out_pixel_size,
                                   multithread=True,
                                   resampleAlg='average',
                                   targetAlignedPixels=True
                                  )

    northup_mem_img = '/vsimem/northup_img.tif'
    warp_ds = gdal.Warp(northup_mem_img,translate_ds,options=warp_option)

    if translate_ds.RasterCount>1: #obs_ort
        update_sensor_azimuth_after_warp(warp_ds)

    for i_band in range(1, warp_ds.RasterCount + 1):
        warp_ds.GetRasterBand(i_band).SetNoDataValue(out_nodata)
    
    gdal.GetDriverByName('GTiff').Delete(destname)
    translate_ds=None
    
    return warp_ds

def estimate_geometa(ds_ref_rot,target_geotrans,source_geotrans,out_pixel_size):
    
    # Add a buffer
    ds_rows=ds_ref_rot.RasterYSize
    ds_cols=ds_ref_rot.RasterXSize
    #print(ds_rows,ds_cols)

    rot_mat_0 = np.array([[target_geotrans[1],target_geotrans[2],target_geotrans[0]],[target_geotrans[4],target_geotrans[5],target_geotrans[3]],[0,0,1]])
    rot_mat_1 = np.array([[source_geotrans[1],source_geotrans[2],source_geotrans[0]],[source_geotrans[4],source_geotrans[5],source_geotrans[3]],[0,0,1]])
    
    in_pixel_size = np.sqrt(target_geotrans[1]**2+target_geotrans[2]**2)
    res_ratio = out_pixel_size / in_pixel_size
    #print("res_ratio",res_ratio,rot_mat_0,rot_mat_1)
    buff_pixel = 1
    new_ul_coord = rot_mat_0@(np.array([[-buff_pixel,-buff_pixel,1],[ds_cols+buff_pixel,-buff_pixel,1],[-buff_pixel,ds_rows+buff_pixel,1],[ds_cols+buff_pixel,ds_rows+buff_pixel,1]]).T) //out_pixel_size * out_pixel_size
    #print("new_ul_coord:",new_ul_coord)
    target_geotrans_new = [ new_ul_coord[0,0], target_geotrans[1]*res_ratio,target_geotrans[2]*res_ratio, new_ul_coord[1,0], target_geotrans[4]*res_ratio,target_geotrans[5]*res_ratio  ]
    #print(target_geotrans_new)

    rot_mat_0_new = np.array([[target_geotrans[1]*res_ratio,target_geotrans[2]*res_ratio,new_ul_coord[0,0]],[target_geotrans[4]*res_ratio,target_geotrans[5]*res_ratio,new_ul_coord[1,0]],[0,0,1]])    
    
    transfer_rot_mat = np.dot(np.linalg.inv(rot_mat_1),rot_mat_0_new)
    
    out_ds_rows = int((ds_rows+2*buff_pixel)/res_ratio)
    out_ds_cols = int((ds_cols+2*buff_pixel)/res_ratio)
    row_list = np.arange(out_ds_rows)
    col_list = np.arange(out_ds_cols)
    col_v, row_v = np.meshgrid(col_list, row_list)
    
    logging.info("Dimension of output image: {}, {}".format(out_ds_rows,out_ds_cols))
    #print(warp_ds.ReadAsArray().shape)
    #ind_grid = np.dstack((col_v-0.5,row_v+0.5,np.ones(col_v.shape)))
    ind_grid = np.dstack((col_v+0.5,row_v+0.5,np.ones(col_v.shape)))
    
    new_ind_grid = np.dot(ind_grid,transfer_rot_mat.T)[:,:,:2].astype(np.int16)
    #new_ind_grid = np.rint((np.dot(ind_grid,transfer_rot_mat.T)[:,:,:2])-0.5).astype(np.int16)
    
    #print(new_ind_grid[:,:,1].min(),new_ind_grid[:,:,1].max())
    #print(new_ind_grid[:,:,0].min(),new_ind_grid[:,:,0].max())
    
    return out_ds_rows,out_ds_cols,target_geotrans_new,new_ind_grid

def generate_mask_ds(in_ds,out_nodata,band_order):
    data_array = in_ds.GetRasterBand(min(band_order,in_ds.RasterCount)).ReadAsArray()
    mask = (data_array > 0.5*out_nodata)
    
    driver = gdal.GetDriverByName('GTIFF')
    out_ds = driver.Create("/vsimem/temp_mask.tif", mask.shape[1], mask.shape[0], 1, gdal.GDT_Float32, options=["INTERLEAVE=BAND"])
    out_ds.SetGeoTransform(in_ds.GetGeoTransform())
    out_ds.SetProjection(in_ds.GetProjection())
    
    out_ds.GetRasterBand(1).WriteArray(mask)
    out_ds.GetRasterBand(1).FlushCache()
    #print(out_ds)
    return out_ds
    

def save_image(out_img, warp_ds,out_ds_rows,out_ds_cols,out_bands,target_geotrans_new, out_proj, new_ind_grid,ds_out=None,start_band=0):
      

    if ds_out is None:
        driver = gdal.GetDriverByName('ENVI')  #('GTIFF')
        if out_bands>20:
            ds_out = driver.Create(out_img, out_ds_cols, out_ds_rows, out_bands, gdal.GDT_Float32, options=["INTERLEAVE=BIL"])
        else:
            ds_out = driver.Create(out_img, out_ds_cols, out_ds_rows, out_bands, gdal.GDT_Float32, options=["INTERLEAVE=BSQ"])
        ds_out.SetGeoTransform(target_geotrans_new)
        ds_out.SetProjection(out_proj)
        
    #    band_offset = 0
    #else:
    #    print("Not None")
    #    band_offset = start_band # from 0
    #if target_geotrans_new[1]*target_geotrans_new[2] == 0:
    #    rotated = False
    #else:
    #    rotated = True

    band_offset = start_band
            
    ind_new_1 = new_ind_grid[:,:,1].ravel()        
    ind_new_0 = new_ind_grid[:,:,0].ravel()
    ind_new_1 = np.maximum(ind_new_1,0)
    ind_new_1 = np.minimum(ind_new_1,warp_ds.RasterYSize-1)
    ind_new_0 = np.maximum(ind_new_0,0)
    ind_new_0 = np.minimum(ind_new_0,warp_ds.RasterXSize-1)  
    
    for ii in range(warp_ds.RasterCount):
        src_band = warp_ds.GetRasterBand(ii+1)
        band_name = src_band.GetDescription()
        data_array = src_band.ReadAsArray().astype(np.float32)
        logging.info("Band {:3d} {} of {}".format(ii+1+band_offset,band_name,ds_out.RasterCount))

        out_band_array  = np.ones((out_ds_rows, out_ds_cols))
        out_band_array = data_array[ind_new_1, ind_new_0]
            
        Gdal_write_band(ds_out, int(ii+1+band_offset), out_band_array.reshape((out_ds_rows, out_ds_cols)), band_name)
    
    return ds_out
    #ds_out=None


def main():

    parser = argparse.ArgumentParser()      
    parser.add_argument('-b', type=int, required=False,help="Band number in the image to warp (index starts at 1)",nargs='+')
    parser.add_argument('-t', type=str, required=True,help="Target image")
    parser.add_argument('-o', type=str, required=True,help="Output directory")
    parser.add_argument('-g', type=str, required=True,help="GCP JSON file")
    parser.add_argument("-r", type=str, required=True,help="Reference rotated image name")

    parser.add_argument('--nodata', type=float, default=-9999, required=False,help="Nodata value in ouput image")
    parser.add_argument('--innodata', type=float, default=-9999, required=False,help="Nodata value in input image")
    parser.add_argument('-p', type=float, default=30, required=False,help="Output pixel size (same unit of imput image)")
    parser.add_argument('--bound', type=float,nargs='+', required=False, help="Order: MinX, MinY, MaxX, MaxY")
    parser.add_argument('--ext', type=str, required=False, default="obs_geocorr_norot", help="Output suffix (default '_obs_geocorr_norot')")
    parser.add_argument('--mask', type=bool, required=False, default=False, help="Output contribution mask (default: False)")
    parser.add_argument('--corr', type=str, required=False, help="Correction factors of each band in json")
    
    args = parser.parse_args()

    tgt_img = args.t
    out_dir = args.o   
    gcp_file = args.g
    rotated_img = args.r
    out_nodata = args.nodata
    in_nodata = args.innodata
    output_suffix = '_'+args.ext
    bool_out_mask=args.mask
    #print(args.b)
    #outputBounds = None #[543900,3743360,  550000, 3748000]
    
    if args.b is not None:
        band2process = args.b  #int(args.b[0])+1
    else:
        band2process = None
        
    if args.bound is not None:
        if len(args.bound)>=4:
            outputBounds = args.bound
        else:
            outputBounds = None
    else:
        outputBounds = None
      
    #logging.info("Bound: {} {} {} {}".format(tuple(args.bound)))
    #quit()
    

    ds = gdal.Open(tgt_img) #, gdal.GA_Update
    ds_0 = gdal.Open(rotated_img)
    gt = ds.GetGeoTransform()
    target_geotrans = ds_0.GetGeoTransform()

    if args.p is not None:
        out_pixel_size = args.p
    else:
        out_pixel_size = np.sqrt(gt[1]**2+gt[2]**2)

    if args.corr is not None:
        corr_factors_list,in_nodata,reflectance_upper_limit, good_band_list,interleave = np.array(load_json_corr_dict(args.corr))
        logging.info("Correction coefficients loaded: no data value {}".format(in_nodata))
    else:
        corr_factors_list = np.array([1.0]*ds.RasterCount)   
        
        
    proj = osr.SpatialReference(wkt=ds.GetProjection())

    coreg_info = load_json_gcp_dict(gcp_file)
    GCPs = coreg_info['GCPList']
    
    mask_ds = generate_mask_ds(ds,in_nodata,29)
    warp_mask_ds = geocorr_warp_bands(mask_ds,GCPs,None,outputBounds,out_nodata,proj,out_pixel_size)
    
    
    source_geotrans = warp_mask_ds.GetGeoTransform()
    
    out_ds_rows,out_ds_cols,target_geotrans_new,new_ind_grid = estimate_geometa(ds_0,target_geotrans,source_geotrans,out_pixel_size)
    #print(target_geotrans_new, source_geotrans)
    
    #quit()
    
    out_mask_img = out_dir +'/mask_'+os.path.basename(tgt_img).split('.tif')[0]+output_suffix
    
    if bool_out_mask:
        ds_out=save_image(out_mask_img, warp_mask_ds,out_ds_rows,out_ds_cols,warp_mask_ds.RasterCount,target_geotrans_new, ds_0.GetProjection(),new_ind_grid)
    ds_out=None
    warp_mask_ds=None
    #quit()
    #out_full_img = out_dir +'/full_'+os.path.basename(tgt_img).split('.tif')[0]+output_suffix
    out_full_img = out_dir +'/'+os.path.basename(tgt_img).split('.tif')[0]+output_suffix
    
    chunk_size=16
    for chunk_order, band_start in enumerate(np.arange(0,ds.RasterCount,chunk_size)):
        #print(chunk_order,band_start)
        #continue
        
        band_chunk_list = list(range(band_start+1,min(band_start+1+chunk_size,ds.RasterCount+1)))
        #print(band_chunk_list)
        warp_ds = geocorr_warp_bands(ds,GCPs,band_chunk_list,outputBounds,out_nodata,proj,out_pixel_size)
        #print(warp_ds.RasterCount)
        #continue
        if chunk_order==0:
            ds_out=save_image(out_full_img, warp_ds,out_ds_rows,out_ds_cols,ds.RasterCount,target_geotrans_new, ds_0.GetProjection(),new_ind_grid,ds_out=None,start_band=0)
        else:
            ds_out=save_image(out_full_img, warp_ds,out_ds_rows,out_ds_cols,ds.RasterCount,target_geotrans_new, ds_0.GetProjection(),new_ind_grid,ds_out=ds_out,start_band=band_start)
            
        warp_ds = None    
    
    # Rotate the geocorrected image back to an angle same as the reference rotated image

    ds_out=None   
    ds = None    
    ds_0 = None
    
    #gdal.GetDriverByName('GTiff').Delete(northup_mem_img)
    logging.info("Finish resampling.")
   

if __name__ == '__main__':
    main()
 