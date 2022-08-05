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

#python apply_geocorr_resize_rot_chunk.py -t /home/ye6/hys_test/gdrive_test/new_obs_ort/f130522t01p00r07_s2sub_rfl_v1a_img -o ./new_obs_ort/ -g ./new_obs_ort/f130522t01p00r07_s2sub_rfl_v1a_img_b29_GCPs.json -p 30 -r ./new_obs_ort/f130522t01p00r07_s2sub_rfl_v1a_img_b29.tif

# python apply_geocorr_resize_rot_chunk2.py -t /home/ye6/hys_test/gdrive_test/new_obs_ort/f130522t01p00r07_s2sub_rfl_v1a_img -o /home/ye6/hys_test/gdrive_test/temp/ -g ./new_obs_ort/f130522t01p00r07_s2sub_rfl_v1a_img_b29_GCPs.json -p 30 -r ./new_obs_ort/f130522t01p00r07_s2sub_rfl_v1a_img_b29.tif --corr /home/ye6/hys_test/gdrive_test/temp/f130522t01p00r07_s2sub_rfl_v1a_img_corr_ceoffs.json


#/home/ye6/hys_test/gdrive_test/results/f180622t01p00r06_corr_v1k1_img_topo_brdf_b29_GCPs.json

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


# gdal to numpy datatype conversion dictionary
dtype_dict = {
             gdal.GDT_Byte:np.uint8,
             gdal.GDT_Int16:np.int16,
             gdal.GDT_Int32:np.int32,
             gdal.GDT_Float32:np.float32,
             gdal.GDT_Float64:np.float64,
             gdal.GDT_UInt16:np.uint16,
             gdal.GDT_UInt32:np.uint32,
             #14:np.int64,
             #15:np.uint64
             }

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

def geocorr_warp_bands(in_ds,GCPs,band2process,outputBounds,in_nodata,out_nodata,proj,original_pixel_size, out_pixel_size):
    destname = "/vsimem/translate_step1.tif"
    
    translate_option = gdal.TranslateOptions(GCPs=GCPs, bandList=band2process)
    translate_ds = gdal.Translate(destname,in_ds,options=translate_option)
    warp_option_step1 = gdal.WarpOptions(
                                   outputBounds=outputBounds,
                                   format="ENVI",
                                   dstNodata=out_nodata, 
                                   dstSRS=proj,
                                   xRes=original_pixel_size,
                                   yRes=original_pixel_size,
                                   multithread=True,
                                   resampleAlg='near',
                                   targetAlignedPixels=False
                                  )

    northup_mem_img = '/vsimem/northup_img_step1.tif'
    northup_mem_img2 = '/vsimem/northup_img_step2.tif'
    warp_ds_step1 = gdal.Warp(northup_mem_img,translate_ds,options=warp_option_step1)

    warp_option_step2 = gdal.WarpOptions(
                                   outputBounds=outputBounds,
                                   format="ENVI",
                                   dstNodata=out_nodata, 
                                   dstSRS=proj,
                                   xRes=out_pixel_size,
                                   yRes=out_pixel_size,
                                   multithread=True,
                                   resampleAlg='average',
                                   targetAlignedPixels=True
                                  )    
    
    warp_ds = gdal.Warp(northup_mem_img2,warp_ds_step1,options=warp_option_step2)    
    
    for i_band in range(1, warp_ds.RasterCount + 1):
        current_band = warp_ds.GetRasterBand(i_band)
        current_band.SetNoDataValue(out_nodata)
    
    gdal.GetDriverByName('GTiff').Delete(destname)
    translate_ds=None
    warp_ds_step1=None
    
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

def generate_mask_ds(in_ds,out_nodata,band_order,reflectance_upper_limit=1.0):
    data_array = in_ds.GetRasterBand(min(band_order,in_ds.RasterCount)).ReadAsArray()
    mask = (data_array > 0.5*out_nodata) & (data_array < reflectance_upper_limit)
    
    driver = gdal.GetDriverByName('GTIFF')
    out_ds = driver.Create("/vsimem/temp_mask.tif", 
                           mask.shape[1], mask.shape[0], 1,
                           gdal.GDT_Float32, 
                           options=["INTERLEAVE=BAND"]
                          )
    out_ds.SetGeoTransform(in_ds.GetGeoTransform())
    out_ds.SetProjection(in_ds.GetProjection())
    
    out_ds.GetRasterBand(1).WriteArray(mask)
    out_ds.GetRasterBand(1).FlushCache()
    #print(out_ds)
    return out_ds

def generate_mask_ds_stack(in_ds,out_nodata,good_band_list,reflectance_upper_limit,interleave, img_name):
    
    good_list = np.where(np.array(good_band_list)==1)[0]
    #print(good_list.dtype)
    #quit()
    #print("Total good band: {}".format(good_list.shape[0]))

    if interleave=='bsq' or interleave=='BSQ':
        for ii in range(good_list.shape[0]):    
            #print("Band:",good_list[ii]+1)
            data_array = in_ds.GetRasterBand(int(good_list[ii])+1).ReadAsArray()
            mask_tmp = (data_array > 0.5*out_nodata) & (data_array < reflectance_upper_limit)

            if ii==0:
                mask = np.ones(data_array.shape) #.astype(np.bool)

            mask = mask * mask_tmp
    else: # 'bil', 'bip'
        mask = np.ones((in_ds.RasterYSize,in_ds.RasterXSize))
        if interleave=='bip' or interleave=='BIP':
            img_shape=(in_ds.RasterYSize,in_ds.RasterXSize,in_ds.RasterCount)
        if interleave=='bil' or interleave=='BIL':
            img_shape=(in_ds.RasterYSize,in_ds.RasterCount,in_ds.RasterXSize)
        print(in_ds.GetRasterBand(1).DataType)
        image_cube = np.memmap(img_name,dtype=dtype_dict[in_ds.GetRasterBand(1).DataType], mode='r',
                                  shape = img_shape, offset=0)
        
        chunk_size=400
        if interleave=='bip' or interleave=='BIP':            
            for chunk_order, line_start in enumerate(np.arange(0,in_ds.RasterYSize,chunk_size)): 
                line_chunk_list = list(range(line_start,min(line_start+1+chunk_size,in_ds.RasterYSize)))
                chunk_data = image_cube[line_chunk_list,:,:][:,:,good_list]
                mask_tmp = np.all((chunk_data > 0.5*out_nodata) & (chunk_data < reflectance_upper_limit),axis=2)
                mask[line_chunk_list,:]=mask_tmp
                
        if interleave=='bil' or interleave=='BIL':
            for chunk_order, line_start in enumerate(np.arange(0,in_ds.RasterYSize,chunk_size)): 
                line_chunk_list = list(range(line_start,min(line_start+1+chunk_size,in_ds.RasterYSize)))
                chunk_data = image_cube[line_chunk_list,:,:][:,good_list,:]
                mask_tmp = np.all((chunk_data > 0.5*out_nodata) & (chunk_data < reflectance_upper_limit),axis=1)
                mask[line_chunk_list,:]=mask_tmp
    
    driver = gdal.GetDriverByName('GTIFF')
    out_ds = driver.Create("/vsimem/temp_mask.tif", 
                           mask.shape[1], mask.shape[0], 1,
                           gdal.GDT_Float32, 
                           options=["INTERLEAVE=BAND"]
                          )
    out_ds.SetGeoTransform(in_ds.GetGeoTransform())
    out_ds.SetProjection(in_ds.GetProjection())
    
    out_ds.GetRasterBand(1).WriteArray(mask)
    out_ds.GetRasterBand(1).FlushCache()
    #print(out_ds)
    return out_ds
    

def save_image(out_img, warp_ds,out_ds_rows,out_ds_cols,out_bands,target_geotrans_new, out_proj, new_ind_grid,out_nodata,ds_out=None,start_band=0,corr_factors_list=None,out_mask=None):
      

    if ds_out is None:
        driver = gdal.GetDriverByName('ENVI')  #('GTIFF')
        if out_bands>20:
            ds_out = driver.Create(out_img, 
                                   out_ds_cols, out_ds_rows, out_bands, 
                                   gdal.GDT_Float32, 
                                   options=["INTERLEAVE=BIL"]
                                  )
        else:
            ds_out = driver.Create(out_img, 
                                   out_ds_cols, out_ds_rows, out_bands, 
                                   gdal.GDT_Float32, 
                                   options=["INTERLEAVE=BSQ"]
                                  )
        ds_out.SetGeoTransform(target_geotrans_new)
        ds_out.SetProjection(out_proj)
    
    #    band_offset = 0
    #else:
    #    print("Not None")
    #    band_offset = start_band # from 0

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
        data_array  =src_band.ReadAsArray().astype(np.float32)
        logging.info("Band {:3d} of {} {} ".format(ii+1+band_offset,ds_out.RasterCount,band_name))
        
        out_band_array  = np.ones((out_ds_rows, out_ds_cols))
        out_band_array = (data_array[ind_new_1, ind_new_0.ravel()]).reshape((out_ds_rows, out_ds_cols))
        
        if corr_factors_list is not None:
            #band_data = current_band.ReadAsArray()
            #print(band_data.shape,out_mask.shape)
            #temp_mask = band_data>0.5*in_nodata
            out_band_array[out_mask] = out_band_array[out_mask]*corr_factors_list[ii]
            out_band_array[~out_mask] = out_nodata
            #print(np.count_nonzero(out_mask),band_data.shape,out_nodata, in_nodata)
            #current_band.WriteArray(band_data)
            #current_band.FlushCache()
            #break        
        
        Gdal_write_band(ds_out, int(ii+1+band_offset), out_band_array, band_name)
    
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
    parser.add_argument('--ext', type=str, required=False, default="geocorr_norot", help="Output suffix (default '_geocorr_norot')")
    parser.add_argument('--corr', type=str, required=False, help="Correction factors of each band in json")

    args = parser.parse_args()

    tgt_img = args.t
    out_dir = args.o   
    gcp_file = args.g
    rotated_img = args.r
    out_nodata = args.nodata
    in_nodata = args.innodata
    output_suffix = '_'+args.ext
    
    file_basename = os.path.basename(tgt_img)
    print(file_basename)
    if file_basename.startswith('f'):
        data_year = file_basename[1:3] # aviris classic
    elif file_basename.startswith('ang'):
        data_year = file_basename[5:7]
    else:
        logging.info("File name is not recognized as an AVIRIS image file.")
        quit()
    
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
    
    original_pixel_size = np.sqrt(gt[1]**2+gt[2]**2)
    proj = osr.SpatialReference(wkt=ds.GetProjection())

    coreg_info = load_json_gcp_dict(gcp_file)
    GCPs = coreg_info['GCPList']
    
    if data_year=='18':
        mask_ds = generate_mask_ds_stack(ds,in_nodata,good_band_list,reflectance_upper_limit,interleave, tgt_img)        
    else:
        mask_ds = generate_mask_ds(ds,in_nodata,29,reflectance_upper_limit)
    
    warp_mask_ds = geocorr_warp_bands(mask_ds,GCPs,None,outputBounds,out_nodata,out_nodata,proj,original_pixel_size,out_pixel_size)    
   
    source_geotrans = warp_mask_ds.GetGeoTransform()
    
    out_ds_rows,out_ds_cols,target_geotrans_new,new_ind_grid = estimate_geometa(ds_0,target_geotrans,source_geotrans,out_pixel_size)
    
    out_mask_img = out_dir +'/mask_'+os.path.basename(tgt_img).split('.tif')[0]+output_suffix
    ds_out_mask=save_image(out_mask_img, warp_mask_ds,out_ds_rows,out_ds_cols,warp_mask_ds.RasterCount,target_geotrans_new, ds_0.GetProjection(),new_ind_grid,out_nodata)
    
    out_mask = (ds_out_mask.GetRasterBand(1).ReadAsArray()==1).astype(np.bool)
    #print(out_mask.shape)
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
        warp_ds = geocorr_warp_bands(ds,GCPs,band_chunk_list,outputBounds,in_nodata,out_nodata,proj,original_pixel_size,out_pixel_size)
        #print(warp_ds.RasterCount)
        #continue
        if chunk_order==0:
            ds_out=save_image(out_full_img, warp_ds,out_ds_rows,out_ds_cols,ds.RasterCount,target_geotrans_new, ds_0.GetProjection(),new_ind_grid,out_nodata,ds_out=None,start_band=0,corr_factors_list=corr_factors_list,out_mask=out_mask)
        else:
            ds_out=save_image(out_full_img, warp_ds,out_ds_rows,out_ds_cols,ds.RasterCount,target_geotrans_new, ds_0.GetProjection(),new_ind_grid,out_nodata,ds_out=ds_out,start_band=band_start,corr_factors_list=corr_factors_list,out_mask=out_mask)
            
        warp_ds = None    
    
    # Rotate the geocorrected image back to an angle same as the reference rotated image

    ds_out=None   
    ds = None    
    ds_0 = None
    ds_out_mask = None
    
    #gdal.GetDriverByName('GTiff').Delete(northup_mem_img)
    logging.info("Finish resampling.")
   

if __name__ == '__main__':
    main()
 