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


# python f:/tmp/anf_test/northup_2_rot_v5.py -img c:/tmp/topo_f130411_r12_ans.tif   -ref I:/hyspiri/sb_r12_ans  -od c:/tmp/

import argparse, os, sys, subprocess
import numpy as np
import tarfile
#from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator #, griddata, 

import logging
logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.DEBUG)

try:    
    from osgeo import gdal
except:
    import gdal

dict_enum_interp = {
    1: 'nearest',
    2: 'linear',
    3: 'slinear',
    4: 'cubic',
    5: 'quintic'
}    
    
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
    
def main():
    
    parser = argparse.ArgumentParser(description = "Rotation tool.")
    parser.add_argument("-img", help="Input image pathname",required=True, type = str)
    parser.add_argument("-od", help="Output directory for all resulting products", required=True, type = str)
    parser.add_argument("-ref", help="Reference rotated image pathname",required=True, type = str)
    parser.add_argument("-interp", help="Interplolation method code (1: 'nearest',2: 'linear',3: 'slinear',4: 'cubic',5: 'quintic')",default=2,required=False, type = int)
    #parser.add_argument("-epsg", help="EPSG", required=True, type = str)

    args = parser.parse_args()
    
    infile = args.img
    refimg = args.ref
    outpath = args.od
    if args.interp in [1,2,3,4,5]:
        interp_method = dict_enum_interp[args.interp]  #'nearest' 'linear' 'cubic'
    else:
        print('default linear')
        interp_method = dict_enum_interp[2]
    
    logging.info('Rotating {} to the orientation of {}'.format(infile, refimg)) #print(infile, refimg)    
    
    if refimg.startswith('/vsitar'): 
        basename = os.path.basename(refimg)
        tar = tarfile.open((refimg.split('.tar.gz')[0]).split('/vsitar/')[1]+'.tar.gz', "r:gz")
        for tarinfo in tar:
            if tarinfo.name.endswith(basename+'.hdr'):
                print(tarinfo.name)
                tar.extract(tarinfo,path=os.path.dirname(outpath))
                #hdrinfo_0=EHH.ENVI_Header(os.path.dirname(outpath)+'/'+basename+'.hdr')
                os.remove(os.path.dirname(outpath)+'/'+basename+'.hdr')
                break
    else:    
        pass #print("regular")
        #hdrinfo_0=EHH.ENVI_Header(os.path.splitext(refimg)[0]+".hdr")
        
    if infile.startswith('/vsitar'): 
        basename = os.path.basename(infile)
        tar = tarfile.open((infile.split('.tar.gz')[0]).split('/vsitar/')[1]+'.tar.gz', "r:gz")
        for tarinfo in tar:
            if tarinfo.name.endswith(basename+'.hdr'):
                print(tarinfo.name)
                tar.extract(tarinfo,path=os.path.dirname(outpath))
                #hdrinfo=EHH.ENVI_Header(os.path.dirname(outpath)+'/'+basename+'.hdr')
                os.remove(os.path.dirname(outpath)+'/'+basename+'.hdr')
                break
    else:    
        pass #print("regular")
        #hdrinfo=EHH.ENVI_Header(os.path.splitext(infile)[0]+".hdr")     
    
   
    ds_0 = gdal.Open(refimg)
    ds_1 = gdal.Open(infile)
    
    target_geotrans = ds_0.GetGeoTransform()
    #target_geotrans = update_rot(rot_angle_0, target_geotrans)    

    source_geotrans = ds_1.GetGeoTransform()
    #source_geotrans = update_rot(rot_angle, source_geotrans)
    
    if target_geotrans[1]*target_geotrans[2]==0:
        # NO rotation is needed, use nearest 
        interp_method = dict_enum_interp[1]
    
    #print(source_geotrans)
    #print(target_geotrans)
    #quit()
    ds_nband=ds_0.RasterCount
    in_ds_nband=ds_1.RasterCount
    ds_rows=ds_0.RasterYSize
    ds_cols=ds_0.RasterXSize
    ds_rows_1=ds_1.RasterYSize
    ds_cols_1=ds_1.RasterXSize
    #print("Ref IMG shape:",ds_rows,ds_cols)
    #print("North up imput IMG shape:",ds_1.RasterYSize,ds_1.RasterXSize)
    #quit()
    rot_mat_0 = np.array([[target_geotrans[1],target_geotrans[2],target_geotrans[0]],[target_geotrans[4],target_geotrans[5],target_geotrans[3]],[0,0,1]])
    rot_mat_1 = np.array([[source_geotrans[1],source_geotrans[2],source_geotrans[0]],[source_geotrans[4],source_geotrans[5],source_geotrans[3]],[0,0,1]])
    #print(np.dot(np.linalg.inv(rot_mat_1),rot_mat_0))
    #rot_mat_1 = np.array([[source_geotrans[1],source_geotrans[2],source_geotrans[0]+0.5*source_geotrans[1]],[source_geotrans[4],source_geotrans[5],source_geotrans[3]+0.5*source_geotrans[5]],[0,0,1]])
    
    #print(np.dot(np.linalg.inv(rot_mat_1),rot_mat_0))
    #quit()
    transfer_rot_mat = np.dot(np.linalg.inv(rot_mat_1),rot_mat_0)
    #transfer_rot_mat_1 = np.dot(np.linalg.inv(rot_mat_0),rot_mat_1)
    
    
    row_list = np.arange(ds_rows)
    col_list = np.arange(ds_cols)
    col_v, row_v = np.meshgrid(col_list, row_list)

    row_list_1 = np.arange(ds_rows_1)+0.5
    col_list_1 = np.arange(ds_cols_1)+0.5
    #col_v_1, row_v_1 = np.meshgrid(col_list_1, row_list_1,indexing='ij')
    
    # v5 v7 ok, v1 v2 v3 v4 v6 v8 not ok
    ind_grid_0 = np.dstack((col_v+0.5,row_v+0.5,np.ones(col_v.shape)))
    #ind_grid_1 = np.dstack((col_v_1+0.5,row_v_1+0.5,np.ones(col_v_1.shape))) # v1 30m v2 15m v5 15m v7 30m
    #ind_grid = np.dstack((col_v,row_v,np.ones(col_v.shape))) # v3 30m v4 15m v6 15m v8 30m
    
    '''
    # Start nearest neighbor searching
    # col_v+0.5,row_v+0.5 15m same as v2
    #
    ind_grid_v2 = np.dstack((col_v+0.5,row_v+0.5,np.ones(col_v.shape)))
    geo_grid_flt_v2 = np.dot(ind_grid_v2,rot_mat_0.T)[:,:,:2]
    dsm_row_list = np.arange(ds_1.RasterYSize)
    dsm_col_list = np.arange(ds_1.RasterXSize)
    dsm_col_v, dsm_row_v = np.meshgrid(dsm_col_list, dsm_row_list)
    grid_dsm = np.dstack((dsm_col_v,dsm_row_v,np.ones(dsm_col_v.shape)))
    geo_grid_dsm = np.dot(grid_dsm,rot_mat_1.T)[:,:,:2]
    src_points =np.concatenate([np.expand_dims(geo_grid_dsm[:,:,0].ravel(),axis=1),
                                np.expand_dims(geo_grid_dsm[:,:,1].ravel(),axis=1)],axis=1)
    #print(src_points.shape)
    tree = cKDTree(src_points,balanced_tree= False)
    dst_points = np.concatenate([geo_grid_flt_v2[:,:,0].ravel()[:,np.newaxis],
                                 geo_grid_flt_v2[:,:,1].ravel()[:,np.newaxis]],
                                 axis=1)
    indexes = tree.query(dst_points,k=1)[1]

    indices_int = np.unravel_index(indexes,(ds_1.RasterYSize,
                                            ds_1.RasterXSize))
    #print(len(indices_int))
    #for b,a in zip(indices_int[0][10:30],indices_int[1][10:30]):
    #    print(a,b)
    # End nearest neighbor searching
    '''
    #quit()
    #new_ind_grid = np.rint((np.dot(ind_grid,transfer_rot_mat.T)[:,:,:2])).astype(np.int16) # v1 v2 v3 v4
    #new_ind_grid = (np.dot(ind_grid,transfer_rot_mat.T)[:,:,:2]).astype(np.int16)  # v5 v6 v7 v8
    
    #new_ind_grid = np.dot(ind_grid_1,transfer_rot_mat_1.T)[:,:,:2]
    #new_ind_grid = new_ind_grid.reshape(-1, new_ind_grid.shape[-1])  # interpolation
    
    new_ind_grid_0 = (np.dot(ind_grid_0,transfer_rot_mat.T)[:,:,:2]) # interpolation
    
    
    #new_ind_grid_flt = np.dot(ind_grid,transfer_rot_mat.T)[:,:,:2]
    #print(new_ind_grid[:,:,1].max(),new_ind_grid[:,:,0].max())
    '''
    zip_ind = zip(new_ind_grid[:,:,1].ravel()[10:], new_ind_grid[:,:,0].ravel()[10:],new_ind_grid_flt[:,:,1].ravel()[10:], new_ind_grid_flt[:,:,0].ravel()[10:],indices_int[0][10:],indices_int[1][10:])
    for order,val in enumerate(zip_ind):
        a,b,c,d,e,f=val
        xxx,yyy=np.unravel_index(order,(ds_0.RasterYSize,
                                ds_0.RasterXSize))
        if (a>379*2) and (a<384*2) and (b>46*2) and (b<51*2):
            print(xxx,yyy,a,b,'|',int(np.rint(c)),int(np.rint(d)),'|',e,f,'|',"{:.2f}".format(c),"{:.2f}".format(d))
    '''
    
    #ind_new_1 = new_ind_grid[:,:,1].ravel()        
    #ind_new_0 = new_ind_grid[:,:,0].ravel()
    #ind_new_1 = np.maximum(ind_new_1,0)
    #ind_new_1 = np.minimum(ind_new_1,ds_1.RasterYSize-1)
    #ind_new_0 = np.maximum(ind_new_0,0)
    #ind_new_0 = np.minimum(ind_new_0,ds_1.RasterXSize-1)     
    
    #quit()
    basename = os.path.basename(infile.split('.tif')[0])
    driver = gdal.GetDriverByName('GTIFF')
    #ds_out = driver.Create(out_fn, out_col, out_row, 2, gdal.GDT_Float32, options=["INTERLEAVE=BIL"])
    ds_out = driver.Create(outpath+'/'+basename+'_rot.tif', ds_cols, ds_rows, in_ds_nband, gdal.GDT_Float32, options=["INTERLEAVE=BAND"])
    ds_out.SetGeoTransform(target_geotrans)
    ds_out.SetProjection(ds_0.GetProjection())

    
    for ii in range(in_ds_nband):
        #print(ii)
        src_band = ds_1.GetRasterBand(ii+1)
        band_name = src_band.GetDescription()
        data_array  =src_band.ReadAsArray().astype(np.float32)
        
        out_band_array  = np.ones((ds_rows, ds_cols))
       
        interp = RegularGridInterpolator((row_list_1,col_list_1), data_array,method=interp_method,bounds_error=False, fill_value=-9999)
        out_band_array = interp((new_ind_grid_0[:,:,1].ravel(),new_ind_grid_0[:,:,0].ravel()))
        #out_band_array = data_array[new_ind_grid[:,:,1].ravel(), new_ind_grid[:,:,0].ravel()] # not using cKDTree
        #out_band_array = data_array[indices_int] # using cKDTree
        Gdal_write_band(ds_out, ii+1, out_band_array.reshape((ds_rows, ds_cols)), band_name)
        #break
    
    ds_0 = None
    ds_1 = None
    ds_out=None 
    
    logging.info("Finish Rotation.")
    
    #subprocess.call("gdalwarp -tr 30 30 -tap -t_srs EPSG:32611 {0}/{1}_rot.tif {0}/{1}_rot_warptap.tif".format(outpath,basename),shell=True)

        
if __name__== "__main__":
    main()
    
        