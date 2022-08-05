
# Setup the environment for co-registration script

### 1. Arosics

There are two main ways of installing [Arosics](https://danschef.git-pages.gfz-potsdam.de/arosics/doc/installation.html). It requires a lot of dependecies. Therefore, it is not easy to install it in a normal way. I did not succeed in the pure **conda** way, but did it with a way mixing with **conda** and **pip** (install the 11 mentioned libraries with conda or conda-forge first, then use pip to install **Arosics**).

### 2. Google Earth Engine Python API

Install [Google Earth Engine Python API](https://developers.google.com/earth-engine/guides/python_install), and authenticate to the Earth Engine servers. It just needs to be done once. The credential file is saved at *$HOME/.config/earthengine/credentials*.

For **pip**, just follow the main page directly.

For **conda**, follow this [instruction](https://developers.google.com/earth-engine/guides/python_install-conda).


If GEE API is used in another machine that the default credentials storage location is not accessible, an alternative is to use service account (https://developers.google.com/earth-engine/guides/service_account). Private key got from the service account should be carried with the script. If the script needs to access the asset of your personel account, you need to grant the reading permission to the service account (something like *XXXXX@XXXX.iam.gserviceaccount.com*) in the GEE code editor. (**This method works after the new GEE credential method update. Tested on Aug 5, 2022.**)

```python

import ee
service_account = 'my-service-account@...gserviceaccount.com'
credentials = ee.ServiceAccountCredentials(service_account, 'privatekey.json')
ee.Initialize(credentials) 

#ee.Initialize() # use default path

```


### 3. Google Drive Python API

This [page](https://developers.google.com/drive/api/v3/quickstart/python) provides the step-by-step instruction for setting up a **Desktop Application**.

It takes some time to setup the part of ["Create a project and enable the API"](https://developers.google.com/workspace/guides/create-project) and the part of ["Create credentials"](https://developers.google.com/workspace/guides/create-credentials)(choose **Desktop Application**). After these steps, a "*credentials.json*" is retrieved and should be kept in the working directory.

The script for downloading file from Google Drive is based on the example script [*quickstart.py*](https://developers.google.com/drive/api/v3/quickstart/python#step_2_configure_the_sample) with some adjustments (change *SCOPES*, change the way of authorization from *run_local_server* to *run_console*, change the way of searching file).

After setting up, user needs to authenticate the script and get a token for each Google API. Authentication should be done once and the retrieved tokens can be used and renewed automatically for future executions. The token file for GEE API is saved in the library folder, while the token file (named '*token.json*' in the script) for Google Drive API is saved with the script.

Create a folder in Google Drive to save the output file temporarily (etc. "**tmp_gee_download**")

The list of my conda environment are attached here.("*command_to_env_list.txt*". The libraries installed by *pip* is not listed in the *conda env* command, so I attach another one for *pip* libraries)


# Script workflow

1. Provide a single band image (e.g. a band in the RED wavelength range). It can be ENVI, GeoTIFF, or any format GDAL can read. The following information will be extracted:
    * EPSG code of the image
    * Boundary coordinates of the image in LAT/LON
    * Acquisition date of the flightline
    * *test_gee2.py*

2. Use Google Earth Engine to generate a reference image, and save it to a folder in Google Drive.
   * Search for Sentinel-2 MSI Band 4 (10m) first, and Landsat-8 OLI Pan band (15m) if no matching image found in the Sentinel-2 archive. 
   * Search for the reference image with the AVIRIS boundary spatially
   * Other filtering criterior include date range (within N days before and after the acquisition date) and cloud cover percentage (<20%)
   * Export image to the specified Google Drive folder in the targeting EPSG projection.
   * The submitted GEE tasks can be viewed online at [GEE Task Manager](https://code.earthengine.google.com/tasks) (login required). Normal Task Status changes in this order: **READY-RUNNING-COMPLETED**.
   * Reference file name for **A.tif** will be **ref_A.tif**
   * *test_gee2.py* (*Two Google API* and related libraries are required)

3. Use Google Drive API to download the reference image to a local path. It (*gdrive_download_file.py*) is called inside *test_gee2.py*. Files will **NOT** be deleted by the script. They have to be removed manually.

4. Use Arosics to match the target (one-band AVIRIS image) and the reference image (one-band image).
   * Rotate the AVIRIS image to make it north up (in memory) (for image with non-zero rotation angle)
   * Generating Ground Control Points (GCPs) using the reference image and the rotated AVIRIS image (zero rotation angle)
   * Convert the image coordinates in rotated image (zero rotation angle) to the image coordinates in the original image (non-zero rotation angle) in the GCPs, and the geographic coordinates are kept the same.
   * PROJ ERROR can be ignored while using multiplecores
   * Save GCPs to JSON file
   * JSON file name for **A.tif** will be **A_GCPs.json**
   * *arosics_test.py* (*arosics* and related libraries are required)

5. Apply warping to any image product from the same flightline with the JSON GCPs using GDAL
   * Any product that are from the same flightline means the image size should be the same (columns and rows). Georeference information should be the same as well.
   * They can be reflectance image, predicted trait map, etc..
   * *apply_geocorr_json2.py* (*only GDAL* is required)


# Run the scripts in order
The first script pulls the reference image to a local folder, then the second script generates the GCPs. The final product is just a JSON file for one flightline.

```bash

$filename='f140603t01p00r08_rfl_v1b_img_topo_brdf_b29'
$indir='/input_dir/'
$outdir='/output_dir/'

python test_gee2.py --inputfiles "$indir"/"$filename".tif --outputdir $outdir

# the output reference image has a prefix of 'ref_', which can be change inside the script
python arosics_test.py -r "$outdir"/"ref"_"$filename".tif -t "$indir"/"$filename".tif -o $outdir  --GCP_Only

```

To correct any related image with the GCPs, run each of these.

```bash
# correct the original image
python apply_geocorr_json2.py  -t "$indir"/"$filename".tif -o $outdir -g "$outdir"/"$filename"_GCPs.json

# correct the first band of the original image with an output boundary, output in 15m pixel size
python apply_geocorr_json2.py  -t "$indir"/"$filename".tif -o $outdir -g "$outdir"/"$filename"_GCPs.json -b 1 --bound 294411 4087259 340026  4107815 -p 15

# or correct other product with the same georeference and the same image dimension
python apply_geocorr_json2.py  -t "$indir"/"$filename"_LMA -o $outdir -g "$outdir"/"$filename"_GCPs.json

```

