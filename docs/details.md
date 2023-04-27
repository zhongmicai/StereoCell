## output of pipeline
* ```registration.tif```  Registration of moving image. 
* ```stitched.tif``` Whole slide image (Moving image)
* ```<chipno>_matrix.tif``` Fixed image generate by gene matrix 
* ```registration_cell_segmentation.tif``` Cell segmentation of registration.tif
* ```registration_tissue_segmentation.tif```  Tissue segmentation of registration.tif
* ```template.txt``` Intermediate files produced during the registration process
* ```<chipno>_***.ipr``` Image process record file
* ```bin_filtering.gem``` Gene matrix filter by tissue region

## How to improve results of cell segmentation

* By tissue segmentation
```shell
cd scripts

# improve results of cell segmentation
python .\mask_merge.py \
-s D:\data\weights\registration_cell_segmentation.tif \
-d D:\data\weights\registration_tissue_segmentation.tif \
-o registration_cell_segmentation_bk.tif
```