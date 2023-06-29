## output of pipeline
* ```stitched_image.tif``` Whole slide image (Moving image)
* ```register_image.tif```  Registration of moving image. 
* ```nuclei_mask.tif``` Cell segmentation of registration.tif
* ```tissue_mask.tif```  Tissue segmentation of registration.tif
* ```nuclei_mask_profile.txt``` Single cell gene matrix with nuclei mask
* ```cell_mask_profile.txt``` Single cell gene matrix with cell labeling

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