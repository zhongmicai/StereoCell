StereoCell allows software flexibility to support individual single-step module operations. 
The corresponding execution is placed in directory __scripts__.

## Stitching
```text
python .\stitching.py \
--tiles_path /data/SS200000135TL_D1
--output_file /result/stitched_image.tif
```
*  ```--tiles_path```  Image data produced by microscopy. 
* ```--output_file``` The output path of the stitched image file. 

## Registration
First, we use the dark line detection algorithm to get alpha (the chip is often placed with an angle (alpha) when experimenter 
using microscope taking images). By comparing the interval between the line on the chip and the corresponding interval on the 
image, the parameter s will be determined. Second, the algorithm will find the center of gravity (i0 and g0) of the image 
and the expression matrix, align the center of gravity, move the image in the area near the center of gravity and flip 
(flip horizontally and not flip sequentially) and rotate (0°, 90°, 180°, and 270° in clockwise direction). The algorithm 
will use hash algorithm to calculate the Hamming distance. The smallest distance corresponding parameters will be used as 
the final registration parameters. We use the dark line which is in specific arrangement on the chip as moving unit to ensure 
the accuracy of the registration. 


```text
python .\registration.py \
--image_file /result/stitched_image.tif \
--output_file /data/register_image.tif \
--gene_exp_data /data/SS200000135TL_D1.gem.gz \
--chip_no SS200000135TL_D1
```
*  ```--image_file```  Image data produced by microscopy. 
* ```--gene_exp_data``` Gene matrix of [Stereo-seq](https://bgi-australia.com.au/stomics). 
* ```--chip_no``` Serial number of STOmics chip. 
* ```--output_file``` The output path of single cell gene matrix produced by StereoCell. 

## Tissue segmentation
```text
python .\segmentation.py \
--type tissue
--image_file /result/register_image.tif \
--output_file /result/tissue_mask.tif \
```
*  ```--image_file```  Stain image, usually a registered image. 
* ```--type``` Is cell or tissue. 
* ```--output_file``` The output path of binarized mask file. 

## Cell segmentation
```text
python .\segmentation.py \
--type nuclei
--image_file /result/register_image.tif \
--output_file /result/tissue_mask.tif \
```
* ```--image_file```  Stain image, usually a registered image. 
* ```--type``` Is cell or tissue. 
* ```--output_file``` The output path of binarized mask file. 

## Cell labeling
```text
python labeling.py \
--mask_path D:\StereoCell\data\cell_mask.tif \
--gem_path D:\StereoCell\data\gene.gem.gz \
--output_path D:\StereoCell\data
```
* ```--image_file``` The path of binarized cell mask file. 
* ```--gene_exp_data``` Gene matrix of [Stereo-seq](https://bgi-australia.com.au/stomics). 
* ```--output_path``` The output path of single cell matrix. 

