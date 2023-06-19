StereoCell allows software flexibility to support individual single-step module operations. 
The corresponding execution is placed in directory __scripts__.

## Stitching
```text
python .\stitching.py \
--input D:\data\test\SS200000135TL_D1\SS200000135TL_D1
--output D:\data\test\stitched.tif
```
*  ```--input```  Image data produced by microscopy. 
* ```--output``` The output path of the stitched image file. 

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
--input D:\data\test\SS200000135TL_D1 \
--output D:\data\test\paper \
--matrix D:\data\test\SS200000135TL_D1.gem.gz \
--chipno SS200000135TL_D1
```
*  ```--input```  Image data produced by microscopy. 
* ```--matrix``` Gene matrix of [Stereo-seq](https://bgi-australia.com.au/stomics). 
* ```--chipno``` Serial number of STOmics chip. 
* ```--output``` The output path of single cell gene matrix produced by StereoCell. 

## Tissue segmentation
```text
python .\segmentation.py \
--type tissue
--input D:\StereoCell\data\image_6467_16800_512_512.tif \
--output D:\StereoCell\data\image_tissue.tif \
```
*  ```--input```  Stain image, usually a registered image. 
* ```--type``` Is cell or tissue. 
* ```--output``` The output path of binarized mask file. 

## Cell segmentation
```text
python .\segmentation.py \
--type cell
--input D:\StereoCell\data\image_6467_16800_512_512.tif \
--output D:\StereoCell\data\image_cell.tif \
```
* ```--input```  Stain image, usually a registered image. 
* ```--type``` Is cell or tissue. 
* ```--output``` The output path of binarized mask file. 

## Cell labeling
```text
python labeling.py \
--mask_path D:\StereoCell\data\cell_mask.tif \
--gem_path D:\StereoCell\data\gene.gem.gz \
--out_path D:\StereoCell\data
```
* ```--mask_path``` The path of binarized cell mask file. 
* ```--gem_path``` Gene matrix of [Stereo-seq](https://bgi-australia.com.au/stomics). 
* ```--out_path``` The output path of single cell matrix. 

