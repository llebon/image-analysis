image-analysis
==============

MATLAB image analysis code used for availability assay experiments in 
"Fringe proteins modulate Notch-ligand cis and trans interactions to specify 
signaling states"

The MATLAB function segImagesDAPI.m takes three images: a segmentation 
image (in this case, an image with illumination in the 
low-wavelength DAPI range), and two corresponding fluorescence images.  

The segImagesDAPI.m function calls the file SegContour 3.m, the 
segmentation routine that identifies cell bodies in an image.  
SegContour3.m returns a mask where pixels in an identified cell are labeled and 
non-cell pixels have a value of 0.  

Next, segImagesDAPI.m calls ourclearborder.m, which eliminates any cells 
identified on the edge of the image.

Finally, the output mask of SegContour3.m is used to identify regions in 
the fluoresence images corresponding to individual cells.  The cell 
properties and fluorescence values for each cell in the image are 
calculated and recorded, and returned in the data structure segmeneted.
