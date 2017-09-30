# Probabilistic Volume Rendering - Feature extraction using 3D Convolutional Sparse Coding
This is the repo for feature extraction in Probabilistic Volume Rendering

Data in *.tif format have been store in folders:
- data_bonsai# CT bonsai	
- data_kiwi	 # MRI kiwi
- data_ldCT	 # Low dose CT
- data_tooth # CT Tooth	

Run main_generate_all.m to generate all sparse maps, response maps (3 levels with hierarchical summation at each level)
Maps are stored in maps_* 

Base implementation was extended from 2D version of MATLAB code from http://brendt.wohlberg.net/software/SPORCO/

If you use this code please refer to these papers: 

```
@article{quan_intelligent_2017,
	title = {An {Intelligent} {System} {Approach} for {Probabilistic} {Volume} {Rendering} using {Hierarchical} 3D {Convolutional} {Sparse} {Coding}},
	volume = {PP},
	issn = {1077-2626},
	doi = {10.1109/TVCG.2017.2744078},
	journal = {IEEE Transactions on Visualization and Computer Graphics},
	author = {Quan, T. M. and Choi, J. and Jeong, H. and Jeong, W. K.},
	year = {2017},
}
```
```
@inproceedings{wohlberg_efficient_2014,
	title = {Efficient convolutional sparse coding},
	doi = {10.1109/ICASSP.2014.6854992},
	booktitle = {2014 {IEEE} {International} {Conference} on {Acoustics}, {Speech} and {Signal} {Processing} ({ICASSP})},
	author = {Wohlberg, B.},
	month = may,
	year = {2014},
	pages = {7173--7177},
}
```
```
@article{wohlberg_efficient_2015, 
	author={B. Wohlberg}, 
	journal={IEEE Transactions on Image Processing}, 
	title={Efficient Algorithms for Convolutional Sparse Representations}, 
	year={2016}, 
	volume={25}, 
	number={1}, 
	pages={301-315}, 
	doi={10.1109/TIP.2015.2495260}, 
	ISSN={1057-7149}, 
	month={Jan},
}
```
