# Fault-Detection-Using-Color
Code for paper - Z. Wang, D. Temel and G. AlRegib, "Fault detection using color blending and color transformations," 2014 IEEE Global Conference on Signal and Information Processing (GlobalSIP), Atlanta, GA, 2014, pp. 999-1003.

<p align="center">
  <img src=/Images/graphical.png/>
</p> 

### Paper
Preprint: https://ghassanalregibdotcom.files.wordpress.com/2016/10/wang2014_globalsip.pdf

IEEE: https://ieeexplore.ieee.org/document/7032271

This is a brief explanation and demonstration of detecting faults using color blending and color transformations.

### Citation
If you find our paper and repository useful, please consider citing our paper:  
```
@INPROCEEDINGS{7032271, 
author={Z. {Wang} and D. {Temel} and G. {AlRegib}}, 
booktitle={2014 IEEE Global Conference on Signal and Information Processing (GlobalSIP)}, 
title={Fault detection using color blending and color transformations}, 
year={2014}, 
pages={999-1003}, 
doi={10.1109/GlobalSIP.2014.7032271}, 
ISSN={}, 
month={Dec},}

```

### Abstract 
In the field of seismic interpretation, univariate data-based maps are commonly used by interpreters, especially for fault detection. In these maps, the contrast between target regions and the background is one of the main factors that affect the accuracy of interpretation. Since univariate data-based maps are not capable of providing a high-contrast representation, to overcome this issue, we turn them into multivariate data-based representations using color blending. We blend neighboring time sections or frames that are viewed in the time direction of migrated seismic volumes as if they corresponded to the red, green, and blue channels of a color image. Furthermore, to extract more reliable structural information, we apply color transformations. Experimental results show that the proposed method improves the accuracy of fault detection by limiting the average distance between detected fault lines and the ground truth into one pixel.

### Keywords
seismic interpretation, color space transformations, color blending, perception-based detection, skeletonization
### Code
* proposed_1576_Color.m” utilizes the proposed method to detect faults at time = 1576ms. Variable “bw3” in this function equals to “./cmpr/proposed_1576_Color.mat”. 
* proposed_1604_Color.m” utilizes the proposed method to detect faults at time = 1604ms. Variable “bw3” in this function equals to “./cmpr/proposed_1604_Color.mat”. This function generates all figures in the paper. 
* proposed_1624_Color.m: utilizes the proposed method to detect faults at time = 1576ms. Variable “bw3” in this function equals to “./cmpr/proposed_1624_Color.mat”.
* proposed_1576/1604/1624_noColor.m: are similar to the previous three m-files, in which the proposed method does not involve color representation.
* Zhang_1576/1604/1624.m: are similar to previous m-files, which utilize Zhang’s method [11] to label faults in time sections at 1576, 1604, and 1624ms. 
* cmpr.m in Folder “cmpr”: compares faults detected by various methods with manually labeled ones. Three methods are included in this function: proposed method with or without color representation and method proposed by Zhang [11]. This function generates all results in Table I. 



### Contact:

Ghassan AlRegib:  alregib@gatech.edu, https://ghassanalregib.com/, 

Dogancan Temel: dcantemel@gmail.com, http://cantemel.com/


