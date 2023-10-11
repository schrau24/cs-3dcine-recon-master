# cs-3dcine-recon

#### MATLAB-based code for reconstructing 3d cine whole-heart data acquired using the AUMC PROUD patch

#### Description:
Compressed sensing (CS) based reconstruction of raw Philips data is performed using MATLAB-based implementations of 
ReconFrame and bart toolboxes. This code is specifically designed for reconstructing data acquired using the PROspective Undersampling in multiple Dimensions [(PROUD) patch](https://www.amc.nl/web/leren/research-62/research-amsterdam-umc/software.htm) built and validated at Amsterdam UMC, location AMC.

#### Required MATLAB toolboxes:
* *Image Processing Toolbox*
* *Parallel Computing Toolbox*

#### Additional toolboxes:
* [bart](https://github.com/mrirecon/bart), either on a linux system or using a [Windows-based install](https://bart-doc.readthedocs.io/en/latest/install.html). 
* [MRecon](https://www.gyrotools.com/gt/index.php/products/reconframe) from Gyrotools.

#### Getting started:
cs-3dcine-recon requires the definition of a reconstruction protocol file (.pro.dat). This file contains all information needed from the reconstruction code. 
To begin, open MATLAB at the *cs-3dcine-recon-master* folder and type:
```
run_defineprofile.m
```
