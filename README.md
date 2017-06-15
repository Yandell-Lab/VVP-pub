# VVP
Variant prioritization / burden test.  Version 1.5

## INSTALL
### DEPENDENCIES  

1. Gnu scientific library (https://www.gnu.org/software/gsl/)  
2. openmp compatible version of gcc.  If your compiler (clang) is not, you can remove the -fopenmp flag in the CMakelists.txt file.  Change the line that looks like: set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lz -lm -Wall -fopenmp -lgsl -lgslcblas") to set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lz -lm -Wall -lgsl -lgslcblas")
3.  make, cmake

### BUILD

In the VVP directory:

`mkdir build`

`cd build`

`cmake ../`

`make`

Make will build 3 executables:  build_background, VVP, VAAST3

Note:  This has been built and run on Mac laptops and Linux servers.  

## EXAMPLE RUNNING VVP

To see available parameters of the executables, run with the -h option.  

Before running VVP, a background must be built.  From the VVP directory:

`cd example`

`../build/build_background -i 1KG_cftr_background.recode.vep.vcf.gz -o 1KG.build -b 2500 -v CSQ,4,6,1,15`

The build_background step produces output to stdout for each of the variants in the background vcf file.  It also creates 5 different output files with extensions .bin, .bit, .max, .chr_offsets.txt, .dist.  These files contained information used by VVP.  

To run prioritize variants using VVP (in the example folder):

`../build/VVP -i target_spiked_simple.vcf.gz -d 1KG.build -v CSQ,4,6,1,15 1> target.spiked.vvp.out`

Then target_spiked.vvp.out contains the vvp output.

### PREPARE VCF FILE FOR ANALYSIS

The VVP pipeline does not support mulitallelic lines, these must first be decomposed.  We recommend using vt decompose to accomplish this task (http://genome.sph.umich.edu/wiki/Vt).

**Mandatory** preprocessing of a vcf file includes **multiallelic decomposition and VEP annotation**.  It is important to decompose **BEFORE** annnotating because of potential annotation collisions.  Our recommended steps are to use vt to decompose and normalize variants followed by VEP annotation.  No special options in VEP are required for the variant annotation.  Testing has been done with VEP v82.    







