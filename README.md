# VVP
Variant prioritization / burden test.  Version 1.5

## INSTALL
### DEPENDENCIES  

1. Gnu scientific library (https://www.gnu.org/software/gsl/)  
2. openmp compatible version of gcc.  If your compiler (clang) is not, you can remove the -fopenmp flag in the Makefile.  Change the line that looks like: CFLAGS = -lz -lm -O3 -lgsl -lgslcblas -fopenmp #-Wall to CFLAGS = -lz -lm -O3 -lgsl -lgslcblas #-fopenmp #-Wall
3. zlib (https://zlib.net)
4. make

### BUILD

In the VVP directory:

`make`

Make will build 2 executables:  build_background and VVP

Note:  This has been built and run on Mac laptops and Linux servers.  

## EXAMPLE RUNNING VVP

To see available parameters of the executables, run with the -h option.  

Before running VVP, a background must be built.  From the VVP directory:

`cd example`

`../build_background -i 1KG_cftr_background.recode.vep.vcf.gz -o 1KG.build -b 2500 -v CSQ,4,6,1,15`

The build_background step produces output to stdout for each of the variants in the background vcf file.  It also creates several different output files including extensions .bin, .chr_offsets.txt, .dist.  These files contained information used by VVP.  

To run prioritize variants using VVP (in the example folder):

`../VVP -i target_spiked_simple.vcf.gz -d 1KG.build -v CSQ,4,6,1,15 1> target.spiked.vvp.out`

target_spiked.vvp.out contains the vvp output.

### PREPARE VCF FILE FOR ANALYSIS

The VVP pipeline does not support mulitallelic lines, these must first be decomposed.  We recommend using vt decompose to accomplish this task (http://genome.sph.umich.edu/wiki/Vt).

**Mandatory** preprocessing of a vcf file includes **multiallelic decomposition and VEP annotation**.  It is important to decompose **BEFORE** annnotating because of potential annotation collisions.  Our recommended steps are to use vt to decompose and normalize variants followed by VEP annotation.  No special options in VEP are required for the variant annotation.  Testing has been done with VEP v82.    

## VVP BACKGROUND
A prebuilt background based on gnomAD (http://gnomad.broadinstitute.org/) for use with VVP can be downloaded here (2.5GB): https://s3-us-west-2.amazonaws.com/gnomad-vvp-background/gnomad.062717.build.tar.gz

## VVP OUTPUT
VVP outputs a tab delimited file with 31 columns.  The columns are the following: 

|column name|description|
|-----------|-----------|
|chr| chromosome |
|start| variant start coord|
|ref| reference allele|
|var| variant allele |
|gene| gene id |
|transcript| transcript id |
|hemi_score| raw variant score for hemizygous genotype |
|hemi_vvp| vvp score for hemizygous genotype |
|nhemi| number of hemizygous indivduals |
|hemi_indvs| list of hemizygous individuals |
|hemi_nocall| number of hemizygous nocalls |
|het_score| raw variant score for heterozygous genotype|
|het_vvp| vvp score for heterozygous genotype |
|nhet| number of heterozygous individuals |
|het_indvs| list of heterozygous individuals |
|het_nocall| number of heterozygous nocalls |
|hom_score| raw variant score for homozygous genotype |
|hom_vvp| vvp score for homozygous genotype|
|nhom| number of homozygous individuals |
|hom_indvs| list of homozygous individuals |
|hom_nocall| number of homozygous nocalls |
|coding_ind| 1 if variant is coding, 0 otherwise |
|indel_ind| 1 if variant is an indel, 0 otherwise |
|aa_score| amino acid weight |
|n_bhemi| number of hemizygous background individuals |
|n_bhet| number of heterozygous background individuals |
|n_bhom| number of homozygous background individuals |
|n_bnocall| number of alleles nocalled in background |
|bit_offset| byte offset to background |
|vid| variant id |
|ll_weight| optional extra weight |






