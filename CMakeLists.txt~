cmake_minimum_required(VERSION 2.8)
project(VVP)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lz -lm -Wall -lgsl -fopenmp") 

set(BB_SOURCE_FILES build_background.c score_variant.c sds.c aa_weight.c parse_vcf.c bit_array.c)
set(VVP_SOURCE_FILES score_variants.c score_variant.c sds.c aa_weight.c parse_vcf.c vvp_lookup.c search_binary_bkgrnd.c)

add_executable(build_background ${BB_SOURCE_FILES})
add_executable(VVP ${VVP_SOURCE_FILES})
