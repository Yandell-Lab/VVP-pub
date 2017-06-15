//
//  vvp_headers.h
//  VVP_C
//
//  Created by steven on 8/12/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#ifndef VVP_C_vvp_headers_h
#define VVP_C_vvp_headers_h

#include "sds.h"
#include "uthash.h"
#include "kvec.h"
#include "khash.h"
#include "bit_array.h"
#include <math.h>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include <pthread.h>

struct config {
    sds target_vcf;
    sds background_vcf;
    sds db_prefix;
    sds vvp_formatted;
    char inheritance_filters;
    char penetrance;
    int mother_sample_index;
    int father_sample_index;
    int proband_sample_index;
    int sibling_sample_index;
    int sibling_affected;
    int nb;
    int nts;
    int nt;
    sds anno_tag;
    int variant_pos;
    int gene_pos;
    int aa_pos;
    int so_pos;
    char iht;
    int format_output;
    int np;
    int only_coding;
    int only_snv;
    size_t n_permutations;
    size_t mat_rows;
};

#ifdef _OPENMP
#include <omp.h>
#endif

#define BUF_SIZE 5000000
#define MAX_SCORES 100000
#define LINE_BYTE_SIZE 35
#define NPERCENTILES 100
#define MAX_VARS 20
#define FEATURE_NAME_LEN 50
#define FEATURE_NAME_LENGTH 50

#endif
