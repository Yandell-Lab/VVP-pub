//
//  vvp_lookup.h
//  VVP_C
//
//  Created by steven on 8/12/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#ifndef __VVP_C__vvp_lookup__
#define __VVP_C__vvp_lookup__

#include "vvp_headers.h"

#define NPERCENTILES 100

struct feature_lookups {
    char feature_name[FEATURE_NAME_LEN];
    float coding_vals[NPERCENTILES];
    float noncoding_vals[NPERCENTILES];
    int n_coding;
    int n_noncoding;
    UT_hash_handle hh;
};

void load_feature_lookups_b(sds lookup_file);

int score_lookup(char * feature_name, float score, int coding);

int score_lookup_b(char * feature_name, float score, int coding);

void destroy_feature_lookups();


#endif /* defined(__VVP_C__vvp_lookup__) */
