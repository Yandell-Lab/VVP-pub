//
//  background_max_scores.h
//  VVP_C
//
//  Created by steven on 8/13/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#ifndef __VVP_C__background_max_scores__
#define __VVP_C__background_max_scores__

#include "vvp_headers.h"

#define BKRND_GENE_NAME_LEN 50

struct bkgrnd_max_scores {
    char gene[BKRND_GENE_NAME_LEN];
    float * max_scores;
    UT_hash_handle hh;
};

void init_bkgrnd_max(char * bkgrnd_max, int nb, char iht);

void init_bkgrnd_max_b(char * bkgrnd_max, int nb, char iht);

struct bkgrnd_max_scores * get_gene_max(char * gene);

void cleanup_bkgrnd_max();

#endif /* defined(__VVP_C__background_max_scores__) */

