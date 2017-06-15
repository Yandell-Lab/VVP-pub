//
//  aa_weight.h
//  VVP_dev_xcode
//
//  Created by STEVEN FLYGARE on 10/11/16.
//  Copyright © 2016 IDbyDNA. All rights reserved.
//

#ifndef aa_weight_h
#define aa_weight_h

#include "vvp_headers.h"
#include "parse_vcf.h"

struct aa_matrix {
    char aa_change[20];
    float score;
    float cons;
    float uncons;
    UT_hash_handle hh;
};

void init_aa_score();
void get_aaw(struct transcript_anno_info ** ttai, sds ref, sds var, float phast);

#endif /* aa_weight_h */
