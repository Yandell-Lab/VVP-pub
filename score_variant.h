//
//  score_variant.h
//  vcf_parser
//
//  Created by steven on 6/25/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#ifndef __vcf_parser__score_variant__
#define __vcf_parser__score_variant__

#include "vvp_headers.h"
#include "parse_vcf.h"

void score_variant_b(struct variant * v, int no_allele_frequency);
void score_variant_t_b(struct variant * v, int nb, int xu, int no_allele_frequency); //nb is with nocalls taken into account, xu is the background allele count

#endif /* defined(__vcf_parser__score_variant__) */

