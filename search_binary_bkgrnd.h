//
//  search_binary_bkgrnd.h
//  VVP_C
//
//  Created by steven on 8/11/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#ifndef __VVP_C__search_binary_bkgrnd__
#define __VVP_C__search_binary_bkgrnd__

#include "vvp_headers.h"

struct chr_offsets {
    char chr[3];
    //char * chr;
    uint64_t byte_start;
    uint64_t byte_end;
    int n_entries;
    UT_hash_handle hh;
};

struct var_info {
    char var_type;
    int length;
    int nhet;
    int nhom;
    int nhemi;
    int nocall;
    uint64_t bit_offset;
};

struct m_var_info {
    struct var_info ** vi;
    int nv;
};

unsigned char * load_bin_db(sds file_prefix, int * n_background);

unsigned char * load_bit_db(sds file_prefix);

struct chr_offsets * load_offsets(sds file_prefix);

struct m_var_info * search_binary_bkgrnd(char * chr, size_t pos, unsigned char * mm_bin, struct chr_offsets * chro);

void destroy_chr_offsets(struct chr_offsets * chro);

#endif /* defined(__VVP_C__search_binary_bkgrnd__) */
