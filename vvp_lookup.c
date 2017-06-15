//
//  vvp_lookup.c
//  VVP_C
//
//  Created by steven on 8/12/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#include "vvp_lookup.h"

static struct feature_lookups * lookups;

static unsigned char * mm_dist;
static uint64_t mm_dist_size;
static int dist_line_size;

void load_feature_lookups_b(sds lookup_file) {
    
    int fdSrc = open(lookup_file, O_RDWR, 0);
    struct stat st;
    fstat(fdSrc, &st);
    mm_dist = (unsigned char *)mmap(NULL, st.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fdSrc, 0);
    if(mm_dist == MAP_FAILED){
        fprintf(stderr, "FATAL:  could not create mmap from %s\n", lookup_file);
        exit(1);
    }
    
    mm_dist_size = st.st_size;
    
    dist_line_size = sizeof(char)*FEATURE_NAME_LEN + sizeof(float)*NPERCENTILES + sizeof(float)*NPERCENTILES + sizeof(size_t) + sizeof(size_t);

}

int score_lookup_b(char * feature_name, float score, int coding){
    
    int min = 0;
    int max = mm_dist_size / dist_line_size; //place pointer at start of final line
    float * percentiles = NULL;
    while (max >= min) {
        //uint mid = min + ((max - min) >> 1); //floor average
        int mid = (min + max) / 2;
        uint64_t tmp_offset = mid*dist_line_size;
        int cmp = strcmp((char *)(mm_dist+tmp_offset), feature_name);
        if ( cmp == 0 ) {
            if (coding == 1) {
                percentiles = (float *)(mm_dist + tmp_offset + sizeof(char)*FEATURE_NAME_LEN);
            }
            else if(coding == 0) {
                percentiles = (float *)(mm_dist + tmp_offset + sizeof(char)*FEATURE_NAME_LEN + sizeof(float)*NPERCENTILES);
            }
            else {
                return -1;
            }
            int i=0;
            for (i = 0; i < NPERCENTILES; i++) {
                if ( score < percentiles[i] || fabsf(score - percentiles[i]) < .01 ) { //return if within .01
                    return i;
                }
            }
            return 100;
        }
        else if (cmp < 0){
            min = mid + 1;
        }
        else {
            max = mid - 1;
        }
    }
    
    return -1;
    
}

void destroy_feature_lookups(){
    struct feature_lookups *s, *tmp;
    HASH_ITER(hh, lookups, s, tmp){
        HASH_DEL(lookups, s);
        free(s);
    }
}
