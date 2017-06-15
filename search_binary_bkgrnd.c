//
//  search_binary_bkgrnd.c
//  VVP_C
//
//  Created by steven on 8/11/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#include "search_binary_bkgrnd.h"

static uint64_t mm_size;

unsigned char * load_bin_db(sds file_prefix, int * n_background){
    
    unsigned char * mm_bin = NULL;
    
    sds bin_file = sdsdup(file_prefix);
    bin_file = sdscat(bin_file, ".bin");
    
    int fdSrc = open(bin_file, O_RDWR, 0);
    struct stat st;
    fstat(fdSrc, &st);
    
    mm_bin = (unsigned char *)mmap(NULL, st.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fdSrc, 0);
    if(mm_bin == MAP_FAILED){
        fprintf(stderr, "FATAL:  could not create mmap from %s\n", bin_file);
        exit(1);
    }
    
    //load memmap into memory
    /*
    size_t filesize = st.st_size;
    size_t page_size = getpagesize();
    unsigned char * buf[page_size];
    size_t pos=0;
    for (pos=0; pos < filesize; pos += page_size) {
        size_t this_page_size = filesize - pos;
        if (this_page_size > page_size){
            this_page_size = page_size;
        }
        memcpy(buf, mm_bin + pos, this_page_size);
    }
    */
    
    mm_size = st.st_size;
    fprintf(stderr, "MMAP size for .bin: %llu\n\n", mm_size);
    
    sdsfree(bin_file);
    
    *n_background = *((int *)mm_bin);
    
    return mm_bin;
}

unsigned char * load_bit_db(sds file_prefix){
    unsigned char * mm_bits = NULL;
    
    sds bit_file = sdsdup(file_prefix);
    bit_file = sdscat(bit_file, ".bit");
    
    int fdSrc = open(bit_file, O_RDWR, 0);
    struct stat st;
    fstat(fdSrc, &st);
    mm_bits = (unsigned char *)mmap(NULL, st.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fdSrc, 0);
    if(mm_bits == MAP_FAILED){
        fprintf(stderr, "FATAL:  could not create mmap from %s\n", bit_file);
        exit(1);
    }
    
    sdsfree(bit_file);
    
    return mm_bits;

}

struct chr_offsets * load_offsets(sds file_prefix){
    
    sds offset_file = sdsdup(file_prefix);
    offset_file = sdscat(offset_file, ".chr_offsets.txt");
    
    FILE * offsets = fopen(offset_file, "r");
    if (! offsets) {
        fprintf(stderr, "FATAL: could not open offsets file %s\n", offset_file);
        exit(1);
    }
    
    char chr[3];
    memset(chr, '\0', 3);
    uint64_t start;
    uint64_t end;
    int n;
    struct chr_offsets * chro = NULL;

    
    while (fscanf(offsets, "%s\t%llu\t%llu\t%d", chr, &start, &end, &n) == 4) {
        
        struct chr_offsets * tc = (struct chr_offsets *)calloc(1, sizeof(struct chr_offsets));
        //tc->chr = (char *)calloc(3, sizeof(char));
        strcpy(tc->chr, chr);
        tc->byte_start = start;
        tc->byte_end = end;
        tc->n_entries = n;
        memset(chr, '\0', 3);
        
        struct chr_offsets * ttc = NULL;
        HASH_FIND_STR(chro, tc->chr, ttc);
        //HASH_FIND(hh, chro, &tc->chr, 3*sizeof(char), ttc);
        if (ttc != NULL) {
            fprintf(stderr, "FATAL: chromosome %s already seen in offsets", chr);
            exit(1);
        }
        else {
            HASH_ADD_STR(chro, chr, tc);
            //HASH_ADD(hh, chro, chr, 3*sizeof(char), tc);
        }
        
    }
    
    sdsfree(offset_file);
    fclose(offsets);
    
    return chro;
    
}

void get_variants_S(struct m_var_info ** vi, unsigned char * mm_bin, uint64_t mm_offset){ // "_S" means side effects
    
    (*vi)->vi[(*vi)->nv] = (struct var_info *)malloc(sizeof(struct var_info));
    mm_offset += 4; //skip start position
    (*vi)->vi[(*vi)->nv]->var_type = *( (char *)(mm_bin + mm_offset) );
    mm_offset += 1;
    (*vi)->vi[(*vi)->nv]->length = *( (int *)(mm_bin + mm_offset) );
    mm_offset += 4;
    (*vi)->vi[(*vi)->nv]->nhet = *( (int *)(mm_bin + mm_offset) );
    mm_offset += 4;
    (*vi)->vi[(*vi)->nv]->nhom = *( (int *)(mm_bin + mm_offset) );
    mm_offset += 4;
    (*vi)->vi[(*vi)->nv]->nhemi = *( (int *)(mm_bin + mm_offset) );
    mm_offset += 4;
    (*vi)->vi[(*vi)->nv]->nocall = *( (int *)(mm_bin + mm_offset) );
    mm_offset += 4;
    (*vi)->vi[(*vi)->nv]->bit_offset = *( (uint64_t *)(mm_bin + mm_offset) );
    (*vi)->nv++; //increment number of variants
}


struct m_var_info * search_binary_bkgrnd(char * chr, size_t pos, unsigned char * mm_bin, struct chr_offsets * chro){
    
    struct m_var_info * vi = (struct m_var_info *)malloc(sizeof(struct m_var_info *));
    vi->vi = (struct var_info **)malloc(sizeof(struct var_info *)*MAX_VARS);
    vi->nv = 0;
    
    struct chr_offsets * ttc = NULL;
    HASH_FIND_STR(chro, chr, ttc);
    //HASH_FIND(hh, chro, &chr, 3*sizeof(char), ttc);
    if (ttc == NULL) {
        fprintf(stderr, "WARNING: chromosome %s not in offsets\n", chr);
        return vi;
    }
    else {
        
        int min = 0;
        int max = ttc->n_entries - 1; //place pointer at start of last entry
        while (max >= min) { //binary search
            int mid = (min + max) / 2;
            uint64_t tmp_offset = ttc->byte_start + mid*LINE_BYTE_SIZE + 2; //+2 because only two chars were written for chromosome placeholder
            int * mm_pos = (int *)(mm_bin + tmp_offset);
            if ( (*mm_pos) == pos) {
                get_variants_S(&vi, mm_bin, tmp_offset);
                
                //look 'down' from found position
                uint64_t d_offset = tmp_offset;
                while (d_offset >= LINE_BYTE_SIZE && *( (int *)(mm_bin + (d_offset - LINE_BYTE_SIZE)) ) == pos) {
                    d_offset -= LINE_BYTE_SIZE;
                    get_variants_S(&vi, mm_bin, d_offset);
                }
                
                //look 'up' from found position
                uint64_t u_offset = tmp_offset + LINE_BYTE_SIZE;
                while (u_offset < mm_size && *( (int *)(mm_bin+u_offset) ) == pos) {
                    get_variants_S(&vi, mm_bin, u_offset);
                    u_offset += LINE_BYTE_SIZE;
                }
                
                return vi;
                
            }
            else if (*mm_pos < pos) {
                min = mid + 1;
            }
            else {
                max = mid - 1;
            }
        }
    }
    
    
    return vi;
}

void destroy_chr_offsets(struct chr_offsets * chro){
    struct chr_offsets *s, *tmp;
    HASH_ITER(hh, chro, s, tmp){
        HASH_DEL(chro, s);
        free(s);
    }
}

