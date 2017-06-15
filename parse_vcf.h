//
//  parse_vcf.h
//  VVP_dev_xcode
//
//  Created by STEVEN FLYGARE on 10/10/16.
//  Copyright © 2016 IDbyDNA. All rights reserved.
//

#ifndef parse_vcf_h
#define parse_vcf_h

#include "vvp_headers.h"
#include "aa_weight.h"

struct vep_field_info {
    uint8_t gene_index;
    uint8_t transcript_index;
    uint8_t seq_ontology_tag_index;
    uint8_t amino_acid_change_index;
    sds annotation_tag_name;
    int ll_weight_index;
};

struct transcript_anno_info {
    char transcript_name[FEATURE_NAME_LENGTH];
    float aaw;
    float het_score;
    int het_vvp;
    float hom_score;
    int hom_vvp;
    float hemi_score;
    int hemi_vvp;
    int coding;
    float llw; //likelihood weight
    sds pref;
    sds pvar;
    kvec_t(sds) anno_tags;
    UT_hash_handle hh;
};

struct gene_transcript {
    char gene_name[FEATURE_NAME_LENGTH];
    struct transcript_anno_info * tai;
    UT_hash_handle hh;
};

struct variant {
    sds chr;
    sds vid;
    size_t pos;
    sds ref;
    sds var;
    int indel;
    struct gene_transcript * gt;
    float phast;
    int nref;  //total number of ref alleles
    int ni; //total number of individuals
    int b_nhet;
    int b_nhom;
    int b_nhemi;
    int b_nocall;
    uint64_t bit_offset;
    sds hemi_indv; //hemizygous individuals (comma separated list)
    sds het_indv; //heterozygous indivdiuals (comma separated list)
    sds hom_indv; //homozygous individuals (comma separated list)
    kvec_t(int) hemi;
    kvec_t(int) hets;
    kvec_t(int) homs;
    kvec_t(int) het_nocalls; //for heterozygous nocalls
    kvec_t(int) hom_nocalls; //for homozygous nocalls
    kvec_t(int) hemi_nocalls; //for hemizygous nocalls
};


void initialize_parse_vcf(uint8_t gene_index, uint8_t transcript_index, uint8_t seq_ontology_tag_index, uint8_t amino_acid_change_index, sds annotation_tag_name, int ll_weight_index);

struct variant * parse_vcf_line(sds line, int no_aa_weight);

void destroy_variant(struct variant * v);

void print_variant(struct variant * v);

#endif /* parse_vcf_h */
