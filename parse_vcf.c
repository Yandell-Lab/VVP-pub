//
//  parse_vcf.c
//  VVP_dev_xcode
//
//  Created by STEVEN FLYGARE on 10/10/16.
//  Copyright © 2016 IDbyDNA. All rights reserved.
//

#include "parse_vcf.h"

#define ERR_GENOTYPE "ERR GENOTYPE"
#define VCF_PROB "VCF FORMAT PROBLEM"
#define ANNO_PROB "ANNOTATION PROBLEM"

static struct vep_field_info vfi;

void initialize_parse_vcf(uint8_t gene_index, uint8_t transcript_index, uint8_t seq_ontology_tag_index, uint8_t amino_acid_change_index, sds annotation_tag_name, int ll_weight_index){
    
    init_aa_score();
    
    vfi.gene_index = gene_index;
    vfi.transcript_index = transcript_index;
    vfi.seq_ontology_tag_index = seq_ontology_tag_index;
    vfi.amino_acid_change_index = amino_acid_change_index;
    vfi.annotation_tag_name = sdsdup(annotation_tag_name);
    vfi.ll_weight_index = ll_weight_index;
    
    sdsfree(annotation_tag_name);
}

void load_gt_info(struct variant ** v, sds * data, int data_len){
    
    int indv_index = 0;
    int i;
    for (i = 9; i < data_len; i++) {
        
        int tmp_count = 0;
        sds * tmp_data = sdssplitlen(data[i], (int)sdslen(data[i]), ":", 1, &tmp_count);
        int gd_count = 0;
        sds * genotype_data = NULL;
        genotype_data = sdssplitlen(tmp_data[0], (int)sdslen(tmp_data[0]), "|", 1, &gd_count);
        if (sdslen(tmp_data[0]) > 1 && gd_count < 2) { //if "|" didn't split, then try "/"
            sdsfreesplitres(genotype_data, gd_count);
            genotype_data = sdssplitlen(tmp_data[0], (int)sdslen(tmp_data[0]), "/", 1, &gd_count);
        }
        sdsfreesplitres(tmp_data, tmp_count);
        
        if (gd_count == 1) { //hemizygous situation
            (*v)->ni++;
            if (strcmp(genotype_data[0], ".") == 0) {
                kv_push(int, (*v)->hemi_nocalls, indv_index);
            }
            else if (strcmp(genotype_data[0], "0") == 0){
                (*v)->nref+=1;
            }
            else if (strcmp(genotype_data[0], "1") == 0){
                kv_push(int, (*v)->hemi, indv_index);
            }
        }
        else if (gd_count == 2){ //diploid call
            (*v)->ni++;
            if (strcmp(genotype_data[0], genotype_data[1]) == 0) { //homozygous call
                if (strcmp(genotype_data[0], ".") == 0) {
                    kv_push(int, (*v)->hom_nocalls, indv_index);
                }
                else if (strcmp(genotype_data[0], "0") == 0){
                    (*v)->nref+=2;
                }
                else if (strcmp(genotype_data[0], "1") == 0){
                    kv_push(int, (*v)->homs, indv_index);
                }
            }
            else { //heterozygous call
                int j;
                for (j = 0; j < 2; j++) {
                    if (strcmp(genotype_data[j], ".") == 0) {
                        kv_push(int, (*v)->het_nocalls, indv_index);
                    }
                    else if (strcmp(genotype_data[j], "0") == 0){
                        (*v)->nref+=1;
                    }
                    else if (strcmp(genotype_data[j], "1") == 0){
                        kv_push(int, (*v)->hets, indv_index);
                    }
                }
                
            }
        }
        else {
            fprintf(stderr, "WARNING:\t%s\tgenotype problem\nchr:%s\tpos:%zu\tcol:%d\n", VCF_PROB, (*v)->chr, (*v)->pos, i);
            //exit(0);
        }
        indv_index++;
        
        sdsfreesplitres(genotype_data, gd_count);
    }
}

void get_aa_change(sds aa_tag, struct transcript_anno_info ** ttai){
    
    if (sdslen(aa_tag) < 1) {
        return;
    }
    int aas = 0;
    sds * aa = sdssplitlen(aa_tag, (int)sdslen(aa_tag), "/", 1, &aas);
    if (aas > 1) {
        (*ttai)->pref = sdsdup(aa[0]);
        (*ttai)->pvar = sdsdup(aa[1]);
    }
    sdsfreesplitres(aa, aas);
}


void check_add_gene_transcript_tags(sds gene_name, sds transcript_name, sds annotation_tags, sds aa_tag, float ll_weight, struct variant ** v) {
    
    struct gene_transcript * tgt = NULL;
    HASH_FIND_STR((*v)->gt, gene_name, tgt); //check for and add gene
    if (tgt == NULL) {
        tgt = (struct gene_transcript *)malloc(sizeof(struct gene_transcript));
        memset(tgt->gene_name, '\0', sizeof(char)*FEATURE_NAME_LENGTH);
        if (sdslen(gene_name) < 1) {
            strncpy(tgt->gene_name, "NONE", 4);
        }
        else {
            strncpy(tgt->gene_name, gene_name, sdslen(gene_name));
        }
        tgt->tai = NULL;
        HASH_ADD_STR((*v)->gt, gene_name, tgt);
    }
    
    struct transcript_anno_info * ttai = NULL;
    HASH_FIND_STR(tgt->tai, transcript_name, ttai); //check for and add transcript
    if (ttai == NULL) {
        ttai = (struct transcript_anno_info *)malloc(sizeof(struct transcript_anno_info));
        memset(ttai->transcript_name, '\0', sizeof(char)*FEATURE_NAME_LENGTH);
        if (sdslen(transcript_name) < 1) {
            strncpy(ttai->transcript_name, "NONE", 4);
        }
        else {
            strncpy(ttai->transcript_name, transcript_name, sdslen(transcript_name));
        }
        
        ttai->aaw = 1.0;
        ttai->llw = ll_weight;
        ttai->het_score = -1.0;
        ttai->hom_score = -1.0;
        ttai->hemi_score = -1.0;
        ttai->het_vvp = -1;
        ttai->hom_vvp = -1;
        ttai->hemi_vvp = -1;
        ttai->coding = 0;
        ttai->pref = NULL;
        ttai->pvar = NULL;
        get_aa_change(aa_tag, &ttai); //get aa weight
        
        kv_init(ttai->anno_tags); //initialize and add annotation tags
        int tags = 0;
        sds * anno_tags = sdssplitlen(annotation_tags, (int)sdslen(annotation_tags), "&", 1, &tags);
        int i;
        for (i=0; i < tags; i++) {
            kv_push(sds, ttai->anno_tags, sdsdup(anno_tags[i]));
        }
        sdsfreesplitres(anno_tags, tags);
        
        HASH_ADD_STR(tgt->tai, transcript_name, ttai);
    }
    else {
        int placeholder = 1;
        //fprintf(stderr, "WARNING:\t%s\nGene:\t%s\nTranscript:\t%s already added.  Multiple annotations per transcript not allowed, will only use first.  Variant location:\t%s\t%zu\n\n", ANNO_PROB, gene_name, transcript_name, (*v)->chr, (*v)->pos);
    }
    
}

void load_annotation_info(sds info_field, struct variant ** v){
    
    int count = 0;
    sds * info_data = sdssplitlen(info_field, (int)sdslen(info_field), ";", 1, &count);
    int i = 0;
    for (i = 0; i < count; i++) {
        int tag_count = 0;
        sds * tag_data = sdssplitlen(info_data[i], (int)sdslen(info_data[i]), "=", 1, &tag_count);
        if (tag_count > 1) {
            if (strcmp(tag_data[0], vfi.annotation_tag_name) == 0) {
                int n_annotations;
                sds * annotations = sdssplitlen(tag_data[1], (int)sdslen(tag_data[1]), ",", 1, &n_annotations);
                int j;
                for (j = 0; j < n_annotations; j++) {
                    int pieces = 0;
                    sds * anno_pieces = sdssplitlen(annotations[j], (int)sdslen(annotations[j]), "|", 1, &pieces);
                    float ll_weight = -1.0;
                    if (vfi.ll_weight_index >= 0) {
                        ll_weight = atof(anno_pieces[vfi.ll_weight_index]);
                    }
                    check_add_gene_transcript_tags(anno_pieces[vfi.gene_index], anno_pieces[vfi.transcript_index], anno_pieces[vfi.seq_ontology_tag_index], anno_pieces[vfi.amino_acid_change_index], ll_weight, v);
                    sdsfreesplitres(anno_pieces, pieces);
                }
                sdsfreesplitres(annotations, n_annotations);
            }
            else if (strcmp(tag_data[0], "PHAST") == 0) {
                    (*v)->phast = atof(tag_data[1]);
            }
        }
        sdsfreesplitres(tag_data, tag_count);
    }
    sdsfreesplitres(info_data, count);
}


struct variant * parse_vcf_line(sds line, int no_aa_weight){
    
    //initialize variant struct values as line is parsed
    struct variant * v = (struct variant *)malloc(sizeof(struct variant));
    v->chr = NULL;
    v->bit_offset = 0;
    sdstrim(line, "\n");
    
    //first split line by tab
    int count = 0;
    sds * data = sdssplitlen(line, (int)sdslen(line), "\t", 1, &count);
    
    if (count < 10) {
        fprintf(stderr, "%s:\t vcf line has fewer than 10 columns, will skip:\t%s:%s, %s, %s\n", VCF_PROB, data[0], data[1], data[4], data[7]);
        sdsfreesplitres(data, count);
        free(v);
        return NULL;
    }
    
    v->chr = sdsdup(data[0]);
    v->pos = atoll(data[1]);
    v->vid = sdsdup(data[2]);
    v->ref = sdsdup(data[3]);
    
    //now find number of alternate alleles, if > 1, then throw error
    int num_va = 0;
    sds * variant_alleles = sdssplitlen(data[4], (int)sdslen(data[4]), ",", 1, &num_va);
    if (num_va > 1) {
        fprintf(stderr, "FATAL:\t vcf line has >1 alternate alleles.  File needs to be decomposed before processing\n");
        fprintf(stderr, "%s", line);
        sdsfreesplitres(variant_alleles, num_va);
        sdsfreesplitres(data, count);
        exit(0);
    }
    
    v->var = sdsdup(variant_alleles[0]);
    sdsfreesplitres(variant_alleles, num_va);
    //indel indicator
    if ((sdslen(v->ref) != sdslen(v->var)) || strcmp(v->ref, "-") == 0 || strcmp(v->var, "-") == 0) {
        v->indel = 1;
    }
    
    //load genotype information
    v->ni = 0;
    v->nref = 0;
    kv_init(v->hets);
    v->het_indv = sdsempty();
    kv_init(v->homs);
    v->hom_indv = sdsempty();
    kv_init(v->hemi);
    v->hemi_indv = sdsempty();
    kv_init(v->het_nocalls);
    kv_init(v->hom_nocalls);
    kv_init(v->hemi_nocalls);
    load_gt_info(&v, data, count);
    
    //load annotation and protein information
    v->gt = NULL;
    v->phast = -1.0;
    load_annotation_info(data[7], &v);
    
    //add amino acid weights
    struct gene_transcript * c, * t;
    HASH_ITER(hh, v->gt, c, t) {
        struct transcript_anno_info * current, * tmp;
        HASH_ITER(hh, c->tai, current, tmp) {
            get_aaw(&current, v->ref, v->var, v->phast);
            if (no_aa_weight == 1) {
                current->aaw = 1.0;
            }
        }
    }

    sdsfreesplitres(data, count);
    
    return v;
}

struct variant * parse_allele_frequency_line(sds line, int no_aa_weight){

    //initialize variant struct values as line is parsed
    struct variant * v = (struct variant *)malloc(sizeof(struct variant));
    v->chr = NULL;
    v->bit_offset = 0;
    sdstrim(line, "\n");
    
    //first split line by tab
    int count = 0;
    sds * data = sdssplitlen(line, (int)sdslen(line), "\t", 1, &count);
    
    if (count < 10) {
        fprintf(stderr, "%s:\t vcf line has fewer than 10 columns, will skip:\t%s:%s, %s, %s\n", VCF_PROB, data[0], data[1], data[4], data[7]);
        sdsfreesplitres(data, count);
        free(v);
        return NULL;
    }
    
    v->chr = sdsdup(data[0]);
    v->pos = atoll(data[1]);
    v->vid = sdsdup(data[2]);
    v->ref = sdsdup(data[3]);
    
    //now find number of alternate alleles, if > 1, then throw error
    int num_va = 0;
    sds * variant_alleles = sdssplitlen(data[4], (int)sdslen(data[4]), ",", 1, &num_va);
    if (num_va > 1) {
        fprintf(stderr, "FATAL:\t vcf line has >1 alternate alleles.  File needs to be decomposed before processing\n");
        fprintf(stderr, "%s", line);
        sdsfreesplitres(variant_alleles, num_va);
        sdsfreesplitres(data, count);
        exit(0);
    }
    
    v->var = sdsdup(variant_alleles[0]);
    sdsfreesplitres(variant_alleles, num_va);
    //indel indicator
    if ((sdslen(v->ref) != sdslen(v->var)) || strcmp(v->ref, "-") == 0 || strcmp(v->var, "-") == 0) {
        v->indel = 1;
    }
    
    //load genotype information
    
    v->ni = 0;
    v->nref = 0;
    kv_init(v->hets);
    v->het_indv = sdsempty();
    kv_init(v->homs);
    v->hom_indv = sdsempty();
    kv_init(v->hemi);
    v->hemi_indv = sdsempty();
    kv_init(v->het_nocalls);
    kv_init(v->hom_nocalls);
    kv_init(v->hemi_nocalls);
    
    //load_gt_info(&v, data, count);
    int total_called = atoi(data[5]);
    int n_alleles_called = atoi(data[6]);
    int nocalled_alleles = total_called - n_alleles_called; //number of nocall alleles
    int n_alleles_alt = atoi(data[7]);
    int n_hom_indvs = atoi(data[8]);
    
    int n_hemi = 0;
    int n_het = 0;
    if (strcmp(v->chr, "X") == 0) {
        n_hemi = atoi(data[9]); //hemizgyous will be the same as the number of male alleles
        int tmp = n_alleles_alt - ((n_hom_indvs * 2) + n_hemi);
        n_het = tmp >= 0 ? tmp : 0; //whatever is left will be het
    }
    else {
        n_het = n_alleles_alt - (n_hom_indvs *2); //number of heterozygous alleles -- same as het individuals
    }
    
    
    int n_hom_alt = n_alleles_alt - (n_het + n_hemi); //number of alleles homozygous for alt allele
    
    v->nref = n_alleles_called - (n_hom_alt + n_het + n_hemi);
    
    //generate random numbers for the individual indexes
    gsl_rng_env_setup();
    gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
    
    int i;
    for (i = 0; i < (n_hom_alt / 2); i++) {
        int z = (int)gsl_rng_uniform_int(r, (total_called/2)+1);
        kv_push(int, v->homs, z);
    }
    for (i = 0; i < n_het; i++) {
        int z = (int)gsl_rng_uniform_int(r, (total_called/2)+1);
        kv_push(int, v->hets, z);
    }
    for (i = 0; i < n_hemi; i++) {
        int z = (int)gsl_rng_uniform_int(r, (total_called/2)+1);
        kv_push(int, v->hemi, z);
    }
    for (i = 0; i < nocalled_alleles; i++) {
        int z = (int)gsl_rng_uniform_int(r, (total_called/2)+1);
        kv_push(int, v->het_nocalls, z);
    }
    
    gsl_rng_free(r);
    
    //load annotation and protein information
    v->gt = NULL;
    v->phast = -1.0;
    load_annotation_info(data[10], &v);
    
    //add amino acid weights
    struct gene_transcript * c, * t;
    HASH_ITER(hh, v->gt, c, t) {
        struct transcript_anno_info * current, * tmp;
        HASH_ITER(hh, c->tai, current, tmp) {
            get_aaw(&current, v->ref, v->var, v->phast);
            if (no_aa_weight == 1) {
                current->aaw = 1.0;
            }
        }
    }
    
    sdsfreesplitres(data, count);
    
    return v;

}



void destroy_variant(struct variant * v){
    sdsfree(v->chr);
    sdsfree(v->vid);
    sdsfree(v->ref);
    sdsfree(v->var);
    kv_destroy(v->hets);
    sdsfree(v->het_indv);
    kv_destroy(v->homs);
    sdsfree(v->hom_indv);
    kv_destroy(v->hemi);
    sdsfree(v->hemi_indv);
    kv_destroy(v->het_nocalls);
    kv_destroy(v->hom_nocalls);
    kv_destroy(v->hemi_nocalls);
    
    struct gene_transcript * c, * t;
    HASH_ITER(hh, v->gt, c, t) {
        struct transcript_anno_info * current, * tmp;
        HASH_ITER(hh, c->tai, current, tmp) {
            int i;
            if(current->pref != NULL) sdsfree(current->pref);
            if(current->pvar != NULL) sdsfree(current->pvar);
            size_t tmpl = c->tai->anno_tags.n;
            for (i=0; i < tmpl; i++) {
                sdsfree(kv_pop(c->tai->anno_tags));
            }
            kv_destroy(c->tai->anno_tags);
            HASH_DEL(c->tai, current);
            free(current);
        }
        HASH_DEL(v->gt, c);
        free(c);
    }
    free(v);
    v = NULL;
}

void print_vec_comma(kvec_t(int) * tmp){
    int i;
    for (i = 0; i < (*tmp).n; i++) {
        fprintf(stdout, "%d,", kv_A(*tmp, i));
    }
}

void print_variant(struct variant * v){
    fprintf(stdout, "chr: %s, pos: %zu\n", v->chr, v->pos);
    fprintf(stdout, "id: %s\n", v->vid);
    fprintf(stdout, "ref: %s, var: %s\n", v->ref, v->var);
    fprintf(stdout, "number of individuals: %d\n", v->ni);
    fprintf(stdout, "number ref alleles: %d\n", v->nref);
    
    int i;
    if (v->hemi.n > 0) {
        fprintf(stdout, "HEMIZYGOUS individuals: %zu\n", v->hemi.n);
        //for (i = 0; i < v->hemi.n; i++) fprintf(stdout, "%d,", kv_A(v->hemi, i));
        //fprintf(stdout, "\n\n");
    }
    if (v->hets.n > 0) {
        fprintf(stdout, "HETEROZYGOUS individuals: %zu\n", v->hets.n);
        //for (i = 0; i < v->hets.n; i++) fprintf(stdout, "%d,", kv_A(v->hets, i));
        //fprintf(stdout, "\n\n");
    }
    if (v->homs.n > 0) {
        fprintf(stdout, "HOMOZYGOUS individuals: %zu\n", v->homs.n);
        //for (i = 0; i < v->homs.n; i++) fprintf(stdout, "%d,", kv_A(v->homs, i));
        //fprintf(stdout, "\n\n");
    }
    if (v->hemi_nocalls.n > 0) {
        fprintf(stdout, "HEMIZYGOUS nocalled individuals: %zu\n", v->hemi_nocalls.n);
        //for (i = 0; i < v->hemi_nocalls.n; i++) fprintf(stdout, "%d,", kv_A(v->hemi_nocalls, i));
        //fprintf(stdout, "\n\n");
    }
    if (v->het_nocalls.n > 0) {
        fprintf(stdout, "HETEROZYGOUS nocalled individuals: %zu\n", v->het_nocalls.n);
        //for (i = 0; i < v->het_nocalls.n; i++) fprintf(stdout, "%d,", kv_A(v->het_nocalls, i));
        //fprintf(stdout, "\n\n");
    }
    if (v->hom_nocalls.n > 0) {
        fprintf(stdout, "HOMOZYGOUS nocalled individuals: %zu\n", v->hom_nocalls.n);
        //for (i = 0; i < v->hom_nocalls.n; i++) fprintf(stdout, "%d,", kv_A(v->hom_nocalls, i));
        //fprintf(stdout, "\n\n");
    }
    
    struct gene_transcript * c, * t;
    HASH_ITER(hh, v->gt, c, t) {
        fprintf(stdout, "GENE: %s\n", c->gene_name);
        struct transcript_anno_info * current, * tmp;
        HASH_ITER(hh, c->tai, current, tmp) {
            fprintf(stdout, "\tTRANSCRIPT: %s\n", current->transcript_name);
            fprintf(stdout, "\t\taaw: %f, het_score: %f, hom_score: %f, coding: %d\n", current->aaw, current->het_score, current->hom_score, current->coding);
            size_t tmpl = current->anno_tags.n;
            for (i = 0; i < tmpl; i++) {
                fprintf(stdout, "\t\tSO TAG: %s\n", kv_A(current->anno_tags, i));
            }
        }
    }
    
}


