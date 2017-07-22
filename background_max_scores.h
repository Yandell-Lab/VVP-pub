//
//  background_max_scores.c
//  VVP_C
//
//  Created by steven on 8/13/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#include "background_max_scores.h"

static struct bkgrnd_max_scores * bms;

void get_indv_scores(struct bkgrnd_max_scores ** tbms, sds indv_scores, char iht){
    
    int count = 0;
    sds * data = sdssplitlen(indv_scores, (int)sdslen(indv_scores), ";", 1, &count);
    int i=0;
    int indv = 0;
    float hom = 0;
    float het1 = 0;
    float het2 = 0;
    float sum = 0;
    for (i=0; i < count; i++) {
        sscanf(data[i], "%d:%f,%f,%f,%f", &indv, &hom, &het1, &het2, &sum);
        if (iht == 'r') { //recessive
            if ((hom > (het1 + het2)) || (het1 <= 0 || het2 <= 0)) {
                (*tbms)->max_scores[indv] = hom;
            }
            else {
                (*tbms)->max_scores[indv] = (het1 + het2);
            }
        }
        else if (iht == 'd') { //dominant
            if (het1 > het2) {
                (*tbms)->max_scores[indv] = het1;
            }
            else {
                (*tbms)->max_scores[indv] = het2;
            }
        }
        else { //no inheritance, choose sum
            (*tbms)->max_scores[indv] = sum;
        }
    }
    
    sdsfreesplitres(data, count);
}

void get_indv_scores_b(struct bkgrnd_max_scores ** tbms, float * scores, int nb, char iht){
    
    int i=0;
    int j=0;
    float hom, het1, het2, sum;
    for (i=0; i < 4*nb; i+=4) {
        hom = scores[i];
        het1 = scores[i+1];
        het2 = scores[i+2];
        sum = scores[i+3];
        
        if (iht == 'r' || iht == 'x') { //recessive
            if (hom > (het1 + het2) || (het1 <= 0 || het2 <= 0)) {
                (*tbms)->max_scores[j] = hom;
            }
            else {
                (*tbms)->max_scores[j] = (het1 + het2);
            }
        }
        else if (iht == 'd') { //dominant
            /*if (hom > het1 && hom > het2) {
                (*tbms)->max_scores[j] = hom;
            }*/
            if (het1 > het2) {
                (*tbms)->max_scores[j] = het1;
            }
            else {
                (*tbms)->max_scores[j] = het2;
            }
        }
        else { //no inheritance, choose sum
            (*tbms)->max_scores[j] = sum;
        }
        j++;
    }
}



#define MAX_BUF 1000000

void init_bkgrnd_max(char * bkgrnd_max, int nb, char iht){
    bms = NULL;
    
    FILE * max_in = fopen(bkgrnd_max, "r");
    if (! max_in) {
        fprintf(stderr, "FATAL: could not open %s for loading\n", bkgrnd_max);
        exit(1);
    }
    
    int line_count = 0;
    
    char * buffer = malloc(sizeof(char)*MAX_BUF);
    while ( fgets(buffer, MAX_BUF, max_in) != NULL)  {
        line_count += 1;
        if (line_count % 1000 == 0) {
            fprintf(stderr, "%d,", line_count);
        }
        sds tmpl = sdsnew(buffer);
        sdstrim(tmpl, "\n");
        int count = 0;
        sds * data = sdssplitlen(tmpl, (int)sdslen(tmpl), "\t", 1, &count);
        if (count != 2) {
            fprintf(stderr, "WARNING: line in max_score wrong format, will be skipped: %s", tmpl);
            sdsfreesplitres(data, count);
            sdsfree(tmpl);
            continue;
        }
        
        struct bkgrnd_max_scores * tbms = (struct bkgrnd_max_scores *)calloc(1, sizeof(struct bkgrnd_max_scores));
        strcpy(tbms->gene, data[0]);
        tbms->max_scores = (float *)calloc(nb, sizeof(float));
        get_indv_scores(&tbms, data[1], iht);
        
        sdsfree(tmpl);
        sdsfreesplitres(data, count);
        
        HASH_ADD_STR(bms, gene, tbms);
    }
    free(buffer);
}

void init_bkgrnd_max_b(char * bkgrnd_max, int nb, char iht){
    bms = NULL;
    
    FILE * max_in = fopen(bkgrnd_max, "rb");
    if (! max_in) {
        fprintf(stderr, "FATAL: could not open %s for loading\n", bkgrnd_max);
        exit(1);
    }
    
    char feature[BKRND_GENE_NAME_LEN];
    float * scores = (float *)malloc(sizeof(float)*4*nb);
    int line_count = 0;
    
    while ( fread(&feature, sizeof(char), BKRND_GENE_NAME_LEN, max_in) )  {
        line_count += 1;
        if (line_count % 10000 == 0) {
            fprintf(stderr, "%d,", line_count);
        }
        //read in float data
        memset(scores, '\0', sizeof(float)*4*nb);
        fread(scores, sizeof(float), 4*nb, max_in);
        
        struct bkgrnd_max_scores * tbms = (struct bkgrnd_max_scores *)calloc(1, sizeof(struct bkgrnd_max_scores));
        strcpy(tbms->gene, feature);
        tbms->max_scores = (float *)calloc(nb, sizeof(float));
        get_indv_scores_b(&tbms, scores, nb, iht);
        
        /* debug info
        if (strcmp(feature, "ENSG00000130283") == 0){
            fprintf(stderr, "\nPRINTING BACKGROUND SCORES\n\n%s\n", feature);
            int scores_index = 0;
            int i=0;
            for (i=0; i < nb; i++) {
                if (tbms->max_scores[i] > 0) {
                    int j = 0;
                    fprintf(stderr, "\t%d", i);
                    for (j = 0; j < 4; j++) {
                        fprintf(stderr, "\t%f", scores[scores_index+j]);
                    }
                    fprintf(stderr, "\t%f", tbms->max_scores[i]);
                    fprintf(stderr, "\n");
                }
                scores_index += 4;
            }
        }
        */
        
        
        
        HASH_ADD_STR(bms, gene, tbms);
        
    }
    
    free(scores);
    
}


struct bkgrnd_max_scores * get_gene_max(char * gene){
    
    struct bkgrnd_max_scores * t = NULL;
    HASH_FIND_STR(bms, gene, t);
    return t;
    
}

void cleanup_bkgrnd_max(){
    
    struct bkgrnd_max_scores *s, *tmp;
    HASH_ITER(hh, bms, s, tmp){
        HASH_DEL(bms, s);
        free(s);
    }

}


