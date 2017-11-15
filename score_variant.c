//
//  score_variant.c
//  vcf_parser
//
//  Created by steven on 6/25/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#include "score_variant.h"

float compute_score(int nb, int nt, int xa, int xu, float aaw, float llw, int no_allele_frequency) {
    
    //nb = nb - nt >= 0 ? nb - nt : 0;
    //xu = xu - xa >= 0 ? xu - xa : 0;
    
    int x = xa + xu;
    int n = nb + nt;
    
    float p = (float)x / (float) n;
    
    if (p < 1e-10 || (1.0 - p) < 1e-10) { //if p is essentially 0 or 1, the score is 0
        return 0.0;
    }
    
    float pu = 0.0;
    if (nb > 0) {
        pu = (float)xu/(float)nb;
        if (pu >= 1.0) { //if everyone in the background has the allele, the score is 0
            return 0.0;
        }
        else if (pu < 1e-10){
            pu = 1e-6;
        }
    }
    if (pu < 1e-10) { //in case everyone in the background is nocalled
        pu = 1e-6;
    }
    
    float pa = (float)xa / (float)nt;
    if (pa >= 1.0) {
        pa = 1.0 - 1e-6;
    }
    else if (pa <= 1e-10) { //error if pa is 0 or negative -- no affecteds with allele
        return -1.0;
    }
    
    double plog = log(p);
    /*if (errno == EDOM || errno == ERANGE) {
     fprintf(stderr, "log(p) failed, p is %f\n", p);
     free(vac);
     return -1.0; //return -1 to mean an error
     }*/
    
    double iplog = log(1.0 - p);
    /*if (errno == EDOM || errno == ERANGE) {
     fprintf(stderr, "log(1.0 - p) failed, p is %f\n", p);
     free(vac);
     return -1.0; //return -1 to mean an error
     }*/
    
    double pulog = log(pu);
    /*if (errno == EDOM || errno == ERANGE) {
     fprintf(stderr, "log(pu) failed, pu is %f\n", pu);
     free(vac);
     return -1.0; //return -1 to mean an error
     }*/
    
    double ipulog = log(1.0 - pu);
    /*if (errno == EDOM || errno == ERANGE) {
     fprintf(stderr, "log(1.0 - pu) failed, pu is %f\n", pu);
     free(vac);
     return -1.0; //return -1 to mean an error
     }*/
    
    double palog = log(pa);
    /*if (errno == EDOM || errno == ERANGE) {
     fprintf(stderr, "log(pa) failed, pa is %f\n", pa);
     free(vac);
     return -1.0; //return -1 to mean an error
     }*/
    
    double ipalog = log(1.0 - pa);
    /*if (errno == EDOM || errno == ERANGE) {
     fprintf(stderr, "log(1.0 - pa) failed, pa is %f\n", pu);
     free(vac);
     return -1.0; //return -1 to mean an error
     }*/
    
    double aalog = log(aaw);
    /*if (errno == EDOM || errno == ERANGE) {
     fprintf(stderr, "log(tv->aaw) failed, tv->aaw is %f\n", tv->aaw);
     free(vac);
     return -1.0; //return -1 to mean an error
     }*/
    
    if (llw < 0) {
        llw = 1.0;
    }
    else if (llw <= 1e-10) {
        llw = 1e-6;
    }
    float log_llw = log(llw);
    
    float numerator = x*plog + (n-x)*iplog;
    float denominator = xu*pulog + (nb - xu)*ipulog + xa*palog + (nt - xa)*ipalog;
    float diff = no_allele_frequency == 0 ? (numerator - denominator) : 0.0;
    float score = -2.0*(log_llw + aalog + diff);
    if (score <= 0.0) {
        return 0.0;
    }
    
    return score;

}

void score_variant_b(struct variant * v, int no_allele_frequency){
    
    int nb = v->nref + v->hemi.n + v->hets.n + 2*(v->homs.n);
    int xu = nb - v->nref;
    
    struct gene_transcript * c, * t;
    HASH_ITER(hh, v->gt, c, t) {
        struct transcript_anno_info * current, * tmp;
        HASH_ITER(hh, c->tai, current, tmp) {
            current->hemi_score = compute_score(nb, 1, 1, xu, current->aaw, current->llw, no_allele_frequency);
            current->het_score = compute_score(nb, 2, 1, xu, current->aaw, current->llw, no_allele_frequency);
            current->hom_score = compute_score(nb, 2, 2, xu, current->aaw, current->llw, no_allele_frequency);
        }
    }
}

void score_variant_t_b(struct variant * v, int nb, int xu, int no_allele_frequency){
    
    struct gene_transcript * c, * t;
    HASH_ITER(hh, v->gt, c, t) {
        struct transcript_anno_info * current, * tmp;
        HASH_ITER(hh, c->tai, current, tmp) {
            current->hemi_score = compute_score(nb, 1, 1, xu, current->aaw, current->llw, no_allele_frequency);
            current->het_score = compute_score(nb, 2, 1, xu, current->aaw, current->llw, no_allele_frequency);
            current->hom_score = compute_score(nb, 2, 2, xu, current->aaw, current->llw, no_allele_frequency);
        }
    }
}
