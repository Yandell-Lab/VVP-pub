//
//  main.c
//  vvp_score
//
//  Created by steven on 8/11/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#include "vvp_headers.h"
#include "search_binary_bkgrnd.h"
#include "parse_vcf.h"
#include "vvp_lookup.h"
#include "score_variant.h"

#define WORK_SIZE 1000000

static sds input_vcf;
static sds db_prefix;
static sds output;
static int ncpus;
static int snv_only;
static int coding_only;
static int no_aa_weights;
static sds anno_tag_name;
static uint8_t gene_index;
static uint8_t transcript_index;
static uint8_t so_tag_index;
static uint8_t aa_index;
static int ll_weight_index;

static int n_background;
static unsigned char * mm_bin;
static struct chr_offsets * chro;

void usage(int exit_code) {
    fprintf(stderr, "Usage: VVP [options] -i <vcf file> -o <output prefix>\n\n");
    fprintf(stderr, "Options: (*mandatory)\n");
    fprintf(stderr, "* -i    filename      Input vcf file. Can be zipped or unzipped.  Can be 'stdin'\n");
    fprintf(stderr, "* -d    filename      database prefix\n");
    fprintf(stderr, "* -v    string        string with comma separated annotation components in info field\n");
    fprintf(stderr, "                      Format: <csq>,<gene index>,<transcript index>,<so_tag_index>,<aa_index>\n");
    fprintf(stderr, "                      Example: CSQ,4,6,1,15\n");
    fprintf(stderr, "-o      filename      fomatted output file name (for use in burden permutation)\n");
    fprintf(stderr, "-n      #             Number of threads to use, default = 1\n");
    fprintf(stderr, "-w      int           Column index (zero based) in annotation tag as extra likelihood weight\n");
    fprintf(stderr, "-x      None          Set to turn off AA scoring -- all AA weights will be set to 1.0\n");
    fprintf(stderr, "-l      flag          Set to ignore indels.  Default is to score indels\n");
    fprintf(stderr, "-c      flag          Set to ignore non-coding variants.  Default is to score non-coding variants.\n\n");
    exit(exit_code);
}

void parse_command_line(int argc, const char * argv[]) {
    int opt;
    int sig;
    sds * tmp_info;
    int tmp_count;
    if (argc > 1 && strcmp(argv[1], "-h") == 0)
        usage(0);
    while ((opt = getopt(argc, argv, "i:d:v:o:n:w:cxl")) != -1) {
        switch (opt) {
            case 'i' :
                input_vcf = sdsnew(optarg);
                break;
            case 'd' :
                db_prefix = sdsnew(optarg);
                break;
            case 'v' :
                tmp_info = sdssplitlen(optarg, (int)strlen(optarg), ",", 1, &tmp_count);
                if (tmp_count != 5) {
                    fprintf(stderr, "ARGUMENT ERROR:\tmust assign five annotation components, here only %d in %s\n", tmp_count, optarg);
                    usage(1);
                }
                anno_tag_name = sdsdup(tmp_info[0]);
                gene_index = atoi(tmp_info[1]);
                transcript_index = atoi(tmp_info[2]);
                so_tag_index = atoi(tmp_info[3]);
                aa_index = atoi(tmp_info[4]);
                sdsfreesplitres(tmp_info, tmp_count);
                break;
            case 'w':
                sig = atoi(optarg);
                if (sig < 0) {
                    fprintf(stderr, "ARGUMENT ERROR:\textra weight index must be >= 0\n");
                    usage(1);
                }
                ll_weight_index = sig;
                break;
            case 'l' :
                snv_only = 1;
                break;
            case 'x' :
                no_aa_weights = 1;
                break;
            case 'c' :
                coding_only = 1;
                break;
            case 'n' :
                sig = atoi(optarg);
                if (sig < 1){
                    fprintf(stderr, "ARGUMENT ERROR:\tnumber of cpus must be set to an integer > 0\n");
                    usage(1);
                }
                ncpus = (int)sig;
                break;
            case 'o' :
                output = sdsnew(optarg);
                break;
            default:
                usage(0);
                break;
        }
    }
    
    if (input_vcf == NULL || sdslen(input_vcf) < 2) {
        fprintf(stderr, "Missing mandatory option -i\n");
        usage(1);
    }
    if (db_prefix == NULL || sdslen(db_prefix) < 2) {
        fprintf(stderr, "Missing mandatory option -d\n");
        usage(1);
    }
    if (anno_tag_name == NULL || sdslen(anno_tag_name) < 2) {
        fprintf(stderr, "Missing mandatory option -v\n");
        usage(1);
    }
    
}

struct var_info * check_background_allele(struct m_var_info * bvs, struct variant * v){
    
    int i = 0;
    for (i = 0; i < bvs->nv; i++) {
        int diff = abs((int)sdslen(v->ref) - (int)sdslen(v->var));
        if (diff == 0) { //SNV or MNP
            
            if (sdslen(v->var) == 1) { //SNV
                if (bvs->vi[i]->var_type == v->var[0]) { //check to see if same SNV
                    return bvs->vi[i];
                }
            }
            
            else if (sdslen(v->var) > 1) { //MNP
                if (bvs->vi[i]->length == sdslen(v->var)) {
                    return bvs->vi[i];
                }
            }
        }
        else { //INDEL
            if (bvs->vi[i]->length == diff) {
                
                if ( ( (int)sdslen(v->ref) > (int)sdslen(v->var) ) && bvs->vi[i]->var_type == 'D') {
                    return bvs->vi[i]; //deletion
                }
                else if( ( (int)sdslen(v->ref) < (int)sdslen(v->var) ) && bvs->vi[i]->var_type == 'I'  ) {
                    return bvs->vi[i]; //insertion
                }
            }
        }
        
    }
    return NULL;
}

void id_variant_to_string(struct variant * v){
    
    size_t i;
    size_t n_indv = 0;
    if (v->hemi.n > n_indv) {
        n_indv = v->hemi.n;
    }
    if (v->hets.n > n_indv) {
        n_indv = v->hets.n;
    }
    if (v->homs.n > n_indv) {
        n_indv = v->homs.n;
    }
    
    
    
    for (i = 0; i < n_indv; i++) {
        if (i < v->hemi.n) {
            v->hemi_indv = sdscatprintf(v->hemi_indv, "%d,", kv_A(v->hemi, i));
        }
        if (i < v->hets.n) {
            v->het_indv = sdscatprintf(v->het_indv, "%d,", kv_A(v->hets, i));
        }
        if (i < v->homs.n) {
            v->hom_indv = sdscatprintf(v->hom_indv, "%d,", kv_A(v->homs, i));
        }
    }
    
    if(sdslen(v->hemi_indv) < 1){
        v->hemi_indv = sdscat(v->hemi_indv, ".");
    }
    if(sdslen(v->het_indv) < 1){
        v->het_indv = sdscat(v->het_indv, ".");
    }
    if(sdslen(v->hom_indv) < 1){
        v->hom_indv = sdscat(v->hom_indv, ".");
    }

}


struct variant * parse_score(sds vcf_line){
    
    struct variant * v = parse_vcf_line(vcf_line, no_aa_weights);
    struct m_var_info * bvs = search_binary_bkgrnd(v->chr, v->pos, mm_bin, chro);
    struct var_info * bvi = check_background_allele(bvs, v);
    if (bvi != NULL) {
        v->b_nhemi = bvi->nhemi;
        v->b_nhet = bvi->nhet;
        v->b_nhom = bvi->nhom;
        v->b_nocall = bvi->nocall;
        v->bit_offset = bvi->bit_offset;
        score_variant_t_b(v, n_background*2 - bvi->nocall, bvi->nhet + 2*bvi->nhom + bvi->nhemi);
    }
    else {
        v->b_nhemi = 0;
        v->b_nhet = 0;
        v->b_nhom = 0;
        v->b_nocall = 0;
        v->bit_offset = 0;
        score_variant_t_b(v, n_background*2, 0);
    }
    
    id_variant_to_string(v);
    
    struct gene_transcript * c, * t;
    HASH_ITER(hh, v->gt, c, t) {
        struct transcript_anno_info * current, * tmp;
        HASH_ITER(hh, c->tai, current, tmp) {
            
            current->hemi_vvp = score_lookup_b(current->transcript_name, current->hemi_score, current->coding);
            current->het_vvp = score_lookup_b(current->transcript_name, current->het_score, current->coding);
            current->hom_vvp = score_lookup_b(current->transcript_name, current->hom_score, current->coding);
        }
    }
    
    size_t i;
    for (i=0; i < bvs->nv; i++) {
        free(bvs->vi[i]);
    }
    free(bvs->vi);
    free(bvs);
    
    return v;
}

void process_vcf_lines(kvec_t(sds) * vcf_lines, struct variant *** variants){
    
    size_t n_lines = kv_size(*vcf_lines);
    size_t i;
    
    #pragma omp parallel for schedule(static)
    for (i = 0; i < n_lines; i++) {
        (*variants)[i] = parse_score(kv_A(*vcf_lines, i));
    }
    
}

/*
int scale_het(int x){
    float b = 0.055;
    return (int)100.0*(1.0 / (1.0 + exp(b*(10.0 - x))));
}*/


int main(int argc, const char ** argv) {
    
    input_vcf = NULL;
    db_prefix = NULL;
    anno_tag_name = NULL;
    gene_index = 0;
    transcript_index = 0;
    so_tag_index = 0;
    aa_index = 0;
    output = NULL;
    ncpus = 1;
    no_aa_weights = 0;
    snv_only = 0;
    coding_only = 0;
    n_background = 0;
    mm_bin = NULL;
    chro = NULL;
    ll_weight_index = -1;
    
    parse_command_line(argc, argv);
    
    #ifdef _OPENMP
    omp_set_num_threads(ncpus);
    #endif

    
    FILE * formatted_output = NULL;
    
    if (output != NULL) {
        formatted_output = fopen(output, "w");
    }
    
    sds dist_output = sdsempty();
    dist_output = sdscatprintf(dist_output, "%s.dist", db_prefix);
    load_feature_lookups_b(dist_output);
    sdsfree(dist_output);
    
    mm_bin = load_bin_db(db_prefix, &n_background); //create memory map of background
    chro = load_offsets(db_prefix); //load byte offsets in memory map
    
    initialize_parse_vcf(gene_index, transcript_index, so_tag_index, aa_index, anno_tag_name, ll_weight_index);
    
    gzFile * gf = NULL;
    
    if (strcmp(input_vcf, "stdin") == 0){
        gf = gzdopen(STDIN_FILENO, "r");
    }
    else {
        gf = gzopen(input_vcf, "r");
    }
    
    if (! gf) {
        fprintf(stderr, "FATAL: vcf file %s cannot be read\n", input_vcf);
        exit(1);
    }
    
    fprintf(stdout, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","chr", "start", "ref", "var", "gene", "transcript", "hemi_score", "hemi_vvp", "nhemi", "hemi_indvs", "hemi_nocall", "het_score", "het_vvp", "nhet", "het_indvs", "het_nocall", "hom_score", "hom_vvp", "nhom", "hom_indvs", "hom_nocall", "coding_ind", "indel_ind", "aa_score", "n_bhemi", "n_bhet", "n_bhom", "n_bnocall", "bit_offset", "vid", "ll_weight");
    
    char * buffer = (char *)malloc(sizeof(char)*BUF_SIZE);
    kvec_t(sds) vcf_lines;
    kv_init(vcf_lines);
    char * l = gzgets(gf, buffer, BUF_SIZE);
    if (buffer[0] != '#') {
        kv_push(sds, vcf_lines, sdsnew(buffer));
    }
    
    while (l != NULL) {
        
        while (kv_size(vcf_lines) < WORK_SIZE) {
            l = gzgets(gf, buffer, BUF_SIZE);
            if (l == NULL) {
                break;
            }
            if (buffer[0] != '#') {
                kv_push(sds, vcf_lines, sdsnew(buffer));
            }
        }
        
        
        if (kv_size(vcf_lines) < 1) {
            break;
        }
        
        size_t i;
        size_t n_lines = kv_size(vcf_lines);
        struct variant ** variants = (struct variant **)calloc(n_lines, sizeof(struct variant *));
        process_vcf_lines(&vcf_lines, &variants);
        
        for (i = 0; i < n_lines; i++) {
            struct variant * tv = variants[i];
            if (tv != NULL) {
                int indel_ind = (int)sdslen(tv->var) - (int)sdslen(tv->ref) == 0 ? 0 : 1;
                if ((tv->hemi.n > 0 || tv->hets.n > 0 || tv->homs.n > 0) && (snv_only == 0 || snv_only > indel_ind)) {
                    struct gene_transcript * c, * t;
                    HASH_ITER(hh, tv->gt, c, t) {
                        struct transcript_anno_info * current, * tmp;
                        HASH_ITER(hh, c->tai, current, tmp) {
                            if (coding_only <= current->coding) {
                                fprintf(stdout, "%s\t%zu\t%s\t%s\t%s\t%s\t", tv->chr, tv->pos, tv->ref, tv->var, c->gene_name, current->transcript_name);
                                fprintf(stdout, "%f\t%d\t%zu\t%s\t%zu\t", current->hemi_score, current->hemi_vvp, tv->hemi.n, tv->hemi_indv, tv->hemi_nocalls.n);
                                fprintf(stdout, "%f\t%d\t%zu\t%s\t%zu\t", current->het_score, current->het_vvp, tv->hets.n, tv->het_indv, tv->het_nocalls.n);
                                fprintf(stdout, "%f\t%d\t%zu\t%s\t%zu\t", current->hom_score, current->hom_vvp, tv->homs.n, tv->hom_indv, tv->hom_nocalls.n);
                                fprintf(stdout, "%d\t%d\t%f\t", current->coding, indel_ind, current->aaw);
                                fprintf(stdout, "%d\t%d\t%d\t%d\t", tv->b_nhemi, tv->b_nhet, tv->b_nhom, tv->b_nocall);
                                fprintf(stdout, "%llu\t%s\t%f\n", tv->bit_offset, tv->vid, current->llw);
                                
                                if (output != NULL) {
                                    int hemi_ind = (tv->b_nhemi > 0) ? 1 : 0;
                                    int het_ind = (tv->b_nhet > 0) ? 1 : 0;
                                    int hom_ind = (tv->b_nhom > 0) ? 1 : 0;
                                    fprintf(formatted_output, "%s\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%d\t%d\t%d\t%llu\t%s\n", tv->chr, current->transcript_name, current->hemi_score, tv->hemi_indv, current->het_score, tv->het_indv, current->hom_score, tv->hom_indv, hemi_ind, het_ind, hom_ind, tv->bit_offset, tv->vid);
                                }
                                
                            }
                        }
                    }
                }
                
                destroy_variant(tv);
            }
            
            sdsfree(kv_pop(vcf_lines));
        
        }
        
    }
        
    kv_destroy(vcf_lines);
    free(buffer);
    
    if (output != NULL) {
        fprintf(stderr, "\nsorting and prepping formatted output for burden calculations (stdout ready for processing)...");
        //prep formatted output file for burden calculations; only works on unix systems
        sds sort_command = sdsnew("sort -k2,2 ");
        sort_command = sdscatprintf(sort_command, "%s > %s.sorted", output, output);
        system(sort_command);
        sdsfree(sort_command);
        
        sds mv_command = sdsnew("mv ");
        mv_command = sdscatprintf(mv_command, "%s.sorted %s", output, output);
        system(mv_command);
        fprintf(stderr, "done\n");
        sdsfree(mv_command);
    }
    
    
    
    
    return 0;
}
