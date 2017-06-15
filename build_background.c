//
//  main.c
//  vcf_parser
//
//  Created by steven on 5/20/15.
//  Copyright (c) 2015 yandell lab. All rights reserved.
//

#include "vvp_headers.h"
#include "parse_vcf.h"
#include "score_variant.h"
#include "bit_array.h"

#define BUF_SIZE 5000000
#define WORK_SIZE 1000000

struct chromosome_info{
    char chromosome[3];
    uint64_t start;
    uint64_t end;
    int n;
    UT_hash_handle hh;
};


struct feature_scores{
    char feature_name[FEATURE_NAME_LENGTH];
    sds gene_name;
    kvec_t(float) coding;
    kvec_t(float) noncoding;
    kvec_t(float) max_hom;
    kvec_t(float) max_het1;
    kvec_t(float) max_het2;
    kvec_t(float) sum_scores;
    UT_hash_handle hh;
};


static sds input_vcf;
static sds output_prefix;
static int ncpus;
static uint8_t indel_only;
static uint8_t coding_only;
static int no_aa_weights;
static int n_background;
static uint64_t byte_offset;
static uint64_t bit_offset;
static sds anno_tag_name;
static uint8_t gene_index;
static uint8_t transcript_index;
static uint8_t so_tag_index;
static uint8_t aa_index;
static int ll_weight_index;

void usage(int exit_code) {
    fprintf(stderr, "Usage: build_background [options] -i <vcf file> -o <output prefix>\n\n");
    fprintf(stderr, "Options: (*mandatory)\n");
    fprintf(stderr, "* -i    filename      Input vcf file. Can be zipped or unzipped\n");
    fprintf(stderr, "* -o    filename      Output prefix for database build\n");
    fprintf(stderr, "* -b    #             Number of individuals in vcf file (Must be EXACT)\n");
    fprintf(stderr, "* -v    string        string with comma separated annotation components in info field\n");
    fprintf(stderr, "                      Format: <csq>,<gene index>,<transcript index>,<so_tag_index>,<aa_index>\n");
    fprintf(stderr, "                      Example: CSQ,4,6,1,15\n");
    fprintf(stderr, "-w      int           Column index (zero based) in annotation tag as extra likelihood weight\n");
    fprintf(stderr, "-n      #             Number of threads to use while parsing, default = 1\n");
    fprintf(stderr, "-x      None          Set to turn off AA scoring -- all AA weights will be set to 1.0\n");
    fprintf(stderr, "-d      None          Set to ignore indels.  Default is to use indels\n");
    fprintf(stderr, "-c      None          Set to ignore non-coding variants.  Default is to use non-coding variants.\n\n");
    exit(exit_code);
}

void parse_command_line(int argc, const char * argv[]) {
    int opt;
    int sig;
    sds * tmp_info;
    int tmp_count;
    if (argc > 1 && strcmp(argv[1], "-h") == 0)
        usage(0);
    while ((opt = getopt(argc, argv, "i:o:b:v:n:w:cxd")) != -1) {
        switch (opt) {
            case 'i' :
                input_vcf = sdsnew(optarg);
                break;
            case 'd' :
                indel_only = 1;
                break;
            case 'c' :
                coding_only = 1;
                break;
            case 'x' :
                no_aa_weights = 1;
                break;
            case 'w':
                sig = atoi(optarg);
                if (sig < 0) {
                    fprintf(stderr, "ARGUMENT ERROR:\textra weight index must be >= 0\n");
                    usage(1);
                }
                ll_weight_index = sig;
                break;
            case 'b' :
                sig = atoi(optarg);
                if (sig < 1) {
                    fprintf(stderr, "ARGUMENT ERROR:\tnumber of individuals in vcf must be > 0\n");
                    usage(1);
                }
                n_background = sig;
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
            case 'n' :
                sig = atoi(optarg);
                if (sig < 1){
                    fprintf(stderr, "ARGUMENT ERROR:\tnumber of cpus must be set to an integer > 0\n");
                    usage(1);
                }
                ncpus = (int)sig;
                break;
            case 'o' :
                output_prefix = sdsnew(optarg);
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
    if (output_prefix == NULL || sdslen(output_prefix) < 2) {
        fprintf(stderr, "Missing mandatory option -o\n");
        usage(1);
    }
    if (n_background < 1) {
        fprintf(stderr, "Missing mandatory option -b\n");
        usage(1);
    }
    if (anno_tag_name == NULL || sdslen(anno_tag_name) < 2) {
        fprintf(stderr, "Missing mandatory option -v\n");
        usage(1);
    }
}

void write_chromosome_offsets(char * output_pre, struct chromosome_info * chr_info){
    
    sds output = sdsempty();
    output = sdscatprintf(output, "%s.chr_offsets.txt", output_pre);
    FILE * chro = fopen(output, "w");
    sdsfree(output);
    if (chro == NULL) {
        fprintf(stderr, "FATAL: could not open file for chromosome offsets: %s", output);
        exit(1);
    }
    
    struct chromosome_info * tci, *tmp;
    HASH_ITER(hh, chr_info, tci, tmp){
        fprintf(chro, "%s\t%llu\t%llu\t%d\n",tci->chromosome,tci->start,tci->end,tci->n);
        HASH_DEL(chr_info, tci);
        free(tci);
    }
    fclose(chro);
    
}

struct variant * parse_score(sds vcf_line){
    
    struct variant * v = parse_vcf_line(vcf_line, no_aa_weights);
    
    score_variant_b(v);
    
    return v;
}

void write_binary(FILE * bin_output, FILE * bit_out, struct variant * v, struct chromosome_info ** chr_info){
    
    struct chromosome_info * tci = NULL;
    HASH_FIND_STR(*chr_info, v->chr, tci);
    if (tci == NULL) {
        struct chromosome_info * ci = (struct chromosome_info *)malloc(sizeof(struct chromosome_info));
        memset(ci->chromosome, '\0', sizeof(char)*3);
        strcpy(ci->chromosome, v->chr);
        ci->start = byte_offset;
        ci->end = byte_offset + LINE_BYTE_SIZE;
        ci->n = 1;
        HASH_ADD_STR(*chr_info, chromosome, ci);
    }
    else {
        tci->end += LINE_BYTE_SIZE;
        tci->n += 1;
    }
    
    sds var_type = sdsnew(v->var); //default is SNV
    
    int length = abs((int)sdslen(v->ref) - (int)sdslen(v->var));
    
    if (length < 1 && sdslen(v->var) > 1) { //case of MNP
        var_type = sdscpy(var_type, "M");
        length = (int)sdslen(v->var);
    }
    else if (length > 0){ //indel
        var_type = sdscpy(var_type, "I");
        if (sdslen(v->var) < sdslen(v->ref)) {
            var_type = sdscpy(var_type, "D");
        }
    }
    
    if (sdslen(v->chr) < 2) {
        fwrite("0", sizeof(char), 1, bin_output);
        fwrite(v->chr, sizeof(char), 1, bin_output);
    }
    else{
        fwrite(v->chr, sizeof(char), 2, bin_output);
    }
    
    fwrite(&(v->pos), sizeof(int), 1, bin_output);
    fwrite(var_type, sizeof(char), 1, bin_output);
    fwrite(&length, sizeof(int), 1, bin_output);
    int nhet = (int)kv_size(v->hets);
    fwrite(&nhet, sizeof(int), 1, bin_output);
    int nhom = (int)kv_size(v->homs);
    fwrite(&nhom, sizeof(int), 1, bin_output);
    int nhemi = (int)kv_size(v->hemi);
    fwrite(&nhemi, sizeof(int), 1, bin_output);
    int nnocall = (int)(kv_size(v->het_nocalls) + kv_size(v->hemi_nocalls) + 2*kv_size(v->hom_nocalls));
    fwrite(&nnocall, sizeof(int), 1, bin_output);
    fwrite(&bit_offset, sizeof(uint64_t), 1, bin_output);
    
    char bits[n_background+1];
    memset(bits, '\0', sizeof(char)*(n_background+1));
    int i = 0;
    for (i = 0; i < kv_size(v->hets); i++) {
        bits[kv_A(v->hets, i)] = '1';
    }
    for (i = 0; i < kv_size(v->homs); i++) {
        bits[kv_A(v->homs, i)] = '1';
    }
    for (i = 0; i < kv_size(v->hemi); i++) {
        bits[kv_A(v->hemi, i)] = '1';
    }
    BIT_ARRAY * ba = bit_array_create(n_background);
    bit_array_from_str(ba, bits);
    bit_offset += bit_array_save(ba, bit_out);
    bit_array_free(ba);
    
    //fprintf(txt_out, "%s\t%d\t%s\t%d\t%d\t%d\t%d\n", v->chr, v->start, var_type, length, v->nhet, v->nhom, v->nocall);
    
    byte_offset += LINE_BYTE_SIZE;
    
    
    sdsfree(var_type);
    
}

void populate_feature_hash(struct feature_scores ** fs, struct variant * v){
    
    struct feature_scores * tfs = NULL;
    
    //populate feature scores hash
    struct gene_transcript * c, * t;
    HASH_ITER(hh, v->gt, c, t) {
        struct transcript_anno_info * current, * tmp;
        HASH_ITER(hh, c->tai, current, tmp) {
            HASH_FIND_STR(*fs, current->transcript_name, tfs);
            if (tfs == NULL) {
                tfs = (struct feature_scores *)malloc(sizeof(struct feature_scores));
                strcpy(tfs->feature_name, current->transcript_name);
                tfs->gene_name = sdsnew(v->gt->gene_name);
                kv_init(tfs->max_het1);
                kv_init(tfs->max_het2);
                kv_init(tfs->max_hom);
                kv_init(tfs->sum_scores);
                kv_init(tfs->coding);
                kv_init(tfs->noncoding);
                int i;
                for (i = 0; i < n_background; i++) {
                    kv_push(float, tfs->max_het1, 0.0);
                    kv_push(float, tfs->max_het2, 0.0);
                    kv_push(float, tfs->max_hom, 0.0);
                    kv_push(float, tfs->sum_scores, 0.0);
                }
                
                HASH_ADD_STR(*fs, feature_name, tfs);
            }
            //print scores to stdout
            fprintf(stdout, "%s\t%zu\t%s\t%s\t%s\t%f\t%f\t%f\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%d\n", v->chr, v->pos, v->ref, v->var, current->transcript_name, current->hemi_score, current->het_score, current->hom_score, v->hemi.n, v->hets.n, v->homs.n, v->hemi_nocalls.n, v->het_nocalls.n, v->hom_nocalls.n, current->coding);
            
            //populate variants into coding / noncoding categories
            if (current->coding == 1) {
                if (current->hemi_score >= 0.0) {
                    kv_push(float, tfs->coding, current->hemi_score);
                }
                if (current->het_score >= 0.0) {
                    kv_push(float, tfs->coding, current->het_score);
                }
                if (current->hom_score >= 0.0) {
                    kv_push(float, tfs->coding, current->hom_score);
                }
            }
            else {
                if (current->hemi_score >= 0.0) {
                    kv_push(float, tfs->noncoding, current->hemi_score);
                }
                if (current->het_score >= 0.0) {
                    kv_push(float, tfs->noncoding, current->het_score);
                }
                if (current->hom_score >= 0.0) {
                    kv_push(float, tfs->noncoding, current->hom_score);
                }
            }
            
            //populate max variant score
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
                    if (current->hemi_score > kv_A(tfs->max_het1, kv_A(v->hemi, i))) {
                        kv_A(tfs->max_het2, kv_A(v->hemi, i)) = kv_A(tfs->max_het1, kv_A(v->hemi, i));
                        kv_A(tfs->max_het1, kv_A(v->hemi, i)) = current->hemi_score;
                    }
                    else if (current->hemi_score > kv_A(tfs->max_het2, kv_A(v->hemi, i))){
                        kv_A(tfs->max_het2, kv_A(v->hemi, i)) = current->hemi_score;
                    }
                }
                if (i < v->hets.n) {
                    if (current->het_score > kv_A(tfs->max_het1, kv_A(v->hets, i))) {
                        kv_A(tfs->max_het2, kv_A(v->hets, i)) = kv_A(tfs->max_het1, kv_A(v->hets, i));
                        kv_A(tfs->max_het1, kv_A(v->hets, i)) = current->het_score;
                    }
                    else if (current->het_score > kv_A(tfs->max_het2, kv_A(v->hets, i))){
                        kv_A(tfs->max_het2, kv_A(v->hets, i)) = current->het_score;
                    }
                }
                if (i < v->homs.n) {
                    if (current->hom_score > kv_A(tfs->max_hom, kv_A(v->homs, i))) {
                        kv_A(tfs->max_hom, kv_A(v->homs, i)) = current->hom_score;
                    }
                }
            }
            
        }
    }
    
}

int float_cmp(const void *a, const void *b){
    float ia = *(const float *)a;
    float ib = *(const float *)b;
    return (ia > ib) - (ia < ib);
}

void calc_percentiles(float * scores, size_t nscores, float ** percentiles){
    
    qsort(scores, nscores, sizeof(float), float_cmp);
    
    int i=0;
    int si = 1;
    while (i < NPERCENTILES && (si-1) < nscores) {
        float tsi = (float)si / (float)nscores;
        float tp = (float)(i+1) / (float)NPERCENTILES;
        if ( tsi >= tp) {
            (*percentiles)[i] = scores[si-1];
            i++;
        }
        else {
            si++;
        }
    }
    (*percentiles)[NPERCENTILES - 1] = scores[nscores - 1]; //max needs to be highest percentile
    
}


void output_percentiles_b(FILE * out, struct feature_scores * fs){
    
    float tpercentiles[NPERCENTILES]; //percentiles are 1-100 inclusively
    float * percentiles = tpercentiles;
    
    fwrite(fs->feature_name, sizeof(char), FEATURE_NAME_LENGTH, out);
    
    memset(percentiles, '\0', sizeof(float)*NPERCENTILES);
    if (fs->coding.n > 0) {
        calc_percentiles(fs->coding.a, fs->coding.n, &percentiles);
    }
    
    fwrite(percentiles, sizeof(float), NPERCENTILES, out); //coding percentiles
    
    
    memset(percentiles, '\0', sizeof(float)*NPERCENTILES);
    if (fs->noncoding.n > 0) {
        calc_percentiles(fs->noncoding.a, fs->noncoding.n, &percentiles);
    }
    
    fwrite(percentiles, sizeof(float), NPERCENTILES, out); //noncoding percentiles
    
    fwrite(&(fs->coding.n), sizeof(size_t), 1, out); //n coding
    
    fwrite(&(fs->noncoding.n), sizeof(size_t), 1, out);//n noncoding
    
}

void output_individual_max_b(FILE * imax_scores, struct feature_scores * s, int nb){ //write binary
    
    int i=0;
    fwrite(s->feature_name, sizeof(char), FEATURE_NAME_LENGTH, imax_scores);
    for (i=0; i < nb; i++) {
        fwrite(&(kv_A(s->max_hom, i)), sizeof(float), 1, imax_scores);
        fwrite(&(kv_A(s->max_het1, i)), sizeof(float), 1, imax_scores);
        fwrite(&(kv_A(s->max_het2, i)), sizeof(float), 1, imax_scores);
        fwrite(&(kv_A(s->sum_scores, i)), sizeof(float), 1, imax_scores);
    }
}

int feature_sort(struct feature_scores * a, struct feature_scores * b) {
    return strcmp(a->feature_name,b->feature_name);
}

void process_vcf_lines(kvec_t(sds) * vcf_lines, struct variant *** variants){
    
    size_t n_lines = kv_size(*vcf_lines);
    size_t i;
    
    #pragma omp parallel for schedule(static)
    for (i = 0; i < n_lines; i++) {
        (*variants)[i] = parse_score(kv_A(*vcf_lines, i));
    }
}

int main(int argc, const char * argv[]) {
    input_vcf = NULL;
    output_prefix = NULL;
    n_background = 0;
    ncpus = 1;
    no_aa_weights = 0;
    indel_only = 0;
    coding_only = 0;
    byte_offset = 0;
    bit_offset = 0;
    anno_tag_name = NULL;
    gene_index = 0;
    transcript_index = 0;
    so_tag_index = 0;
    aa_index = 0;
    ll_weight_index = -1;
    
    size_t lines_processed = 0;
    
    struct chromosome_info * chr_info = NULL; //for binary chromosome info
    struct feature_scores * fs = NULL;
    
    parse_command_line(argc, argv);
    
    #ifdef _OPENMP
    omp_set_num_threads(ncpus);
    #endif
    
    sds binary_background = sdsnew(output_prefix);
    binary_background = sdscat(binary_background, ".bin");
    FILE * binary_out = fopen(binary_background, "wb");
    fwrite(&n_background, sizeof(int), 1, binary_out);
    byte_offset += sizeof(int);
    
    sds bit_background = sdsnew(output_prefix);
    bit_background = sdscat(bit_background, ".bit");
    FILE * bit_out = fopen(bit_background, "wb");
    
    initialize_parse_vcf(gene_index, transcript_index, so_tag_index, aa_index, anno_tag_name, ll_weight_index);
    
    gzFile * gf = gzopen(input_vcf, "r");
    if (! gf) {
        fprintf(stderr, "FATAL: vcf file %s cannot be read\n", input_vcf);
        exit(1);
    }
    
    fprintf(stdout, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","chr", "start", "ref", "var", "transcript", "hemi_score", "het_score", "hom_score", "nhemi", "nhet", "nhom", "hemi_nocall", "het_nocall", "hom_nocall", "coding_ind");
    
    char * buffer = (char *)malloc(sizeof(char)*BUF_SIZE);
    kvec_t(sds) vcf_lines;
    kv_init(vcf_lines);
    char * l = gzgets(gf, buffer, BUF_SIZE);
    
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
        
        size_t n_lines = kv_size(vcf_lines);
        struct variant ** variants = (struct variant **)calloc(n_lines, sizeof(struct variant *));
        process_vcf_lines(&vcf_lines, &variants);
    
        lines_processed += n_lines;
        fprintf(stderr, "# lines processed: %zu\n", lines_processed);
        size_t i;
        for (i = 0; i < n_lines; i++) {
            struct variant * tv = variants[i];
            if (tv != NULL) {
                write_binary(binary_out, bit_out, tv, &chr_info);
                populate_feature_hash(&fs, tv);
                destroy_variant(tv);
                sdsfree(kv_pop(vcf_lines));
            }
            else {
                fprintf(stderr, "WARNING: variant struct NULL, at pos %zu of %zu\n",i,n_lines);
            }
            
        }
        free(variants);
    }
    
        
    kv_destroy(vcf_lines);
    free(buffer);
    
    fprintf(stderr, "DONE READING VCF file, now will populate chr offsets, dist, max files\n");
    
    write_chromosome_offsets(output_prefix, chr_info);
    
    sds dist_output = sdsempty();
    dist_output = sdscatprintf(dist_output, "%s.dist", output_prefix);
    FILE * bkgrnd_dist = fopen(dist_output, "wb");
    sdsfree(dist_output);
    if (bkgrnd_dist == NULL) {
        fprintf(stderr, "FATAL: could not open distribution output file for writing %s", output_prefix);
        exit(1);
    }
    
    sds max_output = sdsempty();
    max_output = sdscatprintf(max_output, "%s.max", output_prefix);
    FILE * imax_scores = fopen(max_output, "wb");
    sdsfree(max_output);
    if (imax_scores == NULL) {
        fprintf(stderr, "FATAL: could not open max output file for writing %s\n", output_prefix);
        exit(1);
    }
    
    HASH_SORT(fs, feature_sort);
    struct feature_scores *s, *tmp;
    HASH_ITER(hh, fs, s, tmp){
        output_percentiles_b(bkgrnd_dist, s);
        output_individual_max_b(imax_scores, s, n_background);
        HASH_DEL(fs, s);
        sdsfree(s->gene_name);
        kv_destroy(s->max_het1);
        kv_destroy(s->max_het2);
        kv_destroy(s->max_hom);
        kv_destroy(s->sum_scores);
        kv_destroy(s->coding);
        kv_destroy(s->noncoding);
        free(s);
    }

    
    return 0;
}

