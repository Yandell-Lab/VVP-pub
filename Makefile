CC = gcc
CFLAGS = -lz -lm -O3 -lgsl -lgslcblas -fopenmp #-Wall
TARGETS = build_background VVP

all: $(TARGETS)

build_background: aa_weight.o bit_array.o score_variant.o parse_vcf.o sds.o build_background.o
	$(CC) $(CFLAGS) aa_weight.o bit_array.o score_variant.o parse_vcf.o sds.o build_background.o -o $@

VVP: aa_weight.o bit_array.o score_variant.o parse_vcf.o sds.o vvp_lookup.o search_binary_bkgrnd.o score_variants.o
	$(CC) $(CFLAGS) aa_weight.o bit_array.o score_variant.o parse_vcf.o sds.o vvp_lookup.o search_binary_bkgrnd.o score_variants.o -o $@ 

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o
	rm -f $(TARGETS)

