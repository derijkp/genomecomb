/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tools.h"
#include "gztools.h"
#include "debug.h"

inline int min ( int a, int b ) { return a < b ? a : b; }

inline int max ( int a, int b ) { return a > b ? a : b; }

#define AMPLICONSBUFFERSIZE 4

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract suplementary alignment (part of chimeric alignment) */
#define BAM_FSUPPL      2048

typedef struct Amplicon {
	DString *chr2;
	int start2;
	int end2;
	DString *name;
} Amplicon;

typedef struct Cigar {
	int size;
	int memsize;
	int *num;
	char *action;
} Cigar;

int cigar_refsize(Cigar *cigar) {
	int count = cigar->size;
	int *num = cigar->num;
	char *action = cigar->action;
	int result = 0;
	while(count--) {
		if (*action == 'M' || *action == 'D' || *action == 'N' || *action == '=' || *action == 'X') {
			result += *num;
		}
		action++;
		num++;
	}
	return result;
}

int parse_cigar(Cigar *cigar,char *string) {
	if (cigar->memsize == 0) {
		cigar->size = 0;
		cigar->memsize = 5;
		cigar->num = malloc(cigar->memsize*sizeof(int));
		cigar->action = malloc(cigar->memsize*sizeof(char));
	}
	cigar->size = 0;
	while (*string && *string != '\t') {
		if (cigar->size >= cigar->memsize) {
			cigar->memsize += 5;
			cigar->num = realloc(cigar->num,cigar->memsize*sizeof(int));
			cigar->action = realloc(cigar->action,cigar->memsize*sizeof(char));
		}
		cigar->num[cigar->size] = (int)strtol(string,&string,10);
		cigar->action[cigar->size++] = *string++;
	}
	return(cigar->size);
}

#define POSQNAME 0
#define POSREFNAME 1
#define POSBEGIN 2
#define POSEND 3
#define POSSTRAND 4
#define POSMAPQUALITY 5
#define POSREF2 6
#define POSBEGIN2 7
#define POSSTRAND2 8
#define POSTLEN 9
#define POSPAIR 10
#define POSPROPERPAIR 11
#define POSUNMAPPED 12
#define POSMATEUNMAPPED 13
#define POSREAD 14
#define POSSECONDARY 15
#define POSQCFAIL 16
#define POSDUPLICATE 17
#define POSSUPPLEMENTARY 18
#define POSCIGAR 19
#define POSSEQLEN 20
#define POSSEQ 21
#define POSQUALITY 22
#define POSOTHER 23

int main(int argc, char *argv[]) {
	FILE *f1=stdin;
	DStringArray *result1=NULL, *header = NULL;
	DString *line1 = NULL,*linekeep = NULL, *headerline = NULL;
	ssize_t read;
	unsigned int curpos=0, flag, header_numfields;
	unsigned int numfields,pos1, count, max, i;
	int pos;
	char *fields[] = {"qname","refname","begin","end","strand","mapquality","ref2","begin2","strand2","tlen","pair","properpair","unmapped","mateunmapped","read","secondary","qcfail","duplicate","supplementary","cigar","seqlen","seq","quality","other"};
	int poss[24];
	char *cur;
	if ((argc != 1)) {
		fprintf(stderr,"Format is: sam2tsv\n");
		exit(EXIT_FAILURE);
	}
	/* allocate */
	line1 = DStringNew(); linekeep=DStringNew();
	/* skip sam header */
	read = DStringGetLine(line1, f1);
	while (read != -1) {
		if (line1->string[0] != '#') break;
		if (line1->size > 2 && line1->string[1] == '@') {
			fprintf(stdout,"%*.*s\n",line1->size-1,line1->size-1,line1->string+1);
		}
		read = DStringGetLine(line1, f1); curpos++;
		pos1++;
	}
	headerline = line1; line1 = DStringNew();
	header_numfields = 1;
	count = headerline->size; cur = headerline->string;
	while (count--) {if (*cur++ == '\t') header_numfields++;}
	header = DStringArrayNew(header_numfields+2);
	DStringSplitTab(headerline, header_numfields, header, 0,NULL);
	max = 0;
	for( i = 0 ; i < 24 ; i++ ) {
		poss[i] = DStringArraySearch(header,fields[i],strlen(fields[i]));
		if (poss[i] == -1) {
			fprintf(stderr,"field \"%s\" not found",fields[i]);
			exit(1);
		}
		if (poss[i] > max) max = poss[i];
	}
	/* we add 2 to max because we need to have the column itself, and an extra space for the remainder */
	result1 = DStringArrayNew(max+2);
	while (1) {
		if (DStringGetTab(line1,f1,max,result1,0,&numfields)) break;
		/* make flag */
		flag = 0;
		if (result1->data[poss[POSPAIR]].string[0] == '1') {flag |= BAM_FPAIRED;}
		if (result1->data[poss[POSPROPERPAIR]].string[0] == '1') {flag |= BAM_FPROPER_PAIR;}
		if (result1->data[poss[POSUNMAPPED]].string[0] == '1') {flag |= BAM_FUNMAP;}
		if (result1->data[poss[POSMATEUNMAPPED]].string[0] == '1') {flag |= BAM_FMUNMAP;}
		if (result1->data[poss[POSREAD]].string[0] == '1') {
			flag |= BAM_FREAD1;
		} else if (result1->data[poss[POSREAD]].string[0] == '2') {
			flag |= BAM_FREAD2;
		}
		if (result1->data[poss[POSSECONDARY]].string[0] == '1') {flag |= BAM_FSECONDARY;}
		if (result1->data[poss[POSQCFAIL]].string[0] == '1') {flag |= BAM_FQCFAIL;}
		if (result1->data[poss[POSDUPLICATE]].string[0] == '1') {flag |= BAM_FDUP;}
		if (result1->data[poss[POSSUPPLEMENTARY]].string[0] == '1') {flag |= BAM_FSUPPL;}

		if (result1->data[poss[POSSTRAND]].string[0] == '-') {flag |= BAM_FREVERSE;}
		if (result1->data[poss[POSSTRAND2]].string[0] == '-') {flag |= BAM_FMREVERSE;}
		/* print */
		DStringputs(result1->data+poss[POSQNAME],stdout); /* qname */
		fputc('\t',stdout);
		fprintf(stdout,"%d",flag); /* flag */
		fputc('\t',stdout);
		DStringputs(result1->data+poss[POSREFNAME],stdout); /* rname */
		fputc('\t',stdout);
		sscanf(result1->data[poss[POSBEGIN]].string,"%d",&pos);
		fprintf(stdout,"%d",pos+1); /* pos */
		fputc('\t',stdout);
		DStringputs(result1->data+poss[POSMAPQUALITY],stdout); /* mapq */
		fputc('\t',stdout);
		DStringputs(result1->data+poss[POSCIGAR],stdout); /* cigar */
		fputc('\t',stdout);
		DStringputs(result1->data+poss[POSREF2],stdout); /* rnext */
		fputc('\t',stdout);
		sscanf(result1->data[poss[POSBEGIN2]].string,"%d",&pos);
		fprintf(stdout,"%d",pos+1); /* pnext */
		fputc('\t',stdout);
		DStringputs(result1->data+poss[POSTLEN],stdout); /* tlen */
		fputc('\t',stdout);
		DStringputs(result1->data+poss[POSSEQ],stdout); /* seq */
		fputc('\t',stdout);
		DStringputs(result1->data+poss[POSQUALITY],stdout); /* qual */
		count=result1->data[poss[POSOTHER]].size;
		if (count) {
			fputc('\t',stdout);
			cur=result1->data[poss[POSOTHER]].string;
			while (count-- && *cur != '\0') {
				if (*cur == ' ') {fputc('\t',stdout);} else {fputc(*cur,stdout);}
				cur++;
			}
		}
		fputc('\n',stdout);
	}
	fclose(f1);
	if (line1) {DStringDestroy(line1);}
	if (linekeep) {DStringDestroy(linekeep);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
