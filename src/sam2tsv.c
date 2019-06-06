/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
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

int cigar_clipleft(Cigar *cigar) {
	int *num = cigar->num;
	char *action = cigar->action;
	if (*action == 'H') {
		action++;
		num++;
	}
	if (*action == 'S') {
		return *num;
	} else {
		return 0;
	}
}

int cigar_clipright(Cigar *cigar) {
	int count = cigar->size;
	int *num = cigar->num+count-1;
	char *action = cigar->action+count-1;
	if (*action == 'H') {
		action--;
		num--;
	}
	if (*action == 'S') {
		return *num;
	} else {
		return 0;
	}
}

int cigar_qsize(Cigar *cigar) {
	int count = cigar->size;
	int *num = cigar->num;
	char *action = cigar->action;
	int result = 0;
	while(count--) {
		if (*action == 'M' || *action == 'I' || *action == 'S' || *action == '=' || *action == 'X') {
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

int main(int argc, char *argv[]) {
	FILE *f1=stdin;
	Cigar cigar;
	DStringArray *result1=NULL;
	DString *line1 = NULL,*linekeep = NULL;
	size_t read;
	unsigned int curpos=0, flag;
	unsigned int numfields,pos1;
	int start,end;
	char *cur;
	if ((argc != 1)) {
		fprintf(stderr,"Format is: sam2tsv\n");
		exit(EXIT_FAILURE);
	}
	/* allocate */
	cigar.memsize = 0;
	line1 = DStringNew(); linekeep=DStringNew();
	/* we add 2 to max because we need to have the column itself, and an extra space for the remainder */
	result1 = DStringArrayNew(10+2);
	/* skip sam header */
	read = DStringGetLine(line1, f1);
	while (read != -1) {
		if (line1->string[0] != '@') break;
		fprintf(stdout,"#%*.*s\n",line1->size,line1->size,line1->string);
		read = DStringGetLine(line1, f1); curpos++;
		pos1++;
	}
	/* init amplicons */
	DStringSplitTab(line1,10,result1,0,&numfields);
	fprintf(stdout,"chromosome\tbegin\tend\tstrand\tqname\tqstart\tqend\tmapquality\tref2\tbegin2\tstrand2\ttlen");
	fprintf(stdout,"\tpair\tproperpair\tunmapped\tmateunmapped\tread\tsecondary\tqcfail\tduplicate\tsupplementary");
	fprintf(stdout,"\tcigar\tseqlen\tseq\tquality\tother\n");
	while (1) {
		pos1++;
		sscanf(result1->data[1].string,"%d",&(flag));
		sscanf(result1->data[3].string,"%d",&start);
		parse_cigar(&cigar,result1->data[5].string);
		start--;
		end = start + cigar_refsize(&cigar);
		/* print result */
		DStringputs(result1->data+2,stdout); /* chromosome/refname */
		fprintf(stdout,"\t%d\t%d",start,end);
		fputc('\t',stdout);
		if (flag & BAM_FREVERSE) {fputc('-',stdout);} else {fputc('+',stdout);}
		fputc('\t',stdout);
		DStringputs(result1->data+0,stdout); /* qname */
		fputc('\t',stdout);
		fprintf(stdout,"%d",cigar_clipleft(&cigar)); /* qstart */
		fputc('\t',stdout);
		end = result1->data[9].size - cigar_clipright(&cigar);
		fprintf(stdout,"%d",end); /* qend */
		fputc('\t',stdout);
		DStringputs(result1->data+4,stdout); /* mapq */
		fputc('\t',stdout);
		DStringputs(result1->data+6,stdout); /* rnext */
		fputc('\t',stdout);
		sscanf(result1->data[7].string,"%d",&start);
		start--;
		fprintf(stdout,"%d",start); /* begin2 */
		fputc('\t',stdout);
		if (flag & BAM_FMREVERSE) {fputc('-',stdout);} else {fputc('+',stdout);}
		fputc('\t',stdout);
		DStringputs(result1->data+8,stdout); /* tlen */
		fputc('\t',stdout);
		if (flag & BAM_FPAIRED) {fputc('1',stdout);} else {fputc('0',stdout);}
		fputc('\t',stdout);
		if (flag & BAM_FPROPER_PAIR) {fputc('1',stdout);} else {fputc('0',stdout);}
		fputc('\t',stdout);
		if (flag & BAM_FUNMAP) {fputc('1',stdout);} else {fputc('0',stdout);}
		fputc('\t',stdout);
		if (flag & BAM_FMUNMAP) {fputc('1',stdout);} else {fputc('0',stdout);}
		fputc('\t',stdout);
		if (flag & BAM_FREAD1) {fputc('1',stdout);} else if (flag & BAM_FREAD2) {fputc('2',stdout);} else {fputc('0',stdout);}
		fputc('\t',stdout);
		if (flag & BAM_FSECONDARY) {fputc('1',stdout);} else {fputc('0',stdout);}
		fputc('\t',stdout);
		if (flag & BAM_FQCFAIL) {fputc('1',stdout);} else {fputc('0',stdout);}
		fputc('\t',stdout);
		if (flag & BAM_FDUP) {fputc('1',stdout);} else {fputc('0',stdout);}
		fputc('\t',stdout);
		if (flag & BAM_FSUPPL) {fputc('1',stdout);} else {fputc('0',stdout);}
		fputc('\t',stdout);
		DStringputs(result1->data+5,stdout); /* cigar */
		fputc('\t',stdout);
		fprintf(stdout,"%d",result1->data[9].size); /* seqlen */
		fputc('\t',stdout);
		DStringputs(result1->data+9,stdout); /* seq */
		fputc('\t',stdout);
		DStringputs(result1->data+10,stdout); /* qual */
		fputc('\t',stdout);
		cur=result1->data[11].string;
		while (*cur != '\0') {
			if (*cur == '\t') {fputc(' ',stdout);} else {fputc(*cur,stdout);}
			cur++;
		}
		fputc('\n',stdout);
		if (DStringGetTab(line1,f1,10,result1,0,&numfields)) break;
	}
	fclose(f1);
	if (line1) {DStringDestroy(line1);}
	if (linekeep) {DStringDestroy(linekeep);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
