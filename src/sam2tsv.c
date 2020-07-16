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
#include "hash.h"
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
	Hash_table *hashtable;
	FILE *f1=stdin;
	Cigar cigar;
	DStringArray *result1=NULL;
	DStringArray *fieldsa=NULL;
	DString *line1 = NULL,*linekeep = NULL, *header = NULL;
	size_t read;
	unsigned int curpos=0, flag;
	unsigned int numfields,pos1;
	int start,end,fpos,contig=0;
	char *fields = NULL, *cur,*keep;
	if (argc < 1) {
		fprintf(stderr,"Format is: sam2tsv ?field? ...\n");
		exit(EXIT_FAILURE);
	}
	/* fields data */
	fprintf(stdout,"#filetype tsv/samfile\n#fileversion	0.99\n#fields	table\n#fields	field	number	type	description\n#fields	chromosome	1	String	Chromosome/Contig\n#fields	begin	1	Integer	Begin of feature (0 based - half open)\n#fields	end	1	Integer	End of feature (0 based - half open)\n#fields	type	1	String	Type of feature (snp,del,ins,...)\n#fields	ref	1	String	Reference sequence, can be a number for large features\n#fields	alt	1	String	Alternative sequence, can be a number for large features\n#fields	name	1	String	name of feature\n#fields	quality	1	Float	Quality score of feature\n#fields	filter	1	String	Filter value\n#fields	alleleSeq1	1	String	allele present on first chromosome/haplotype\n#fields	alleleSeq2	1	String	allele present on second chromosome/haplotype\n#fields	sequenced	1	String	sequenced status: v = variant, r = reference (i.e. not this variant), u = unsequenced\n#fields	zyg	1	String	Zygosity status: m = homozygous, t = heterozygous, r = reference, o = other variant, v = variant but genotype unspecified, c = compound (i.e. genotype has this variant and other variant), u = unsequenced\n#fields	phased	1	Integer	Phased status: 0 if not phased, other integer if phased\n#fields	genotypes	H	Integer	Genotypes\n#fields	alleledepth_ref	1	Integer	Allelic depths for the ref allele\n#fields	alleledepth	A	Integer	Allelic depths for the alt alleles in the order listed\n#fields	frequency	A	Float	Allele Frequency\n");
	header = DStringNewFromChar("chromosome\tbegin\tend\tstrand\tqname\tqstart\tqend\tmapquality\tref2\tbegin2\tstrand2\ttlen\tpair\tproperpair\tunmapped\tmateunmapped\tread\tsecondary\tqcfail\tduplicate\tsupplementary\tcigar\tseqlen\tseq\tquality\tother");
	if (argc == 1) {
		hashtable = NULL;
		fpos = 0;
	} else {
		char *field,*type,*descr,*typestring;
		int descrsize,i,size;
		hashtable = hash_init_size(500);
		/* intentionally not using 0, so dstring_hash_get can return NULL if not present */
		fpos=1;
		cur = fields;
		i = 2;
		while (i <= argc) {
			cur = argv[i-1];
			field = cur;
			while (*cur != '\0' && *cur != ':' && *cur != ' ') cur++;
			if (*cur != '\0') cur++;
			type = cur;
			while (*cur != '\0' && *cur != ':' && *cur != ' ') cur++;
			if (*cur != '\0') cur++;
			descr = cur;
			while (*cur != '\0' && *cur != ':') cur++;
			descrsize = cur - descr;
			if (*cur != '\0') cur++;
			dstring_hash_set(hashtable,DStringNewFromCharS(field,type-field-1),(void *)(long int)fpos++);
			if (*type == 'i') {
				typestring = "Integer";
			} else if (*type == 'Z') {
				typestring = "String";
			} else if (*type == 'f') {
				typestring = "Float";
			} else {
				typestring = "String";
			}
			size = (type-field-1);
			fprintf(stdout,"#fields\t%*.*s\t1\t%s\t%*.*s\n",
				size,size,field,
				typestring,
				descrsize,descrsize,descr);
			DStringAppend(header,"\t");
			DStringAppendS(header,field,type-field-1);
			i++;
		}
		fieldsa = DStringArrayNew(fpos+1);

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
	DStringAppend(header,"\n");
	DStringputs(header,stdout);
	while (1) {
		DString *key = DStringNew();
		char *field;
		int i,other;
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
		i = 1;
		while (i <= fpos) {
			DStringSetS(fieldsa->data+i,"",0);
			i++;
		}
		keep = cur;
		other = 0;
		while (*cur != '\0') {
			char *value;
			void *v = NULL;
			field = cur;
			while (*cur != ':' && *cur != '\0') cur++;
			DStringSetS(key,field,cur-field);
			if (hashtable != NULL) {
				v = dstring_hash_get(hashtable,key);
				/* skip type */
				if (*cur != '\0') cur++;
				while (*cur != ':' && *cur != '\0') cur++;
				if (*cur != '\0') cur++;
			}
			/* process value */
			if (v != NULL) {
				i = (int)(long int)v;
				value = cur;
				while (*cur != '\t' && *cur != '\0') cur++;
				DStringSetS(fieldsa->data + i,value,cur-value);
				if (*cur != '\0') cur++;
			} else {
				if (other) fputc(' ',stdout);
				other = 1;
				cur = keep;
				while (*cur != '\t' && *cur != '\0') {
					fputc(*cur,stdout);
					cur++;
				}
				if (*cur != '\0') cur++;
			}
			keep = cur;
		}
		i = 1;
		while (i <= fpos) {
			fputc('\t',stdout);
			DStringputs(fieldsa->data + i,stdout); /* qual */
			i++;
		}
		fputc('\n',stdout);
		if (DStringGetTab(line1,f1,10,result1,0,&numfields)) break;
	}
	FCLOSE(f1);
	if (line1) {DStringDestroy(line1);}
	if (linekeep) {DStringDestroy(linekeep);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
