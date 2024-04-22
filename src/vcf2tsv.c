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
#include <math.h>
#include <ctype.h>
#include "tools.h"
#include "debug.h"
#include "hash.h"

/* 
for later: conversion of genotype list fiels
If A is the allele in REF and B,C,... are the alleles as ordered in ALT, the ordering of genotypes for the likelihoods is given by:
F(j/k) = (k*(k+1)/2)+j. 
In other words, for biallelic sites the ordering is: AA,AB,BB; 
for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc.
*/

typedef struct AltVar {
	DString *type;
	int begin;
	int end;
	char *ref;
	int refsize;
	char *alt;
	int altsize;
	DString *buffer;
	int curallele;
	int refout;
} AltVar;

typedef struct OBufferBucket {
	DString *string;
	int begin;
	int end;
	int typepos;
} OBufferBucket;

typedef struct OBuffer {
	OBufferBucket **buckets;
	int memsize;
	int cursize;
} OBuffer;

/* i know it is ugly in globals, but dont have the time to transfer all, or put in a nice struct at do at the moment */
static OBuffer *obuffer = NULL;
static AltVar *altvars = NULL;
static DString *line, *string, *temp;
static DString *reftype,*snptype, *deltype, *instype, *subtype, *duptype, *invtype, *cnvtype, *transtype, *bndtype, *unknowntype;
static DString *geno, *outinfo, *formatfieldsnumber, *infofieldsnumber;
static DStringArray *header, *format, *info, *samples;
static DStringArray *formatfields , *infofields;
static DString *ref, *alt, *id;
static DString *num;
DString *svlens=NULL,*svtype=NULL;
char *cursvlen=NULL;
static char *typelist;
static int *order = NULL;
static int numalleles,altvarsmax=0;
/* svendpos and svend captures position and value of "END" in info fields, it is also used for gvcf <NON_REF> calls */
static int linenr,svendpos=-1,svlenpos=-1,svtypepos=-1,chr2pos=-1,pos2pos=-1,strandspos=-1;
int svend=-1, svlen=0;
static int ADpos = -1, MIN_DPpos = -1, DPpos = -1, GQpos = -1, GQXpos = -1;

#define a_chrom(array) (DStringArrayGet(array,0))
#define a_pos(array) (DStringArrayGet(array,1))
#define a_id(array) (DStringArrayGet(array,2))
#define a_ref(array) (DStringArrayGet(array,3))
#define a_alt(array) (DStringArrayGet(array,4))
#define a_qual(array) (DStringArrayGet(array,5))
#define a_filter(array) (DStringArrayGet(array,6))
#define a_info(array) (DStringArrayGet(array,7))
#define a_format(array) (DStringArrayGet(array,8))
#define a_geno(array) (DStringArrayGet(array,9))

OBuffer *createbuffer(OBuffer *obuffer, int size) {
	if (obuffer == NULL) {
		obuffer = (OBuffer *)malloc(sizeof(OBuffer));
		obuffer->memsize = 0;
		obuffer->cursize = 0;
		obuffer->buckets = (OBufferBucket **)malloc(sizeof(OBufferBucket *));
	}
	if (size <= obuffer->memsize) {
		return obuffer;
	}
	obuffer->buckets = (OBufferBucket **)realloc(obuffer->buckets,size*sizeof(OBufferBucket *));
	while (obuffer->memsize < size) {
		OBufferBucket *bucket = (OBufferBucket *)malloc(sizeof(OBufferBucket));
		bucket->begin = 0;
		bucket->string = DStringNew();
		/* DStringSetSize(bucket->string,1000); */
		obuffer->buckets[obuffer->memsize] = bucket;
		obuffer->memsize++;
	}
	return obuffer;
}

OBufferBucket *obuffer_getbucket(OBuffer *obuffer) {
	OBufferBucket *bufferbucket;
	if (obuffer->cursize == obuffer->memsize) {
		createbuffer(obuffer,obuffer->cursize+2);
	}
	bufferbucket = obuffer->buckets[obuffer->cursize];
	obuffer->cursize++;
	return bufferbucket;
}

DString *extractID(DString *string,DString *id) {
	char *str, *end;
	str = strstr(string->string,"ID=");
	if (str == NULL) {
		DStringSetS(id,"",0);
		return id;
	}
	str+= 3;
	end = strstr(str,",");
	DStringSetS(id,str,end-str);
	return id;
}

char extractNumber(DString *string) {
	char *str;
	str = strstr(string->string,"Number=");
	if (str == NULL) {
		return '.';
	}
	str+= 7;
	return str[0];
}

void getfield(DString *result,char *list, int pos) {
	char *cur = list;
	int size = 0;
	while (pos > 0) {
		if (*cur == '\0') break;
		if (*cur == ';') break;
		if (*cur == ',') {
			pos--;
		}
		cur++;
	}
	result->string = cur;
	while (*cur != '\0' && *cur != ',' && *cur != ';') {
		cur++;
		size++;
	}
	result->size = size;
}

char numberfromid(DString *id, char num, char *typelist) {
	char *cur, *idc, def;
	int idsize;
	if (num != '.') {
		return num;
	}
	cur = typelist;
	def = *cur;
	/* if typelist is empty, just return num (.) */
	if (*cur++ == '\0') {return num;}
	/* if typelist only contains a default, just return it */
	if (*cur++ == '\0') {return def;}
	idc = id->string;
	idsize = id->size;
	while (*cur != '\0') {
		if (*idc != *cur || idsize == 0) {
			if (idsize == 0 && *cur == ' ') {
				return *(++cur);
			}
			while (*cur != '\0' && *cur != ' ') cur++;
			if (*cur == '\0') break;
			cur++;
			if (*cur == '\0') break;
			cur++;
			if (*cur == '\0') break;
			cur++;
			idc = id->string;
			idsize = id->size;
		} else {
			cur++; idc++; idsize--;
		}
	}
	return def ;
}

void changetoupper(DString *ds) {
	char *cur = ds->string;
	int count = ds->size;
	while (count--) {
		if (*cur > 96) {
			*cur = toupper(*cur);
		}
		cur++;
	}
}

void printallele(DString *bufferstring,int size,char *string,int limit,int refout) {
	if (size == 0) {
		DStringAppendS(bufferstring,"\t",1);
	} else if (refout) {
		DStringAppendS(bufferstring,"\t",1);
		DStringAppendS(bufferstring,string,1);
	} else if (string == NULL || (limit > 0 && size >= limit)) {
		DStringPrintf(bufferstring,"\t%d", size);
	} else {
		DStringAppendS(bufferstring,"\t",1);
		DStringAppendS(bufferstring,string,size);
	}
}

void fprintallele(FILE *fo,char *allele,int size) {
	if (size == 0) {
		return;
	} else if (size > 20 || allele == NULL) {
		fprintf(fo,"%d", size);
	} else if (size == 9 && strncmp(allele,"<NON_REF>",9) == 0) {
		fputc_unlocked('.',fo);
	} else {
		fprintf(fo,"%*.*s", size,size,allele);
	}
}

int output_comparator(const void *p1, const void *p2)
{
	/* 
	<0 The element pointed by p1 goes before the element pointed by p2
	0  The element pointed by p1 is equivalent to the element pointed by p2
	>0 The element pointed by p1 goes after the element pointed by p2
	*/
	const OBufferBucket *bucket1 = *(const OBufferBucket **)p1;
	const OBufferBucket *bucket2 = *(const OBufferBucket **)p2;
/*
fprintf(stderr,"bucket1: begin:%d end:%d typepos: %d \n",bucket1->begin,bucket1->end,bucket1->typepos); fflush(stderr);
fprintf(stderr,"bucket1: string:%*.*s\n",bucket1->string->size,bucket1->string->size,bucket1->string->string); fflush(stderr);
fprintf(stderr,"bucket2: begin:%d end:%d typepos: %d \n",bucket2->begin,bucket2->end,bucket2->typepos); fflush(stderr);
fprintf(stderr,"bucket2: string:%*.*s\n",bucket2->string->size,bucket2->string->size,bucket2->string->string); fflush(stderr);
*/
	if (bucket1->begin < bucket2->begin) {
		return -1;
	} else if (bucket1->begin > bucket2->begin) {
		return 1;
	} else if (bucket1->end < bucket2->end) {
		return -1;
	} else if (bucket1->end > bucket2->end) {
		return 1;
	}
	return naturalcompare(
		bucket1->string->string + bucket1->typepos, bucket2->string->string + bucket2->typepos,
		bucket1->string->size - bucket1->typepos, bucket2->string->size - bucket2->typepos
	);
}

void process_line_parse_info(DStringArray *linea) 	{
	DString *info = a_info(linea);
	char *cur,*curend;
	int i,pos;
	for (i = 0 ; i< infofields->size ; i++) {
		outinfo[i].size = -1;
	}
	cur = info->string;
	curend = cur;
	svlens=(DString *)NULL; svtype=(DString *)NULL; cursvlen=(char *)NULL;
	while (1) {
		while (*curend != '=' && *curend != ';' && *curend != '\0') curend++;
		if (curend == cur) {
			cur = ++curend;
			continue;
		}
		pos = DStringArraySearch(infofields,cur,curend-cur);
		if (*curend == '=') {curend++;}
		cur = curend;
		if (*curend) {while (*curend != ';' && *curend != '\0') curend++;}
		if (pos != -1) {
			outinfo[pos].string = cur;
			outinfo[pos].size = curend-cur;
			if (pos == svtypepos) {
				svtype = outinfo + svtypepos;
			} else if (pos == svendpos) {
				/* vcf 1 based, but for SVs the pos is given before variant and END is defined as POS + length of REF allele - 1 */
				/* so no -1 to get correct endpoint in half open */
				svend=atoi(outinfo[svendpos].string);
			} else if (pos == svlenpos) {
				svlens = outinfo + svlenpos;
			}
		}
		if (*curend == '\0') break;
		cur = ++curend;
	}
	if (svlens != NULL)  {
		cursvlen = svlens->string;
	} else {
		cursvlen = NULL;
	}
}

int process_line_parse_alts(DStringArray *linea,DStringArray *alts,int refout,char locerror,int skiprefindels) {
	int l1,l2,curallele,chrpos,pos,i;
	int snpfound = 0;
	numalleles = alts->size;
	if (numalleles > altvarsmax) {
		altvars = (AltVar *)realloc(altvars,numalleles*sizeof(AltVar));
		altvarsmax = numalleles;
	}
	chrpos = atoi(a_pos(linea)->string);
	chrpos--; /* vcf 1 based */
	for (curallele = 0 ; curallele < numalleles; curallele++) {
		DString *altallele = DStringArrayGet(alts,curallele);
		AltVar *altvar = altvars+curallele;
		char *curref, *curalt;
		int firstremoved = 0;
		/* in vcf genotypes alt alleles number from 1 (0 is ref)*/
		altvar->curallele = curallele + 1;
		/* determine type for this altallele */
		if (cursvlen != NULL)  {
			svlen = atoi(cursvlen); while (*cursvlen != ',' && *cursvlen != '\0') cursvlen++;
		}
		curref = ref->string;
		pos = chrpos;
		l1 = ref->size;
		curalt = altallele->string;
		l2 = altallele->size;
		/* if ((l2 == 0 || *curalt != '<') && (l2 < 2 || (curalt[1] != '<'))) */
		/* remove first base at start if equal (vcf start base for indels) */
		if (l1 > 0 && l2 > 0 && (
			*curref == *curalt
			|| (*curref == 'N' && *curalt >= 'A' && *curalt <= 'Z')
			|| (*curalt == 'N' && *curref >= 'A' && *curref <= 'Z')
		)) {
			pos++;
			curref++; l1--;
			curalt++; l2--;
			firstremoved = 1;
		}
		/* remove extra bases at end (from other overlapping alleles); do first to keep left aligned */
		while (l1 > 0 && l2 > 0 && curref[l1-1] == curalt[l2-1]) {
			l1--; l2--;
		}
		/* remove equal bases at start -> potentially change sub to indel */
		while (*curref == *curalt && l1 > 0 && l2 > 0) {
			/* remove same base before indels */
			pos++;
			curref++; l1--;
			curalt++; l2--;
		}
		altvar->ref = curref; altvar->refsize = l1;
		altvar->alt = curalt; altvar->altsize = l2; altvar->buffer = NULL;
		altvar->refout = 0;
		if (*curalt == '.' || (*curalt == '<' && strcmp(curalt+1,"NON_REF>") == 0)) {
			AltVar *temp = altvars;
			if (numalleles > 1) {
				if (refout) {
					altvar->refout = 1;
					int count = numalleles-1;
					while(count--) {
						if (temp->begin == chrpos && temp->end == chrpos+1) {
							/* don't use if overlapping snp exists */
							altvar->refout = -1; break;
						}
						temp++;
					}
				} else {
					altvar->refout = -1;
				}
				altvar->type = snptype;
				altvar->begin = pos;
				altvar->end = pos+1;
				altvar->altsize = 1;
				altvar->alt[0] = '.';
			} else if (skiprefindels && l1 != 1) {
				/* sam varall sometimes contains ref INDEL (alt =.)
				   they overlap with "snp" refs, causing less good results when used as a varall
				   this option to remove these
				*/
				altvar->refout = -1;
			} else {
				altvar->type = reftype;
				altvar->begin = pos;
				altvar->end = pos+l1;
				altvar->altsize = 1;
				altvar->alt[0] = '.';
				if (svend != -1) {
					if (svend < pos) {
						if (locerror == 'e') {
							fprintf(stderr,"END position %d < begin %d\n",svend,pos);
							exit(1);
						} else if (locerror == 'c') {
							svend = pos;
						}
					}
					altvar->end = svend;
					altvar->refsize = svend - altvar->begin;
					if (altvar->refsize != 1) {altvar->ref = NULL;}
				}
				if (altvar->type == reftype && altvar->end - altvar->begin == 1) {
					altvar->type = snptype;
				}
			}
		} else if ((l2 == 0 || *curalt != '<') && svtype == NULL) {
			if (l1 == 1 && l2 == 1) {
				altvar->type = snptype;
				altvar->begin = pos;
				altvar->end = pos + 1;
				if (pos < snpfound) {snpfound = pos;}
			} else if (l1 == 0 && l2 != 0) {
				altvar->type = instype;
				altvar->begin = pos;
				altvar->end = pos;
			} else if (l1 != 0 && l2 == 0) {
				altvar->type = deltype;
				altvar->begin = pos;
				altvar->end = pos+l1;
			} else {
				altvar->type = subtype;
				altvar->begin = pos;
				altvar->end = pos+l1;
			}
		} else {
			/* sv */
			if (altallele->string[0] == '<') {
				/* vcf uses base before actual variant, but this notation has no matching base in alt to shift it */
				pos++; 
			}
			altvar->begin = pos;
			if (svend != -1) {
				if (svend < pos) {
					if (locerror == 'e') {
						fprintf(stderr,"END position %d < begin %d\n",svend,pos);
						exit(1);
					} else if (locerror == 'c') {
						svend = pos;
					}
				}
				if (svtype != NULL && svtype->size == 3 && strncmp(svtype->string,"INS",3) == 0) {
					/* some output of sniffles has (incorrectly) END != pos for insertions
					 -> fix
					*/
					svend = pos;
				}
				altvar->end = svend;
				altvar->refsize = svend - pos;
			} else {
				altvar->end = pos;
				altvar->refsize = 0;
			}
			altvar->ref = NULL;
			if (l2 > 3 && strncmp(curalt+1,"DEL",3) == 0) {
				altvar->type = deltype;
				altvar->altsize = 0;
			} else if (l2 > 3 && strncmp(curalt+1,"INS",3) == 0) {
				altvar->type = instype;
				if (altvar->end > altvar->begin) {
					altvar->end = altvar->begin;
					altvar->refsize = 0;
				}
				if (svlen > 0) {
					altvar->altsize = svlen;
					altvar->alt = NULL;
				}
			} else if (l2 > 3 && strncmp(curalt+1,"DUP",3) == 0) {
				altvar->type = duptype;
				altvar->altsize = altvar->refsize + svlen;
				altvar->alt = NULL;
			} else if (l2 > 3 && strncmp(curalt+1,"INV",3) == 0) {
				altvar->type = invtype;
				altvar->altsize = 1;
				altvar->alt = "i";
			} else if (l2 > 3 && strncmp(curalt+1,"CNV",3) == 0) {
				altvar->type = cnvtype;
				altvar->altsize = 4;
				altvar->alt = "cnv";
			} else if (l2 > 3 && (strncmp(curalt+1,"TRA",3) == 0 || strncmp(curalt+1,"CTX",3) == 0 || strncmp(curalt+1,"BND",3) == 0)) {
				char *curchr; int cursize;
				altvar->end = altvar->begin;
				altvar->refsize = 0;
				if (curalt[1] == 'B') {
					altvar->type = bndtype;
					if (pos2pos != -1) {
						if (outinfo[pos2pos].size > 0) {
							svend=atoi(outinfo[pos2pos].string);
						}
					}
				} else {
					altvar->type = transtype;
				}
				if (altvar->buffer == NULL) {altvar->buffer = DStringNew();}
				curchr = outinfo[chr2pos].string; cursize = outinfo[chr2pos].size;
				if (cursize > 3 && (curchr[0] == 'c' || curchr[0] == 'C') && curchr[1] == 'h' && curchr[2] == 'r') {
					curchr += 3; cursize -= 3;
				}
				if (strandspos != -1 && outinfo[strandspos].size != -1) {
					/*
					------A------|-----B----- chr1
					------C------|-----D----- chr2
					AD +-	t[p[	piece extending to the right of p is joined after t
					AC ++	t]p]	reverse comp piece extending left of p is joined after t
					CB -+	]p]t	piece extending to the left of p is joined before t
					BD --	[p[t	reverse comp piece extending right of p is joined before t
					*/
					char *s = outinfo[strandspos].string;
					if (s[0] == '+') {
						if (s[1] == '-') {
							DStringPrintf(altvar->buffer,".[%*.*s:%d[",cursize,cursize,curchr,svend-1);
						} else {
							DStringPrintf(altvar->buffer,".]%*.*s:%d]",cursize,cursize,curchr,svend-1);
						}
					} else {
						if (s[1] == '+') {
							DStringPrintf(altvar->buffer,"]%*.*s:%d].",cursize,cursize,curchr,svend-1);
						} else {
							DStringPrintf(altvar->buffer,"[%*.*s:%d[.",cursize,cursize,curchr,svend-1);
						}
					}
				} else {
					DStringPrintf(altvar->buffer,"[%*.*s:%d[",cursize,cursize,curchr,svend-1);
				}
				altvar->altsize = altvar->buffer->size;
				altvar->alt = altvar->buffer->string;
			} else if (svtype != NULL) {
				char *pos,dir = 'r';
				if (svtype->size == 3 && strncmp(svtype->string,"BND",3) == 0) {
					altvar->type = bndtype;
					altvar->refsize = 0;
				} else {
					for (i=0 ; i < svtype->size ; i++) {svtype->string[i] = tolower(svtype->string[i]);}
					altvar->type = svtype;
				}
				pos = strchr(altvar->alt,'[');
				if (pos == NULL) {
					pos = strchr(altvar->alt,']');
					dir = 'l';
				}
				if (pos != NULL) {
					int size, loc, sizeleft;
					pos = strchr(altvar->alt,':');
					size = pos-altvar->alt;
					loc = atoi(pos+1);
					if (dir == 'r') {loc--;}
					if (altvar->buffer == NULL) {altvar->buffer = DStringNew();}
					pos++;
					while(*pos >= 48 && *pos <= 57) pos++;
					sizeleft = altvar->altsize - (pos - altvar->alt);
					if (firstremoved) {
						DStringPrintf(altvar->buffer,".%*.*s:%d%*.*s",size,size,altvar->alt,loc,sizeleft,sizeleft,pos);
					} else {
						DStringPrintf(altvar->buffer,"%*.*s:%d%*.*s.",size,size,altvar->alt,loc,sizeleft,sizeleft,pos);
					}
					altvar->altsize = altvar->buffer->size;
					altvar->alt = altvar->buffer->string;
				} else {
					altvar->altsize = l2;
					altvar->alt = curalt;
				}
			} else {
				/* END without SVTYPE, nor any known variant type in ALT */
				altvar->type = unknowntype;
				altvar->altsize = l2;
				altvar->alt = curalt;
			}
		}
	}
	return(chrpos);
}

void process_line_unsplit(FILE *fo,DStringArray *linea,int excludename,int excludefilter,int refout,char locerror,int skiprefindels) {
	DStringArray *lineformat,*alts;
	DString *ds,*type;
	char zyg, refch;
	int curallele,pos,i;
	int l1,l2,len,igeno,isample,diffchar,diff=0,begin,end;
	svend=-1; svlen=0;
	lineformat = DStringArrayFromChar(a_format(linea)->string,':');
	/* set genos [lrange $line 9 end] */
	ref = a_ref(linea);
	changetoupper(ref);
	changetoupper(a_alt(linea));
	alts = DStringArrayFromChar(a_alt(linea)->string,',');
	/* prepare */
	/* parse info first (need END,SVTYPE and SVLEN for structural variants) */
	/* initialises (global) outinfo, svlens, svtype, cursvlen */
	process_line_parse_info(linea);
	if (svlens != NULL || svtype != NULL ||svend != -1) {
		AltVar *altvar1,*altvar;
		/* structural variants or "END" present, other processing needed*/
		/* determine alt alleles (type, begin, etc ..., store in altvars) */
		/*   altvars: data on alt vars in VarAlt structure, order in vcf list) */
		process_line_parse_alts(linea,alts,refout,locerror,skiprefindels);
		altvar1=altvars+0;
		fprintf(fo,"%s\t%d\t%d\t%*.*s\t", a_chrom(linea)->string, altvar1->begin, altvar1->end, altvar1->type->size, altvar1->type->size, altvar1->type->string);
		if (altvar1->refout) {
			fprintallele(fo,altvar1->ref,altvar1->refsize);
			fputc_unlocked('\t',fo);
			fputc_unlocked('.',fo);
		} else {
			fprintallele(fo,altvar1->ref,altvar1->refsize);
			fputc_unlocked('\t',fo);
			fprintallele(fo,altvar1->alt,altvar1->altsize);
			for (curallele = 1 ; curallele < numalleles; curallele++) {
				altvar = altvars+curallele;
				if (altvar->refsize != altvar1->refsize) {
					fprintf(stderr,"error in alt alleles (different refs on same line not allowed in unsplit): "); DStringPrintTab(stderr,line);	fprintf(stderr,"\n"); exit(1);
				}
				fputc_unlocked(',',fo);
				fprintallele(fo,altvar->alt,altvar->altsize);
			}
		}
	} else {
		/* no structural variants present, needs other (simpler) processing */
		/* will use l1 and l2 to determine type, make sub if one of the alleles differs */
		numalleles = alts->size;
		/* check first base */
		type = NULL;
		pos = atoi(a_pos(linea)->string);
		refch = ref->string[0];
		l1 = ref->size;
		l2 = 0;
		diffchar = 0;
		for (curallele = 0 ; curallele < numalleles; curallele++) {
			DString *temp = DStringArrayGet(alts,curallele);
			if (temp->size > l2) {l2 = temp->size;}
			if (temp->string[0] != refch) {diffchar = 1;}
		}
		if (l1 == 1 && l2 == 1) {
			diff = 0;
			begin = pos - 1;
			end = pos;
			type = snptype;
		} else if (diffchar) {
			diff = 0;
			begin = pos - 1;
			end = pos + l1 - 1;
			type = subtype;
		} else {
			diff = 1;
			begin = pos;
			end = pos + l1 - 1;
			if (l1 == 1) {
				type = instype;
			} else if (l2 == 1) {
				type = deltype;
			} else {
				type = subtype;
			}
		}
		fprintf(fo,"%s\t%d\t%d\t%s", a_chrom(linea)->string, begin, end, type->string);
		fputc_unlocked('\t',fo);
		fprintallele(fo,ref->string+diff,ref->size-diff);
		if (DStringArrayGet(alts,0)->size == 0 && diff > 0) {
			fprintf(stderr,"error in alt alleles: "); DStringPrintTab(stderr,line);	fprintf(stderr,"\n"); exit(1);
		}
		ds = DStringArrayGet(alts,0);
		fputc_unlocked('\t',fo);
		fprintallele(fo,ds->string+diff,ds->size-diff);
		for (curallele = 1 ; curallele < numalleles; curallele++) {
			DString *temp = DStringArrayGet(alts,curallele);
			if (temp->size == 0 && diff > 0) {
				fprintf(stderr,"error in alt alleles: "); DStringPrintTab(stderr,line);	fprintf(stderr,"\n"); exit(1);
			}
			fputc_unlocked(',',fo);
			fprintallele(fo,temp->string+diff,temp->size-diff);
		}
	}
	if (!excludename) fprintf(fo,"\t%s", a_id(linea)->string);
	fprintf(fo,"\t%s", a_qual(linea)->string);
	if (!excludefilter) fprintf(fo,"\t%s", a_filter(linea)->string);
	if (header->size >= 9) {
		NODPRINT("==== Determine geno fields order ====")
		order = realloc(order,formatfields->size*sizeof(int));
		for (i = 0 ; i < formatfields->size ; i++) {
			order[i] = DStringArraySearch(lineformat,DStringArrayGet(formatfields,i)->string,DStringArrayGet(formatfields,i)->size);
		}
		NODPRINT("==== Process genos ====")
		igeno = 9; /* genos start at col 9 */
		for (isample = 0 ; isample < samples->size ; isample++) {
			DStringArray *genoa;
			DString *genotype;
			char *genotypestring,*genotypecur;
			int phased,i,a1,a2;
			geno = DStringArrayGet(linea,igeno++);
			if (geno->size == 1 && geno->string[0] == '.') {
				fprintf(fo,"\t?\t?\t?\t?\t?");
				for (i = 1 ; i < formatfields->size; i++) {
					fputc_unlocked('\t',fo);
					if (formatfieldsnumber->string[i] == 'R') {
						fputc_unlocked('\t',fo);
					}
				}
			} else {
				genoa = DStringArrayFromChar(geno->string,':');
				if (order[0] == -1) {
					genotypestring = "0/0";
				} else {
					genotype = DStringArrayGet(genoa,order[0]);
					genotypestring = genotype->string;
				}
				genotypecur = genotypestring;
				while (*genotypecur != '|' && *genotypecur != '/' && *genotypecur != '\0') {
					genotypecur++;
				}
				if (*genotypecur == '|') {phased = 1;} else {phased = 0;}
				/* print out alleleSeq1 and alleleSeq2 */
				/* ----------------------------------- */
				fputc_unlocked('\t',fo);
				if (genotypestring[0] == '.') {
					fputc_unlocked('?',fo);
					a1 = -2;
				} else if (genotypestring == genotypecur) {
					/* empty genotype */
					fputc_unlocked('-',fo);
					a1 = -1;
				} else {
					a1 = atol(genotypestring);
					if (svlens != NULL || svtype != NULL || svend != -1) {
						if (a1 == 0) {
							fprintallele(fo,altvars[0].ref,altvars[0].refsize);
						} else {
							AltVar *temp = altvars+a1-1;
							if (a1 > numalleles) {
								fprintf(stderr,"allele %d does not exist in line %s\n",a1,a_format(linea)->string);
								exit(1);
							}
							fprintallele(fo,temp->alt,temp->altsize);
						}
					} else if (a1 == 0) {
						fprintallele(fo,ref->string+diff,ref->size-diff);
					} else {
						ds = DStringArrayGet(alts,a1-1);
						if (ds->size == 0 && diff > 0) {
							fprintf(stderr,"error in alt alleles: "); DStringPrintTab(stderr,line);	fprintf(stderr,"\n"); exit(1);
						}
						fprintallele(fo,ds->string+diff,ds->size-diff);
					}
				}
				if (*genotypecur != '\0') {genotypecur++;}
				fputc_unlocked('\t',fo);
				if (*genotypecur == '.') {
					fputc_unlocked('?',fo);
					a2 = -2;
				} else if (*genotypecur == '\0') {
					fputc_unlocked('-',fo);
					a2 = -1;
				} else {
					a2 = atol(genotypecur);
					if (svlens != NULL || svtype != NULL || svend != -1) {
						if (a2 == 0) {
							fprintallele(fo,altvars[0].ref,altvars[0].refsize);
						} else {
							AltVar *temp = altvars+a2-1;
							if (a1 > numalleles) {
								fprintf(stderr,"allele %d does not exist in line %s\n",a1,a_format(linea)->string);
								exit(1);
							}
							fprintallele(fo,temp->alt,temp->altsize);
						}
					} else if (a2 == 0) {
						fprintallele(fo,ref->string+diff,ref->size-diff);
					} else {
						ds = DStringArrayGet(alts,a2-1);
						if (ds->size == 0 && diff > 0) {
							fprintf(stderr,"error in alt alleles: "); DStringPrintTab(stderr,line);	fprintf(stderr,"\n"); exit(1);
						}
						fprintallele(fo,ds->string+diff,ds->size-diff);
					}
				}
				/* print out zyg */
				/* ------------- */
				if ((a1 < 0 && a2 < 0) || (a1 < 0 && a2 == 0) || (a1 == 0 && a2 < 0)) {
					zyg = 'u';
				} else if (a1 < 0 || a2 < 0) {
					zyg = 'v';
				} else if (a1 > 0) {
					if (a2 == a1) {
						zyg = 'm';
					} else if (a2 > 0) {
						zyg = 'c';
					} else {
						zyg = 't';
					}
				} else if (a2 > 0) {
					zyg = 't';
				} else {
					zyg = 'r';
				}
				fprintf(fo,"\t%c",zyg);
				fprintf(fo,"\t%d",phased);
				/* print out genotypes */
				/* ------------------- */
				fputc_unlocked('\t',fo);
				genotypecur = genotypestring;
				while (*genotypecur != '\0') {
					if (*genotypecur == '|') {
						fputc_unlocked(',',fo);
					} else if (*genotypecur == '/') {
						fputc_unlocked(';',fo);
					} else if (*genotypecur == '.') {
						fputc_unlocked('?',fo);
					} else {
						fputc_unlocked(*genotypecur,fo);
					}
					genotypecur++;
				}
				/* print out rest of the fields */
				/* ---------------------------- */
				for (i = 1 ; i < formatfields->size; i++) {
					if (order[i] <= 0 || order[i] >= genoa->size) {
						fprintf(fo,"\t");
						if (formatfieldsnumber->string[i] == 'R') {
							fprintf(fo,"\t");
						}
					} else if (formatfieldsnumber->string[i] == 'R') {
						DString *value;
						char *cur,test=',';
						int count;
						value = DStringArrayGet(genoa,order[i]);
						cur = value->string; count = value->size;
						fputc_unlocked('\t',fo);
						while(count--) {
							if (*cur == test) {
								fputc_unlocked('\t',fo); cur++; test = '\0';
							} else {
								fputc_unlocked(*cur++,fo);
							}
						}
					} else if (formatfieldsnumber->string[i] == 'A') {
						DString *value;
						value = DStringArrayGet(genoa,order[i]);
						fprintf(fo,"\t%*.*s",value->size,value->size,value->string);
					} else {
						fprintf(fo,"\t%s",DStringArrayGet(genoa,order[i])->string);
					}
				}
				if (genoa != NULL) DStringArrayDestroy(genoa);
			}
		}
	}
	/* info output */
	{
		char *cur,*curend;
		DString *info = a_info(linea);
		for (i = 0 ; i< infofields->size ; i++) {
			outinfo[i].size = -1;
		}
		cur = info->string;
		curend = cur;
		while (1) {
			while (*curend != '=' && *curend != ';' && *curend != '\0') curend++;
			if (curend == cur) {
				cur = ++curend;
				continue;
			}
			pos = DStringArraySearch(infofields,cur,curend-cur);
			if (*curend == '=') {curend++;}
			cur = curend;
			if (*curend) {while (*curend != ';' && *curend != '\0') curend++;}
			if (pos != -1) {
				outinfo[pos].string = cur;
				outinfo[pos].size = curend-cur;
			}
			if (*curend == '\0') break;
			cur = ++curend;
		}
		for (i = 0 ; i< infofields->size ; i++) {
			if (i == svendpos || i == svlenpos || i == svtypepos) continue;
			len  = outinfo[i].size;
			if (len == -1) {
				fprintf(fo,"\t");
				if (infofieldsnumber->string[i] == 'R') {
					fprintf(fo,"\t");
				}
			} else if (len == 0) {
				if (infofieldsnumber->string[i] == 'R') {
					fprintf(fo,"\t");
				}
				fprintf(fo,"\t1");
			} else if (infofieldsnumber->string[i] == 'R') {
				char *cur,test=',';
				int count;
				cur = outinfo[i].string; count = outinfo[i].size;
				fputc_unlocked('\t',fo);
				while(count--) {
					if (*cur == test) {
						fputc_unlocked('\t',fo); cur++; test = '\0';
					} else {
						fputc_unlocked(*cur++,fo);
					}
				}
			} else if (infofieldsnumber->string[i] == 'A') {
				fprintf(fo,"\t%*.*s",outinfo[i].size,outinfo[i].size,outinfo[i].string);
			} else {
				fprintf(fo,"\t%*.*s",len,len,outinfo[i].string);
			}
		}
		fprintf(fo,"\n");
	}
	DStringArrayDestroy(lineformat);
}

int process_line_split(OBuffer *obuffer,DStringArray *linea,int excludename,int excludefilter,DString *genotypelist,int refout,char locerror,int skiprefindels) {
	DStringArray *alts = NULL;
	DStringArray *lineformat = NULL;
	char *zyg;
	int len,curallele,chrpos,curpos,i,igeno,isample;
	svend=-1; svlen=0;
	/* determine type, ref, alt, ... for different alleles */
	lineformat = DStringArrayFromChar(a_format(linea)->string,':');
	/* set genos [lrange $line 9 end] */
	ref = a_ref(linea);
	changetoupper(ref);
	changetoupper(a_alt(linea));
	alts = DStringArrayFromChar(a_alt(linea)->string,',');
	/* prepare */
	/* parse info first (need END,SVTYPE and SVLEN for structural variants) */
	/* initialises (global) outinfo, svlens, svtype, cursvlen */
	process_line_parse_info(linea);
	/* determine alt alleles (type, begin, etc ..., store in altvars) */
	/*   altvars: data on alt vars in VarAlt structure, order in vcf list) */
	/*   numalleles: number of alternative alleles */
	chrpos = process_line_parse_alts(linea,alts,refout,locerror,skiprefindels);
	/* go over the different alleles, and print out */
	for (curpos = 0 ; curpos < numalleles; curpos++) {
		OBufferBucket *bufferbucket;
		DString *bufferstring;
		AltVar *altvar;
		int ADpos = -1, MIN_DPpos = -1, GQXpos = -1;
		altvar = altvars+curpos;
		curallele = altvar->curallele;
		/* "ref" allele is not printed out	if refout == -1 */
		if (altvar->refout == -1) continue;
		NODPRINT("==== Print variant info to buffer ====")
		bufferbucket = obuffer_getbucket(obuffer);
		bufferbucket->begin = altvar->begin;
		bufferbucket->end = altvar->end;
		bufferbucket->typepos = 0;
		bufferstring = bufferbucket->string;
		DStringSetS(bufferstring,"",0);
		DStringPrintf(bufferstring,"%s\t%d\t%d\t", a_chrom(linea)->string, altvar->begin, altvar->end);
		bufferbucket->typepos = bufferstring->size;
		DStringAppendS(bufferstring,altvar->type->string,altvar->type->size);
		if (altvar->refout) {
			DStringAppendS(bufferstring,"\t",1);
			DStringAppendS(bufferstring,altvar->ref,1);
			DStringAppendS(bufferstring,"\t.",2);
		} else {
			printallele(bufferstring,altvar->refsize,altvar->ref,20,0);
			printallele(bufferstring,altvar->altsize,altvar->alt,0,0);
		}
		if (!excludename) {
			DStringAppendS(bufferstring,"\t",1);
			DStringAppendS(bufferstring,a_id(linea)->string,a_id(linea)->size);
		}
		DStringAppendS(bufferstring,"\t",1);
		if (altvar->refout) {
			DStringAppendS(bufferstring,".",1);
		} else {
			DStringAppendS(bufferstring,a_qual(linea)->string,a_qual(linea)->size);
		}
		if (!excludefilter) {
			DStringAppendS(bufferstring,"\t",1);
			DStringAppendS(bufferstring,a_filter(linea)->string,a_filter(linea)->size);
		}
		if (header->size >= 9) {
			NODPRINT("==== Determine geno fields order ====")
			order = realloc(order,formatfields->size*sizeof(int));
			for (i = 0 ; i < formatfields->size ; i++) {
				DString *string = DStringArrayGet(formatfields,i);
				order[i] = DStringArraySearch(lineformat,string->string,string->size);
				if (string->size == 2 && string->string[0] == 'A' && string->string[1] == 'D') {
					ADpos = i;
				} else if (string->size == 6 && strncmp(string->string,"MIN_DP",6) == 0) {
					MIN_DPpos = i;
				} else if (string->size == 3 && string->string[0] == 'G' && string->string[1] == 'Q' && string->string[2] == 'X') {
					GQXpos = i;
				}
			}
			NODPRINT("==== Process genos ====")
			igeno = 9; /* genos start at col 9 */
			for (isample = 0 ; isample < samples->size ; isample++) {
				DStringArray *genoa;
				DString *genotype;
				char *genotypestring,*genotypecur;
				int phased,i,a1,a2;
				DStringSetS(genotypelist,"",0);
				geno = DStringArrayGet(linea,igeno++);
				if (geno->size == 1 && geno->string[0] == '.') {
					DStringAppendS(bufferstring,"\t?\t?\t?\t?\t?",10);
					for (i = 1 ; i < formatfields->size; i++) {
						DStringAppendS(bufferstring,"\t",1);
						if (formatfieldsnumber->string[i] == 'R') {
							DStringAppendS(bufferstring,"\t",1);
						}
					}
				} else {
					genoa = DStringArrayFromChar(geno->string,':');
					if (order[0] == -1) {
						genotypestring = "0/0";
					} else {
						genotype = DStringArrayGet(genoa,order[0]);
						genotypestring = genotype->string;
					}
					genotypecur = genotypestring;
					while (*genotypecur != '|' && *genotypecur != '/' && *genotypecur != '\0') {
						genotypecur++;
					}
					if (*genotypecur == '|') {phased = 1;} else {phased = 0;}
					/* print out alleleSeq1 and alleleSeq2 */
					/* ----------------------------------- */
					if (genotypestring[0] == '.') {
						DStringAppendS(bufferstring,"\t?",2);
						a1 = -2;
						DStringAppendS(genotypelist,"?",1);
					} else if (genotypestring == genotypecur) {
						/* empty genotype */
						DStringAppendS(bufferstring,"\t-",2);
						a1 = -1;
						DStringAppendS(genotypelist,"?",1);
					} else {
						a1 = atol(genotypestring);
						if (a1 == 0) {
							printallele(bufferstring,altvar->refsize,altvar->ref,20,altvar->refout);
							DStringAppendS(genotypelist,"0",1);
						} else if (a1 == curallele) {
							printallele(bufferstring,altvar->altsize,altvar->alt,0,altvar->refout);
							DStringAppendS(genotypelist,"1",1);
						} else {
							AltVar *temp = altvars+a1-1;
							if ((temp->begin >= altvar->end && temp->end != altvar->begin) || (temp->end <= altvar->begin && temp->begin != altvar->end)) {
								/* if this genotype is not overlapping with current alt (e.g. an insertion is located after a snp), this position is reference */
								a1 = 0;
								printallele(bufferstring,altvar->refsize,altvar->ref,20,altvar->refout);
								DStringAppendS(genotypelist,"0",1);
							} else if (temp->type == altvar->type && temp->begin == altvar->begin && temp->end == altvar->end) {
								printallele(bufferstring,temp->altsize,temp->alt,0,altvar->refout);
								DStringAppendS(genotypelist,"2",1);
							} else {
								DStringAppendS(bufferstring,"\t@",2);
								DStringAppendS(genotypelist,"2",1);
							}
						}
					}
					if (*genotypecur == '|') {
						DStringAppendS(genotypelist,",",1);
					} else if (*genotypecur == '/') {
						DStringAppendS(genotypelist,";",1);
					} else {
						DStringAppendS(genotypelist,"?",1);
					}
					if (*genotypecur != '\0') {genotypecur++;}
					if (*genotypecur == '.') {
						DStringAppendS(bufferstring,"\t?",2);
						a2 = -2;
						DStringAppendS(genotypelist,"?",1);
					} else if (*genotypecur == '\0') {
						DStringAppendS(bufferstring,"\t-",2);
						a2 = -1;
						DStringAppendS(genotypelist,"?",1);
					} else {
						a2 = atol(genotypecur);
						if (a2 == 0) {
							printallele(bufferstring,altvar->refsize,altvar->ref,20,altvar->refout);
							DStringAppendS(genotypelist,"0",1);
						} else if (a2 == curallele) {
							printallele(bufferstring,altvar->altsize,altvar->alt,0,altvar->refout);
							DStringAppendS(genotypelist,"1",1);
						} else {
							AltVar *temp = altvars+a2-1;
							if ((temp->begin >= altvar->end && temp->end != altvar->begin) || (temp->end <= altvar->begin && temp->begin != altvar->end)) {
								/* if this genotype is not overlapping with current alt (e.g. an insertion is located after a snp), this position is reference */
								a2 = 0;
								printallele(bufferstring,altvar->refsize,altvar->ref,20,altvar->refout);
								DStringAppendS(genotypelist,"0",1);
							} else if (temp->type == altvar->type && temp->begin == altvar->begin && temp->end == altvar->end) {
								printallele(bufferstring,temp->altsize,temp->alt,0,altvar->refout);
								DStringAppendS(genotypelist,"2",1);
							} else {
								DStringAppendS(bufferstring,"\t@",2);
								DStringAppendS(genotypelist,"2",1);
							}
						}
					}
					if ((a1 < 0 && a2 < 0) || (a1 < 0 && a2 == 0) || (a1 == 0 && a2 < 0)) {
						zyg = "u";
					} else if (a1 < 0 || a2 < 0) {
						zyg = "v";
					} else if (a1 == curallele) {
						if (a2 == a1) {
							zyg = "m";
						} else if (a2 > 0) {
							zyg = "c";
						} else {
							zyg = "t";
						}
					} else if (a2 == curallele) {
						if (a1 > 0) {
							zyg = "c";
						} else {
							zyg = "t";
						}
					} else if (a1 > 0 || a2 > 0) {
						zyg = "o";
					} else {
						zyg = "r";
					}
					DStringAppendS(bufferstring,"\t",1);
					DStringAppendS(bufferstring,zyg,1);
					DStringAppendS(bufferstring,"\t",1);
					if (phased) {
						DStringAppendS(bufferstring,"1",1);
					} else {
						DStringAppendS(bufferstring,"0",1);
					}
					DStringAppendS(bufferstring,"\t",1);
					DStringAppendS(bufferstring,genotypelist->string,genotypelist->size);
					while (*genotypecur >= 48 && *genotypecur <= 57) {
						genotypecur++;
					}
					while (*genotypecur != '\0') {
						if (*genotypecur == '|') {
							DStringAppendS(bufferstring,",",1);
							genotypecur++;
						} else if (*genotypecur == '/') {
							DStringAppendS(bufferstring,";",1);
							genotypecur++;
						} else if (*genotypecur >= 48 && *genotypecur <= 57) {
							a1 = atol(genotypecur);
							if (a1 == 0) {
								DStringAppendS(bufferstring,"0",1);
							} else if (a1 == curallele) {
								DStringAppendS(bufferstring,"1",1);
							} else {
								DStringAppendS(bufferstring,"2",1);
							}
							while (*genotypecur >= 48 && *genotypecur <= 57) {
								genotypecur++;
							}
						} else {
							DStringAppendS(bufferstring,genotypecur,1);
							genotypecur++;
						}
					}
					for (i = 1 ; i < formatfields->size; i++) {
						if (order[i] <= 0 || order[i] >= genoa->size) {
							/* field is not given */
							if (i == DPpos && ADpos != -1 && order[ADpos] != -1) {
								/* if depth is not given, get depth from AD field: sum af all alleles counts */
								DString *value;
								char *cur, *end;
								int dp = 0;
								value = DStringArrayGet(genoa,order[ADpos]);
								cur = value->string; end = value->string + value->size;
								while (cur < end) {
									dp += atoi(cur);
									while (cur++) {
										if (*cur == ',') {
											cur++;
											break;
										} else if (cur >= end) {
											break;
										}
									}
								}
								DStringPrintf(bufferstring,"\t%d",dp);
							} else if (i == DPpos && MIN_DPpos != -1 && order[MIN_DPpos] != -1) {
								/* if depth is not given, get depth from MIN_DP field: sum af all alleles counts */
								DString *value;
								char *cur, *end;
								int dp = 0;
								value = DStringArrayGet(genoa,order[MIN_DPpos]);
								cur = value->string; end = value->string + value->size;
								while (cur < end) {
									dp += atoi(cur);
									while (cur++) {
										if (*cur == ',') {
											cur++;
											break;
										} else if (cur >= end) {
											break;
										}
									}
								}
								DStringPrintf(bufferstring,"\t%d",dp);
							} else if (i == GQpos && GQXpos != -1 && order[GQXpos] != -1) {
								/* if GQ (genoqual) is not given, get it from GQX */
								DString *value;
								value = DStringArrayGet(genoa,order[GQXpos]);
								DStringAppendS(bufferstring,"\t",1);
								DStringAppendS(bufferstring,value->string,value->size);
							} else {
								DStringAppendS(bufferstring,"\t",1);
								if (formatfieldsnumber->string[i] == 'R') {
									DStringAppendS(bufferstring,"\t",1);
								}
							}
						} else if (formatfieldsnumber->string[i] == 'R') {
							DString result, *value;
							value = DStringArrayGet(genoa,order[i]);
							getfield(&result,value->string,0);
							if (i == ADpos && result.size == 1 && result.string[0] == '0' && altvar->alt[0] == '.') {
								/* fix for some versions of gatk, that incorrectly set AD to 0 for reference call */
								DStringAppendS(bufferstring,"\t",1);
							} else {
								DStringAppendS(bufferstring,"\t",1);
								DStringAppendS(bufferstring,result.string,result.size);
							}
							getfield(&result,value->string,curallele);
							DStringAppendS(bufferstring,"\t",1);
							DStringAppendS(bufferstring,result.string,result.size);
						} else if (formatfieldsnumber->string[i] == 'A') {
							DString result, *value;
							value = DStringArrayGet(genoa,order[i]);
							getfield(&result,value->string,curallele-1);
							DStringAppendS(bufferstring,"\t",1);
							DStringAppendS(bufferstring,result.string,result.size);
						} else {
							DStringAppendS(bufferstring,"\t",1);
							DStringAppend(bufferstring,DStringArrayGet(genoa,order[i])->string);
						}
					}
					if (genoa != NULL) DStringArrayDestroy(genoa);
				}
			}
		}
		/* info output */
		for (i = 0 ; i< infofields->size ; i++) {
			if (i == svendpos || i == svlenpos || i == svtypepos) continue;
			len  = outinfo[i].size;
			if (len == -1) {
				DStringAppendS(bufferstring,"\t",1);
				if (infofieldsnumber->string[i] == 'R') {
					DStringAppendS(bufferstring,"\t",1);
				}
			} else if (len == 0) {
				if (infofieldsnumber->string[i] == 'R') {
					DStringAppendS(bufferstring,"\t",1);
				}
				DStringAppendS(bufferstring,"\t1",2);
			} else if (infofieldsnumber->string[i] == 'R') {
				DString result;
				getfield(&result,outinfo[i].string,0);
				DStringAppendS(bufferstring,"\t",1);
				DStringAppendS(bufferstring,result.string,result.size);
				getfield(&result,outinfo[i].string,curallele);
				DStringAppendS(bufferstring,"\t",1);
				DStringAppendS(bufferstring,result.string,result.size);
			} else if (infofieldsnumber->string[i] == 'A') {
				DString result;
				getfield(&result,outinfo[i].string,curallele-1);
				DStringAppendS(bufferstring,"\t",1);
				DStringAppendS(bufferstring,result.string,result.size);
			} else {
				DStringAppendS(bufferstring,"\t",1);
				DStringAppendS(bufferstring,outinfo[i].string,len);
			}
		}
		DStringAppendS(bufferstring,"\n",1);
	}
	DStringArrayDestroy(lineformat);
	DStringArrayDestroy(alts);
	return(chrpos);
}

void process_print_buffer(FILE *fo,OBuffer *obuffer,int limit) {
	OBufferBucket **curbucket = obuffer->buckets, **targetbucket;
	OBufferBucket *tempbucket;
	int count = obuffer->cursize;
	/* sort output */
	if (count == 0) return;
	if (count > 1) {
		qsort((void *)obuffer->buckets,obuffer->cursize,sizeof(OBufferBucket *),output_comparator);
	}
	while (count) {
		if ((*curbucket)->begin > limit) break;
		DStringputs((*curbucket)->string,fo);
		curbucket++; count--;
	}
	/* move remaining lines in buffer from the end to the start */
	/* moving not needed if nothing was printed */
	if (count == obuffer->cursize) return;
	obuffer->cursize = count;
	/* moving not needed if buffer is empty */
	if (count == 0) return;
	targetbucket = obuffer->buckets;
	while (count--) {
		tempbucket = *targetbucket;
		*targetbucket = *curbucket;
		*curbucket = tempbucket;
		targetbucket++; curbucket++;
	}
}

int main(int argc, char *argv[]) {
	Hash_table *conv_formata, *donefields, *keepfields = NULL, *doneinfo = NULL;
	FILE *fd = NULL, *fo = NULL, *fh = NULL;
	DStringArray *headerfields, *linea;
	DString *genotypelist=DStringNew(), *prevchr = DStringNew(),*dsbuffer = DStringNew();
	char *tempfile = NULL, *meta = NULL, *splitstring;
	char locerror = 'e';
	int *order = NULL;
	int split,refout,read,i,j,maxtab, min, excludefilter = 0, excludename = 0, infopos = 0, skiprefindels = 0;
	line=NULL; string=NULL; temp=NULL;
	reftype=DStringNewFromChar("ref");
	snptype=DStringNewFromChar("snp"); deltype=DStringNewFromChar("del"); instype=DStringNewFromChar("ins"); subtype=DStringNewFromChar("sub");
	duptype=DStringNewFromChar("dup"); invtype=DStringNewFromChar("inv"); cnvtype=DStringNewFromChar("cnv"); transtype=DStringNewFromChar("trans");
	bndtype=DStringNewFromChar("bnd"); unknowntype=DStringNewFromChar("unk");
	geno = NULL; outinfo = NULL; formatfieldsnumber=NULL; infofieldsnumber=NULL;
	header=NULL; format=NULL; info=NULL; samples=NULL;
	formatfields=NULL ; headerfields=NULL; infofields=NULL; linea=NULL;
	ref=NULL; alt=DStringNew(); id=DStringNew();
	num=DStringNew();
	linenr=0; numalleles=1;
	svtypepos = -1; svendpos = -1; svlenpos = -1;
	DStringSetS(num,".",1);
	altvars = (AltVar *)malloc(5*sizeof(AltVar)); altvarsmax = 5;
	/* actually initialise the buffer here */
	obuffer = createbuffer(obuffer,1);
	line = DStringNew();
	if (argc < 2) {
		fprintf(stderr,"Format is: vcf2tsv split typelist ?infile? ?outfile? ?removefields? ?refout? ?keepfields? ?locerror? ?skiprefindels?\n");
		exit(EXIT_FAILURE);
	}
	splitstring = argv[1];
	if (splitstring[0] == 'o') {
		split = 0;
	} else {
		split = atoi(splitstring);
	}
	if (argc >= 3) {
		typelist = argv[2];
	} else {
		typelist = ". AD R RPA R AC A AF A";
	}
	if (argc < 4 || (argv[3][0] == '-' && argv[3][1] == '\0')) {
		fd = stdin;
	} else {
		fd = fopen64_or_die(argv[3],"r");
		if (fd == NULL) {fprintf(stderr,"file %s does not exists\n",argv[2]);exit(1);}
	}
	if (argc < 5 || (argv[4][0] == '-' && argv[4][1] == '\0')) {
		fo = stdout;
	} else {
		fo = fopen64_or_die(argv[4],"w");
	}
	NODPRINT("==== Reading header ====")
	DStringGetLine(line, fd);
	if (line->size < 16 || strncmp(line->string,"##fileformat=VCF",16) != 0) {
		fprintf(stderr,"error: input is not a vcf file\n");
		exit(1);
	}
	temp = DStringNew();
	/* create hash for fields already used */
	donefields = hash_init();
	/* create hash for conversion of field names */
	conv_formata = hash_init();
	dstring_hash_set(conv_formata,DStringNewFromChar("AD"),(void *)DStringNewFromChar("alleledepth"));
	dstring_hash_set(conv_formata,DStringNewFromChar("GT"),(void *)DStringNewFromChar("genotype"));
	dstring_hash_set(conv_formata,DStringNewFromChar("DP"),(void *)DStringNewFromChar("coverage"));
	dstring_hash_set(conv_formata,DStringNewFromChar("FT"),(void *)DStringNewFromChar("gfilter"));
	dstring_hash_set(conv_formata,DStringNewFromChar("GL"),(void *)DStringNewFromChar("loglikelihood"));
	dstring_hash_set(conv_formata,DStringNewFromChar("GQ"),(void *)DStringNewFromChar("genoqual"));
	dstring_hash_set(conv_formata,DStringNewFromChar("PS"),(void *)DStringNewFromChar("phaseset"));
	dstring_hash_set(conv_formata,DStringNewFromChar("HQ"),(void *)DStringNewFromChar("haploqual"));
	dstring_hash_set(conv_formata,DStringNewFromChar("AN"),(void *)DStringNewFromChar("totalallelecount"));
	dstring_hash_set(conv_formata,DStringNewFromChar("AC"),(void *)DStringNewFromChar("allelecount"));
	dstring_hash_set(conv_formata,DStringNewFromChar("AF"),(void *)DStringNewFromChar("frequency"));
	dstring_hash_set(conv_formata,DStringNewFromChar("AA"),(void *)DStringNewFromChar("Ancestralallele"));
	dstring_hash_set(conv_formata,DStringNewFromChar("DB"),(void *)DStringNewFromChar("dbsnp"));
	dstring_hash_set(conv_formata,DStringNewFromChar("H2"),(void *)DStringNewFromChar("Hapmap2"));
	/* create hash for fields to remove */
	if (argc >= 6) {
		char *start = argv[5];
		char *cur = argv[5];
		while (1) {
			int len;
			while (*cur != ' ' && *cur != '\0') {
				cur++;
			}
			len = cur-start;
			if (len == 6 && strncmp(start,"filter",6) == 0) {
				excludefilter = 1;
			} else if (len == 4 && strncmp(start,"name",4) == 0) {
				excludename = 1;
			} else if (len == 2 && strncmp(start,"DP",2) == 0) {
				fprintf(stderr,"Cannot remove field DP\n");
				exit(1);
			} else {
				dstring_hash_set(conv_formata,DStringNewFromCharS(start,len),(void *)DStringNew());
			}
			if (*cur == '\0') break;
			start = ++cur;
		}
	}
	if (argc >= 7) {
		refout = atoi(argv[6]);
	} else {
		refout = 0;
	}
	if (argc >= 8) {
		char *start = argv[7];
		char *cur = argv[7];
		if (*start == '*') {
			keepfields = NULL;
		} else {
			excludefilter = 1; excludename = 1;
			keepfields = hash_init();
			while (1) {
				int len;
				while (*cur != ' ' && *cur != '\0') {
					cur++;
				}
				len = cur-start;
				if (len == 6 && strncmp(start,"filter",6) == 0) {
					excludefilter = 0;
				} else if (len == 4 && strncmp(start,"name",4) == 0) {
					excludename = 0;
				} else {
					dstring_hash_set(keepfields,DStringNewFromCharS(start,len),(void *)DStringNew());
				}
				if (*cur == '\0') break;
				start = ++cur;
			}
			dstring_hash_set(keepfields,DStringNewFromCharS("END",3),(void *)DStringNew());
		}
	}
	if (argc >= 9) {
		locerror = argv[8][0];
	} else {
		locerror = 'e';
	}
	if (argc >= 10) {
		skiprefindels = atoi(argv[9]);
	} else {
		skiprefindels = 0;
	}
	if (argc >= 11) {
		tempfile = argv[10];
	}
	if (tempfile == NULL) {
		tempfile = tempfilename();
	}
	if (argc >= 12) {
		meta = argv[11];
	} else {
		meta="";
	}
	/* prepare file for making header */
	fh = fopen64_or_die(tempfile,"w");
	fprintf(fh,"%s\n",line->string);
	info = DStringArrayNew(10);
	format = DStringArrayNew(10);
	while ((read = DStringGetLine(line, fd)) != -1) {
		linenr++;
		/* fprintf(fo,"%s\n",line->string); */
		fprintf(fh,"%s\n",line->string);
		if (line->string[0] == '#') {
			if (line->string[1] == '#') {
				if (line->size > 9 && strncmp("FORMAT=",line->string+2,7) == 0) {
					if (DStringArraySearch(format,line->string+9,strlen(line->string+9)) == -1) {
						DStringArrayAppend(format,line->string+9,strlen(line->string+9));
					}
				} else if (line->size > 7 && strncmp("INFO=",line->string+2,5) == 0) {
					if (DStringArraySearch(info,line->string+7,strlen(line->string+7)) == -1) {
						DStringArrayAppend(info,line->string+7,strlen(line->string+7));
					}
				}
			} else {
				/* check for non-tab header */
				int ok=0 ; char *c = line->string;
				while (*c) {
					if (*c == '\t') {ok ++; break;}
					c++;
				}
				if (ok >= 8) {
					header = DStringArrayFromCharM(line->string+1,"\t");
				} else {
					header = DStringArrayFromCharM(line->string+1,"\t ");
				}
				break;
			}
		}
	}
	FCLOSE(fh);
	DStringPrintf(dsbuffer,"cg vcfheader2tsv -stack 1 -split %s -typelist '%s' -showheader 0 -meta '%s' '%s'",splitstring,typelist,meta,tempfile);
	system(dsbuffer->string);
	fprintf(fo,"# ----\n");
	fprintf(fo,"chromosome\tbegin\tend\ttype\tref\talt");
	if (!excludename) fprintf(fo,"\tname");
	fprintf(fo,"\tquality");
	if (!excludefilter) fprintf(fo,"\tfilter");
	if (header->size >= 10) {
		samples = DStringArrayRange(header,9,header->size-1);
	} else {
		samples = DStringArrayNew(0);
	}
	if (header->size >= 9) {
		int fieldpos = 0;
		NODPRINT("==== Parsing format ====")
		formatfields = DStringArrayNew(10);
		DStringArrayAppend(formatfields,"GT",2);
		formatfieldsnumber = DStringNew();
		DStringAppendS(formatfieldsnumber,"G",1);
		fieldpos++;
		headerfields = DStringArrayNew(10);
		DStringArrayAppend(headerfields,"alleleSeq1",10);
		DStringArrayAppend(headerfields,"alleleSeq2",10);
		DStringArrayAppend(headerfields,"zyg",3);
		DStringArrayAppend(headerfields,"phased",6);
		DStringArrayAppend(headerfields,"genotypes",9);
		for (i = 0 ; i < format->size ; i ++) {
			DString *ds;
			id = extractID(DStringArrayGet(format,i),id);
			if (strcmp(id->string,"GT") == 0) continue;
			if (strcmp(id->string,"DP") == 0) {DPpos = fieldpos;}
			if (strcmp(id->string,"AD") == 0) {ADpos = fieldpos;}
			if (strcmp(id->string,"MIN_DP") == 0) {MIN_DPpos = fieldpos;}
			if (strcmp(id->string,"GQX") == 0) {GQXpos = fieldpos;}
			if (strcmp(id->string,"GQ") == 0) {GQpos = fieldpos;}
			ds = (DString *)dstring_hash_get(conv_formata,id);
			if (ds == NULL) {ds = id;}
			if (ds->size == 0) {continue;}
			if (keepfields != NULL) {
				DString *kds = (DString *)dstring_hash_get(keepfields,ds);
				if (kds == NULL) {continue;}
			}
			if ((DString *)dstring_hash_get(donefields,ds) != NULL) {
				continue;
			}
			dstring_hash_set(donefields,DStringDup(ds),(void *)"");
			DStringArrayAppend(formatfields,id->string,id->size);
			fieldpos++;
			num->string[0] = extractNumber(DStringArrayGet(format,i));
			num->string[0] = numberfromid(id,num->string[0],typelist);
			DStringAppendS(formatfieldsnumber,num->string,1);
			if (num->string[0] == 'R') {
				DStringSetS(temp,ds->string,ds->size);
				DStringAppendS(temp,"_ref",4);
				DStringArrayAppend(headerfields,temp->string,temp->size);
			}
			DStringArrayAppend(headerfields,ds->string,ds->size);
		}
		NODPRINT("==== Parsing header/samples ====")
		if (samples->size <= 1) {
			for (i = 0 ; i < headerfields->size ; i++) {
				DString *temp = DStringArrayGet(headerfields,i);
				fprintf(fo,"\t%*.*s",temp->size,temp->size,temp->string);
			}
		} else {
			for (j = 0 ; j < samples->size ; j++) {
				for (i = 0 ; i < headerfields->size ; i++) {
					DString *ds = DStringArrayGet(headerfields,i);
					DStringSetS(temp, ds->string, ds->size);
					DStringAppendS(temp, "-", 1);
					DStringAppendS(temp, DStringArrayGet(samples,j)->string,DStringArrayGet(samples,j)->size);
					fprintf(fo,"\t%*.*s",temp->size,temp->size,temp->string);
				}
			}
		}
	}
	NODPRINT("\n\n==== Parsing info ====")
	/* infofields will contain array of info field names */
	infofields = DStringArrayNew(10);
	/* infofieldsnumber will contain array of info field number, = "R","A",".", ... */
	infofieldsnumber = DStringNew();
	infopos = 0;
	doneinfo = hash_init();
	for (i = 0 ; i < info->size ; i ++) {
		DString *ds;
		id = extractID(DStringArrayGet(info,i),id);
		if ((DString *)dstring_hash_get(doneinfo,id) != NULL) {
			continue;
		}
		dstring_hash_set(doneinfo,DStringDup(id),(void *)"");
		if (id->string[0] == 'D' && id->string[1] == 'P' && id->string[2] == '\0') {
			fprintf(fo,"\ttotalcoverage");
		} else {
			ds = (DString *)dstring_hash_get(conv_formata,id);
			if (ds == NULL) {ds = id;}
			if (ds->size == 0) {
				continue;
			} else if (id->size == 3 && strncmp(id->string,"END",3) == 0) {
				svendpos = infopos;
			} else if (id->size == 5 && strncmp(id->string,"SVLEN",5) == 0) {
				svlenpos = infopos;
			} else if (id->size == 6 && strncmp(id->string,"SVTYPE",6) == 0) {
				svtypepos = infopos;
			} else {
				if (id->size == 7 && strncmp(id->string,"STRANDS",7) == 0) {
					strandspos = infopos;
				} else if (id->size == 4 && strncmp(id->string,"CHR2",4) == 0) {
					chr2pos = infopos;
				} else if (id->size == 4 && strncmp(id->string,"POS2",4) == 0) {
					pos2pos = infopos;
				} else if (keepfields != NULL) {
					/* only put keepers in header */
					DString *kds = (DString *)dstring_hash_get(keepfields,ds);
					if (kds == NULL) {continue;}
				}
				num->string[0] = extractNumber(DStringArrayGet(info,i));
				num->string[0] = numberfromid(id,num->string[0],typelist);
				if (dstring_hash_get(donefields,ds) == NULL) {
					if (num->string[0] == 'R') {
						fprintf(fo,"\t%*.*s_ref",ds->size,ds->size,ds->string);
					}
					fprintf(fo,"\t%*.*s",ds->size,ds->size,ds->string);
				} else {
					if (num->string[0] == 'R') {
						fprintf(fo,"\tinfo_%*.*s_ref",ds->size,ds->size,ds->string);
					}
					fprintf(fo,"\tinfo_%*.*s",ds->size,ds->size,ds->string);
				}
			}
		}
		DStringArrayAppend(infofields,id->string,id->size);
		infopos++;
		num->string[0] = extractNumber(DStringArrayGet(info,i));
		num->string[0] = numberfromid(id,num->string[0],typelist);
		DStringAppendS(infofieldsnumber,num->string,1);
	}
	fprintf(fo,"\n");
	if (samples->size) {
		maxtab = 9+samples->size;
		min = 9+samples->size;
	} else {
		maxtab = 8;
		min = 8;
	}
	outinfo =  (DString *)malloc(infofields->size*sizeof(DString));
	linea = DStringArrayNew(maxtab+2);
	NODPRINT("==== Parsing data ====")
	while ((read = DStringGetTab(line,fd,maxtab,linea,1,NULL)) != -1) {
		if (linea->size < min) {
			fprintf(stderr,"not enough fields in line:\n");
			DStringArrayPuts(linea,"\t",stderr);
			fprintf(stderr,"\n");
			exit(1);
		}
		linenr++;
		if (split) {
			int pos;
			if (DStringCompare(a_chrom(linea), prevchr) != 0) {
				/* print out buffer fully on change of chromosome */
				process_print_buffer(fo,obuffer,2147483647);
				DStringCopy(prevchr,a_chrom(linea));
			}
			pos = process_line_split(obuffer,linea,excludename,excludefilter,genotypelist,refout,locerror,skiprefindels);
			process_print_buffer(fo,obuffer,pos-1);
		} else {
			process_line_unsplit(fo,linea,excludename,excludefilter,refout,locerror,skiprefindels);
		}
	}
	process_print_buffer(fo,obuffer,2147483647);
	if (fd != stdin) {FCLOSE(fd);}
	/* free dstrings */
	DStringDestroy(line);
	DStringDestroy(snptype);
	DStringDestroy(instype);
	DStringDestroy(deltype);
	DStringDestroy(subtype);
	DStringDestroy(string);
	DStringDestroy(temp);
	DStringDestroy(alt);
	DStringDestroy(id);
	/* free arrays */
	DStringArrayDestroy(linea);
	DStringArrayDestroy(format);
	DStringArrayDestroy(samples);
	DStringArrayDestroy(formatfields);
	if (headerfields) DStringArrayDestroy(headerfields);
	hash_destroy(conv_formata,(hash_free_func *)DStringDestroy,(hash_free_func *)DStringDestroy);
	/* free other */
	free(order);
	exit(EXIT_SUCCESS);
}
