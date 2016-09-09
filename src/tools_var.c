/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "debug.h"

typedef struct VariantPos {
	int chr;
	int start;
	int end;
	int type;
	int ref;
	int alt;
	int seq;
	int zyg;
	int a1;
	int a2;
	int max;
	int id;
} VariantPos;

typedef struct Variant {
	DString *chr;
	int start;
	int end;
	DString *type;
	DString *ref;
	DString *alt;
	int id;
} Variant;

/* returns 0 if equal, 1 if equal not including alt, 2 if different */
int varchecksort(Variant *prev,Variant *var,char *filename,int *nextpos) {
	return checksort(prev->chr,&(prev->start),&(prev->end),prev->type,prev->alt,var->chr,var->start,var->end,var->type,var->alt,filename,nextpos);
}

void varputs(Variant var,FILE *f) {
	if (var.chr->size > 3) {
		char *a = var.chr->string;
		if ((a[0] == 'C' || a[0] == 'c') && (a[1] == 'H' || a[1] == 'h') && (a[2] == 'R' || a[2] == 'r')) {
			charputs(a+3,var.chr->size-3,f);
		} else {
			DStringputs(var.chr,f);
		}
	} else {
		DStringputs(var.chr,f);
	}
	fprintf(f,"\t%d\t%d\t",var.start,var.end);
	DStringputs(var.type,f);
	putc_unlocked('\t',f);
	DStringputs(var.ref,f);
	putc_unlocked('\t',f);
	DStringputs(var.alt,f);
	if (var.id != -2) {
		if (var.id == -1) {
			fprintf(f,"\t");
		} else {
			fprintf(f,"\t%d",var.id);
		}
	}
	putc_unlocked('\n',f);
}

void result2var(DStringArray *result,VariantPos varpos, Variant *var) {
	var->chr = result->data+varpos.chr;
	sscanf(result->data[varpos.start].string,"%d",&(var->start));
	sscanf(result->data[varpos.end].string,"%d",&(var->end));
	var->type = result->data+varpos.type;
	var->ref = result->data+varpos.ref;
	if (varpos.alt != -1) {
		var->alt = result->data+varpos.alt;
	}
	if (varpos.id == -1) {
		var->id = -1;
	} else if (varpos.id == -2) {
		var->id = -2;
	} else {
		if (result->data[varpos.id].size != 0) {
			sscanf(result->data[varpos.id].string,"%d",&(var->id));
		} else {
			var->id = -1;
		}
	}
}

int varcompare(Variant *var1, Variant *var2, int split) {
	int comp;
	comp = DStringLocCompare(var1->chr,var2->chr);
	if  (comp != 0) {return comp<0?-2:2;}
	comp = var1->start - var2->start;
	if (comp != 0) {return comp<0?-2:2;}
	comp = var1->end - var2->end;
	if (comp != 0) {return comp<0?-2:2;}
	comp = DStringCompare(var1->type,var2->type);
	if (comp != 0) {return comp<0?-2:2;}
	if (split) {
		comp = DStringCompare(var1->alt,var2->alt);
		if (comp != 0) {return comp<0?-1:1;}
	}
	return 0;
}

void varpos_init(VariantPos *varpos) {
	varpos->max = 0;
	varpos->start = -1;
	varpos->end = -1;
	varpos->type = -1;
	varpos->ref = -1;
	varpos->alt = -1;
	varpos->seq = -1;
	varpos->zyg = -1;
	varpos->a1 = -1;
	varpos->a2 = -1;
	varpos->id = -1;
}

int varpos_max(VariantPos *varpos) {
	int i;
	i = varpos->chr;
	if (varpos->start > i) {i = varpos->start;} ;
	if (varpos->end > i) {i = varpos->end;} ;
	if (varpos->type > i) {i = varpos->type;} ;
	if (varpos->ref > i) {i = varpos->ref;} ; 
	if (varpos->alt > i) {i = varpos->alt;} ;
	if (varpos->seq > i) {i = varpos->seq;} ;
	if (varpos->zyg > i) {i = varpos->zyg;} ;
	if (varpos->a1 > i) {i = varpos->a1;} ;
	if (varpos->a2 > i) {i = varpos->a2;} ;
	if (varpos->id > i) {i = varpos->id;} ;
	varpos->max = i;
	return i;
}

/*
int varpos_findheaderfield_string(char *header,int *lengths,int headersize,const char *array[]) {
	DString *teststring;
	char *string;
	register char *test;
	register int i,len,j,alen;
	alen = (sizeof (array) / sizeof (const char *));
	for (j = 0 ; j < alen ; j++) {
		test = header;
		string = array[j]
		len = strlen(string);
		for (i = 0 ; i < headersize ; i++) {
			if (len == lengtsa[i] && strncmp(test,string,len)) {
				return(i);
			}
			test += lengtsa[i] + 1;
		}
	}
	return(-1);
}
*/

int varpos_findheaderfield(DStringArray *header,const char *array[],int alen) {
	const char *string;
	register int i,len,j;
	for (j = 0 ; j < alen ; j++) {
		string = array[j];
		len = strlen(string);
		for (i = 0 ; i < header->size ; i++) {
			if (header->data[i].size == len && strncmp(header->data[i].string,string,len) == 0) {
				return(i);
			}
		}
	}
	return(-1);
}

void varpos_fromheader(VariantPos *varpos,DStringArray *header) {
	const char * chra[] = {"chromosome","chrom","chr","chr1","genoName","contig"};
	const char * begina[] = {"begin","start","end1","chromStart","genoStart","tStart","txStart"};
	const char * enda[] = {"end","start2","chromEnd","genoEnd","tEnd","txEnd"};
	const char * typea[] = {"type"};
	const char * refa[] = {"ref","reference"};
	const char * alta[] = {"alt","alternative"};
	const char * seqa[] = {"sequenced"};
	const char * zyga[] = {"zyg"};
	const char * a1a[] = {"alleleSeq1"};
	const char * a2a[] = {"alleleSeq2"};
	varpos_init(varpos);
	varpos->chr = varpos_findheaderfield(header,chra,(sizeof (chra) / sizeof (const char *)));
	varpos->start = varpos_findheaderfield(header,begina,(sizeof (begina) / sizeof (const char *)));
	varpos->end = varpos_findheaderfield(header,enda,(sizeof (enda) / sizeof (const char *)));
	varpos->type = varpos_findheaderfield(header,typea,(sizeof (typea) / sizeof (const char *)));
	varpos->ref = varpos_findheaderfield(header,refa,(sizeof (refa) / sizeof (const char *)));
	varpos->alt = varpos_findheaderfield(header,alta,(sizeof (alta) / sizeof (const char *)));
	varpos->seq = varpos_findheaderfield(header,seqa,(sizeof (seqa) / sizeof (const char *)));
	varpos->zyg = varpos_findheaderfield(header,zyga,(sizeof (zyga) / sizeof (const char *)));
	varpos->a1 = varpos_findheaderfield(header,a1a,(sizeof (a1a) / sizeof (const char *)));
	varpos->a2 = varpos_findheaderfield(header,a2a,(sizeof (a2a) / sizeof (const char *)));
}

