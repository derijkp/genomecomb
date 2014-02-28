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
	int a1;
	int a2;
	int max;
} VariantPos;

typedef struct Variant {
	DString *chr;
	int start;
	int end;
	DString *type;
	DString *ref;
	DString *alt;
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
	} else {
	}
}

int varcompare(Variant *var1, Variant *var2, int split) {
	int comp;
	comp = DStringLocCompare(var1->chr,var2->chr);
	if  (comp != 0) {return comp;}
	comp = var1->start - var2->start;
	if (comp != 0) {return comp;}
	comp = var1->end - var2->end;
	if (comp != 0) {return comp;}
	comp = DStringCompare(var1->type,var2->type);
	if (comp != 0) {return comp;}
	if (split) {
		comp = DStringCompare(var1->alt,var2->alt);
		if (comp != 0) {return comp;}
	}
	return 0;
}

int varpos_max(VariantPos *varpos) {
	int i;
	varpos->max = varpos->chr;
	if (varpos->start > i) {i = varpos->start;} ;
	if (varpos->end > i) {i = varpos->end;} ;
	if (varpos->type > i) {i = varpos->type;} ;
	if (varpos->ref > i) {i = varpos->ref;} ; 
	if (varpos->alt > i) {i = varpos->alt;} ;
	if (varpos->a1 > i) {i = varpos->a1;} ;
	if (varpos->a2 > i) {i = varpos->a2;} ;
	varpos->max = i;
	return i;
}
