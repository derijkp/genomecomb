/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "tools_var.h"
#include "debug.h"

/* returns 0 if equal, 1 if equal not including alt, 2 if different */
int varchecksort(Variant *prev,Variant *var,char *filename,int *nextpos) {
	return checksort(prev->chr,&(prev->start),&(prev->end),prev->type,prev->alt,var->chr,var->start,var->end,var->type,var->alt,filename,nextpos,0);
}

void varputs_chr(DString *chr,FILE *f) {
	if (chr->size > 3) {
		char *a = chr->string;
		if ((a[0] == 'C' || a[0] == 'c') && (a[1] == 'H' || a[1] == 'h') && (a[2] == 'R' || a[2] == 'r')) {
			charputs(a+3,chr->size-3,f);
		} else {
			DStringputs(chr,f);
		}
	} else {
		DStringputs(chr,f);
	}
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
	if (varpos.type != -1) {var->type = result->data+varpos.type;} else {var->type = DStringEmtpy();}
	if (varpos.ref != -1) {var->ref = result->data+varpos.ref;} else {var->ref = DStringEmtpy();}
	if (varpos.a1 != -1) {var->a1 = result->data+varpos.a1;} else {var->a1 = NULL;}
	if (varpos.a2 != -1) {var->a2 = result->data+varpos.a2;} else {var->a2 = NULL;}
	if (varpos.alt != -1) {var->alt = result->data+varpos.alt;} else {var->alt = DStringEmtpy();}
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

int regcompare(Variant *var1, Variant *var2) {
	int comp;
	comp = DStringLocCompare(var1->chr,var2->chr);
	if  (comp != 0) {return comp<0?-2:2;}
	comp = var1->start - var2->start;
	if (comp != 0) {return comp<0?-2:2;}
	comp = var1->end - var2->end;
	if (comp != 0) {return comp<0?-2:2;}
	return 0;
}

void var_init(Variant *var) {
	var->start = -1;
	var->end = -1;
	var->type = NULL;
	var->ref = NULL;
	var->alt = NULL;
	var->a1 = NULL;
	var->a2 = NULL;
	var->id = -1;
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

VarFile *OpenVarfile(char *filename,int split) {
	VarFile *varfile = malloc(sizeof(VarFile));
	varfile->linenr = 0;
	varfile->file = filename;
	varfile->split = split;
	varfile->f = gz_open(filename);
	varfile->prevline = NULL;
	varfile->line = DStringNew();
	varfile->headerline = DStringNew();
	varfile->prevline = DStringNew();
	gz_skip_header(varfile->f,varfile->headerline,&varfile->numfields,&varfile->pos);
	varfile->header = DStringArrayNew(varfile->numfields+2);
	varfile->result = DStringArrayNew(varfile->numfields+2);
	varfile->prevresult = DStringArrayNew(varfile->numfields+2);
	DStringSplitTab(varfile->headerline, varfile->numfields, varfile->header, 0,NULL);
	varpos_fromheader(&varfile->varpos,varfile->header);
	varpos_max(&(varfile->varpos));
	varfile->max = varfile->varpos.max;
	varfile->var = (Variant *)malloc(sizeof(Variant));
	var_init(varfile->var);
	varfile->prevvar = (Variant *)malloc(sizeof(Variant));
	var_init(varfile->prevvar);
	varfile->error = 0;
	return varfile;
}

void Varfile_checkbasicfields(VarFile *varfile) {
	VariantPos *varpos;
	varpos = &(varfile->varpos);
	if (varpos->chr == -1) {
		fprintf(stderr,"field chromosome not found in %s\n",varfile->file); exit(EXIT_FAILURE);
	}
	if (varpos->start == -1) {
		fprintf(stderr,"field begin not found in %s\n",varfile->file); exit(EXIT_FAILURE);
	}
	if (varpos->end == -1) {
		fprintf(stderr,"field end not found in %s\n",varfile->file); exit(EXIT_FAILURE);
	}
	if (varpos->type == -1) {
		fprintf(stderr,"field type not found in %s\n",varfile->file); exit(EXIT_FAILURE);
	}
	if (varpos->ref == -1) {
		fprintf(stderr,"field ref not found in %s\n",varfile->file); exit(EXIT_FAILURE);
	}
	if (varpos->alt == -1) {
		fprintf(stderr,"field alt not found in %s\n",varfile->file); exit(EXIT_FAILURE);
	}
}

void Varfile_checkbasicfields_basic(VarFile *varfile) {
	VariantPos *varpos;
	varpos = &(varfile->varpos);
	if (varpos->chr == -1) {
		fprintf(stderr,"field chromosome not found in %s\n",varfile->file); exit(EXIT_FAILURE);
	}
	if (varpos->start == -1) {
		fprintf(stderr,"field begin not found in %s\n",varfile->file); exit(EXIT_FAILURE);
	}
	if (varpos->end == -1) {
		fprintf(stderr,"field end not found in %s\n",varfile->file); exit(EXIT_FAILURE);
	}
}

Variant *varfile_next(VarFile *varfile,int check) {
	Variant *tempvar;
	DString *templine;
	DStringArray *tempresult;
	unsigned int numfields;
	int compcheck = 0;
	/* keep previous var
           =================
	   as var points to line and result, we need to keep these as wel
	   in order to avoid many allocations and copying, we recycle them
	*/
	if (varfile->error) {return NULL;}
	templine = varfile->prevline;
	varfile->prevline = varfile->line;
	varfile->line = templine;
	tempvar = varfile->prevvar;
	varfile->prevvar = varfile->var;
	varfile->var = tempvar;
	tempresult = varfile->prevresult;
	varfile->prevresult = varfile->result;
	varfile->result = tempresult;
	/* get next line
	   ============= */
	varfile->error = gz_DStringGetTab(varfile->line,varfile->f,varfile->max,varfile->result,1,&numfields);
	if (varfile->error) {return NULL;}
	check_numfieldserror(numfields,varfile->numfields,varfile->line,varfile->file,&varfile->pos);
	result2var(varfile->result,varfile->varpos,varfile->var);
	varfile->linenr ++;
	if (check && varfile->linenr > 1) {
		compcheck = varchecksort(varfile->prevvar,varfile->var,varfile->file,NULL);
		if (!varfile->split && compcheck < 2) {
			Variant *prev = varfile->prevvar;
			fprintf(stderr,"error in \"%s\": file uses split alleles (\"%s %d %d %s\" occurs more than once and you are not running with the -split option)",
				varfile->file,prev->chr->string,prev->start,prev->end,prev->type->string);
			exit(EXIT_FAILURE);
		}
	}
	return varfile->var;
}

void CloseVarfile(VarFile *varfile) {
	gz_close(varfile->f);
	DStringDestroy(varfile->line);
	DStringDestroy(varfile->prevline);
	DStringArrayDestroy(varfile->result);
	DStringArrayDestroy(varfile->prevresult);
	free(varfile->var);
	free(varfile->prevvar);
	free(varfile);
}

