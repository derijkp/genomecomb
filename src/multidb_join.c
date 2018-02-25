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

#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct Todo {
	DString *filename;
	FILE *f;
	int max;
	DString *line;
	DStringArray *result;
	unsigned int pos;
	unsigned int numfields;
	int error;
	VariantPos varpos;
	Variant var;
	int seqpos;
	int keepsize;
	int *keepposs;
} Todo;

int main(int argc, char *argv[]) {
	FILE *ftodo,*fvarsnew,*fvarsinsert,*fgeno,*f;
	char *fvarsnewfile,*fvarsinsertfile,*fgenofile;
	DStringArray *result1=NULL;
	DString *filename;
	Todo *todolist = NULL,*todo,*tvars;
	Variant *var;
	DString *line = NULL;
	char *todofile;
	char varid[20];
	unsigned int numfields,pos1,pos2;
	int fvarsnewcount = 0,fvarsinsertcount = 0,fgenocount = 0;
	int split = 1,first;
	int comp, size, newvarid;
	register int i,j;

	if ((argc != 8)) {
		fprintf(stderr,"Format is: multi_join todofile number split varsfile.new varsfile.insert genofile.insert newvaridstart\n");
		exit(EXIT_FAILURE);
	}
	todofile = argv[1];
	size = atoi(argv[2]);
	split = atoi(argv[3]);
	fvarsnewfile = argv[4];
	fvarsinsertfile = argv[5];
	fgenofile = argv[6];
	newvarid = atoi(argv[7]);
	ftodo = fopen64_or_die(todofile,"r");
	fvarsnew = fopen64_or_die(fvarsnewfile,"a");
	fvarsinsert = fopen64_or_die(fvarsinsertfile,"a");
	fgeno = fopen64_or_die(fgenofile,"a");
	todolist = (Todo *)malloc(size*sizeof(Todo));
	tvars = todolist;
	todo = todolist;
	line = DStringNew();
	/* The following allocation is not destroyed at end as it may point to something else */
	/* This will leak mem, but as the prog is finished anyway ... */
	result1 = DStringArrayNew(13);
	for (i = 0 ; i < size; i++) {
		DStringGetTab(line,ftodo,11,result1,1,&numfields);
		if (numfields < 11) {
			fprintf(stderr,"only %d values on line %d of todofile, should be 11\n", numfields, i+1);
			exit(EXIT_FAILURE);
		}
		todo->filename = DStringNew();
		DStringSetS(todo->filename,result1->data[0].string,result1->data[0].size);
		todo->f = fopen64_or_die(result1->data[0].string,"r");
		varpos_init(&(todo->varpos));
		todo->varpos.chr = atoi(result1->data[1].string);
		todo->varpos.start = atoi(result1->data[2].string);
		todo->varpos.end = atoi(result1->data[3].string);
		todo->varpos.type = atoi(result1->data[4].string);
		todo->varpos.ref = atoi(result1->data[5].string);
		todo->varpos.alt = atoi(result1->data[6].string);
		todo->varpos.id = atoi(result1->data[7].string);
		varpos_max(&(todo->varpos));
		todo->max = atoi(result1->data[8].string);
		todo->seqpos = atoi(result1->data[9].string);
		/* keepsize contains the size of the next line, the number of positions to keep in the output */
		todo->keepsize = atoi(result1->data[10].string);
		/* get keepposs, allocate todo->result first and use it also to get poss here (allocate large enough for both) */
		todo->line = DStringNew();
		todo->result = DStringArrayNew(MAX(todo->max+2,todo->keepsize+3));
		DStringGetTab(line,ftodo,todo->keepsize+1,todo->result,0,&numfields);
		todo->keepposs = (int *)malloc(todo->keepsize*sizeof(int));
		for (j = 0 ; j < todo->keepsize ; j++) {
			todo->keepposs[j] = atoi(todo->result->data[j].string);
		}
		skip_header(todo->f,todo->line,&(todo->numfields),&(todo->pos));
		todo->error = DStringGetTab(todo->line,todo->f,todo->max,todo->result,1,&numfields); pos2++;
		if (!todo->error) {
			check_numfieldserror(numfields,todo->numfields,todo->line,todo->filename->string,&(todo->pos));
			result2var(todo->result,todo->varpos,&(todo->var));
		}
		todo++;
	}
	fclose(ftodo);
	var = &(tvars->var);
	while (!tvars->error) {
		pos1++;
		DStringputs(var->chr,fvarsnew);
		fprintf(fvarsnew,"\t%d\t%d\t",var->start,var->end);
		DStringputs(var->type,fvarsnew);
		putc_unlocked('\t',fvarsnew);
		DStringputs(var->ref,fvarsnew);
		putc_unlocked('\t',fvarsnew);
		DStringputs(var->alt,fvarsnew);
		if (var->id < 0) {
			/* this is a new var, also write in insert */
			var->id = newvarid++;
			DStringputs(var->chr,fvarsinsert);
			fprintf(fvarsinsert,"\t%d\t%d\t",var->start,var->end);
			DStringputs(var->type,fvarsinsert);
			putc_unlocked('\t',fvarsinsert);
			DStringputs(var->ref,fvarsinsert);
			putc_unlocked('\t',fvarsinsert);
			DStringputs(var->alt,fvarsinsert);
			fprintf(fvarsinsert,"\t%d\n",var->id);
			fvarsinsertcount++;
		}
		sprintf(varid,"%d",var->id);
		fprintf(fvarsnew,"\t%s\n",varid);
		fvarsnewcount++;
		todo = todolist+1;
		for (i = 1 ; i < size; i++) {
			if (todo->error) {
				comp = -1;
			} else {
				comp = varcompare(var,&(todo->var),split);
			}
			if (comp > 0) {
				fprintf(stderr,"all variants should be in the variant file, so this should not happen");
				exit(EXIT_FAILURE);
			} else if (comp == 0) {
				/* variant is present in the file */
				register int *cur = todo->keepposs;
				j = todo->keepsize;
				first = 1;
				while (j--) {
					if (*cur == -1) {
						putc_unlocked('\t',fgeno);
					} else if (*cur == -2) {
						if (!first) {
							putc_unlocked('\n',fgeno);
							fgenocount++;
						}
						first = 0;
						j--; cur++;
						fprintf(fgeno,"%s\t%d",varid,*cur);
						if (todo->seqpos == -1) {
							putc_unlocked('\t',fgeno);
							putc_unlocked('v',fgeno);
						}
					} else {
						putc_unlocked('\t',fgeno);
						DStringputs(todo->result->data+*cur,fgeno);
					}
					cur++;
				}
				putc_unlocked('\n',fgeno);
				fgenocount++;
				todo->error = DStringGetTab(todo->line,todo->f,todo->max,todo->result,1,&numfields);
				if (!todo->error) {
					check_numfieldserror(numfields,todo->numfields,todo->line,todo->filename->string,&(todo->pos));
					result2var(todo->result,todo->varpos,&(todo->var));
				}
			} else {
				/* variant is not present in the file, do not write anything to the geno file */
			}
			todo++;
		}
		tvars->error = DStringGetTab(tvars->line,tvars->f,tvars->varpos.max,tvars->result,1,&numfields);
		if (tvars->error) break;
		check_numfieldserror(numfields,tvars->numfields,tvars->line,tvars->filename->string,&(tvars->pos));
		result2var(tvars->result,tvars->varpos,var);
	}
	todo = todolist;
	for (i = 0 ; i < size; i++) {
		fclose(todo->f);
		todo++;
	}
	fclose(fvarsnew);
	fclose(fvarsinsert);
	fclose(fgeno);
	filename = DStringNew();
	#
	DStringSet(filename,fvarsnewfile);
	DStringAppend(filename,".maxid");
	f = fopen64_or_die(filename->string,"w");
	fprintf(f,"%s",varid);
	fclose(f);
	#
	DStringSet(filename,fvarsnewfile);
	DStringAppend(filename,".count");
	f = fopen64_or_die(filename->string,"w");
	fprintf(f,"%d",fvarsnewcount);
	fclose(f);
	#
	DStringSet(filename,fvarsinsertfile);
	DStringAppend(filename,".count");
	f = fopen64_or_die(filename->string,"w");
	fprintf(f,"%d",fvarsinsertcount);
	fclose(f);
	#
	DStringSet(filename,fgenofile);
	DStringAppend(filename,".count");
	f = fopen64_or_die(filename->string,"w");
	fprintf(f,"%d",fgenocount);
	fclose(f);
	#
	exit(EXIT_SUCCESS);
}
