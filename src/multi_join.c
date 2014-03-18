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
#include "tools_var.h"
#include "debug.h"

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
	FILE *f;
	DStringArray *result1=NULL;
	Todo *todolist = NULL,*todo,*tvars;
	Variant *var;
	DString *line = NULL;
	char *todofile;
	unsigned int numfields,pos1,pos2;
	int split = 1;
	int comp, size;
	register int i,j;

	if ((argc != 4)) {
		fprintf(stderr,"Format is: multi_join todofile number split\n");
		exit(EXIT_FAILURE);
	}
	todofile = argv[1];
	f = fopen64_or_die(todofile,"r");
	size = atoi(argv[2]);
	split = atoi(argv[3]);
	todolist = (Todo *)malloc(size*sizeof(Todo));
	tvars = todolist;
	todo = todolist;
	line = DStringNew();
	/* The following allocation is not destroyed at end as it may point to something else */
	/* This will leak mem, but as the prog is finished anyway ... */
	result1 = DStringArrayNew(12);
	for (i = 0 ; i < size; i++) {
		DStringGetTab(line,f,10,result1,1,&numfields);
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
		varpos_max(&(todo->varpos));
		todo->max = atoi(result1->data[7].string);
		todo->seqpos = atoi(result1->data[8].string);
		/* keepsize contains the size of the next line, the number of positions to keep in the output */
		todo->keepsize = atoi(result1->data[9].string);
		/* het keepposs, allocate todo->result first and use to get poss: it will be large enough; max must be > keepsize */
		todo->line = DStringNew();
		todo->result = DStringArrayNew(todo->max+2);
		DStringGetTab(line,f,todo->keepsize+1,todo->result,0,&numfields);
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
	fclose(f);
	var = &(tvars->var);
	while (!tvars->error) {
		pos1++;
		DStringputs(var->chr,stdout);
		fprintf(stdout,"\t%d\t%d\t",var->start,var->end);
		DStringputs(var->type,stdout);
		putc_unlocked('\t',stdout);
		DStringputs(var->ref,stdout);
		putc_unlocked('\t',stdout);
		DStringputs(var->alt,stdout);
		todo = todolist+1;
		for (i = 1 ; i < size; i++) {
			if (todo->error) {
				comp = -1;
			} else {
				comp = varcompare(var,&(todo->var),split);
			}
			if (comp > 0) {
				fprintf(stderr,"This should not happen");
				exit(EXIT_FAILURE);
			} else if (comp == 0) {
				register int *cur = todo->keepposs;
				if (todo->seqpos == -1) {
					putc_unlocked('\t',stdout);
					putc_unlocked('v',stdout);
				}
				j = todo->keepsize;
				while (j--) {
					putc_unlocked('\t',stdout);
					DStringputs(todo->result->data+*cur,stdout);
					cur++;
				}
				todo->error = DStringGetTab(todo->line,todo->f,todo->max,todo->result,1,&numfields);
				if (!todo->error) {
					check_numfieldserror(numfields,todo->numfields,todo->line,todo->filename->string,&(todo->pos));
					result2var(todo->result,todo->varpos,&(todo->var));
				}
			} else {
				if (todo->seqpos == -1) {
					putc_unlocked('\t',stdout);
					putc_unlocked('?',stdout);
				}
				j = todo->keepsize;
				while (j--) {
					putc_unlocked('\t',stdout);
					putc_unlocked('?',stdout);
				}
			}
			todo++;
		}
		putc_unlocked('\n',stdout);
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
	exit(EXIT_SUCCESS);
}
