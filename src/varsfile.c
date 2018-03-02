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
#include "gztools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	VarFile *varallvars;
	VariantPos *varpos;
	Variant *var = NULL;
	if (argc != 2) {
		fprintf(stderr,"Format is: varsfile file\n");
		exit(EXIT_FAILURE);
	}
	varallvars = OpenVarfile(argv[1],0);
	varpos = &(varallvars->varpos);
	varfile_next(varallvars,0);
	fprintf(stdout,"chromosome\tbegin\tend\ttype\tref\talt\n");
	while (!varallvars->error) {
		var=varallvars->var;
		if (varpos->chr != -1) {
			DStringputs(var->chr,stdout);
		}
		putc_unlocked('\t',stdout);
		if (varpos->start != -1) {
			fprintf(stdout,"%d",var->start);
		}
		putc_unlocked('\t',stdout);
		if (varpos->end != -1) {
			fprintf(stdout,"%d",var->end);
		}
		putc_unlocked('\t',stdout);
		if (varpos->type != -1) {
			DStringputs(var->type,stdout);
		}
		putc_unlocked('\t',stdout);
		if (varpos->ref != -1) {
			DStringputs(var->ref,stdout);
		}
		putc_unlocked('\t',stdout);
		if (varpos->alt != -1) {
			DStringputs(var->alt,stdout);
		}
		putc_unlocked('\n',stdout);
		varfile_next(varallvars,0);
	}	
	CloseVarfile(varallvars);
	exit(EXIT_SUCCESS);
}
