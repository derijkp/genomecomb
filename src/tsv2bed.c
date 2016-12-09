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
#include "gztools.h"
#include "debug.h"

int findheaderfield(DStringArray *header,const char *string) {
	register int i,len;
	len = strlen(string);
	for (i = 0 ; i < header->size ; i++) {
		if (header->data[i].size == len && strncmp(header->data[i].string,string,len) == 0) {
			return(i);
		}
	}
	return(-1);
}

int main(int argc, char *argv[]) {
	const char* fields[] = {"chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizesblockStarts"};
	VarFile *varfile;
	unsigned int numfields;
	int *poss,numpos,p,pi;
	int max,i;
	if (argc < 2) {
		fprintf(stderr,"Format is: tsv2bed tsvfile ?field?...\n");
		exit(EXIT_FAILURE);
	}
	varfile = OpenVarfile(argv[1],1);
	numpos = argc - 2;
	if (numpos < 3) {numpos = 3;}
	poss = (int *)malloc(numpos*sizeof(int));
	poss[0] = varfile->varpos.chr;
	poss[1] = varfile->varpos.start;
	poss[2] = varfile->varpos.end;
	i = 2;
	pi = 0;
	while (i < argc) {
		const char *string = argv[i];
		if (string[0] == '\0') {
			if (pi < 3 && poss[pi] != -1) {
				i++; pi++;
				continue;
			}
			string = fields[pi];
		}
		p = findheaderfield(varfile->header,string);
		if (p == -1) {
			fprintf(stderr,"Field %s not found\n",string);
			exit(EXIT_FAILURE);
		}
		poss[pi] = p;
		i++; pi++;
	}
	max = 0;
	for(i = 0 ; i < numpos ; i++) {
		if (poss[i] > max) {max = poss[i];}
	}
	varfile->max = max;
	while (1) {
		varfile->error = gz_DStringGetTab(varfile->line,varfile->f,varfile->max,varfile->result,1,&numfields);
		if (varfile->error) break;
		DStringputs(varfile->result->data + poss[0],stdout);
		putc_unlocked('\t',stdout);
		DStringputs(varfile->result->data + poss[1],stdout);
		putc_unlocked('\t',stdout);
		DStringputs(varfile->result->data + poss[2],stdout);
		i = 3;
		while (i < numpos) {
			putc_unlocked('\t',stdout);
			if (poss[i] != -1) {
				DStringputs(varfile->result->data + poss[i],stdout);
			}
			i++;
		}
		putc_unlocked('\n',stdout);
	}
	CloseVarfile(varfile);
	exit(EXIT_SUCCESS);
}
