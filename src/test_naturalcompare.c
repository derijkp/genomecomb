/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include "tools.h"
#include "debug.h"

int main(int argc, char *argv[]) {
	char comp;
	int i,j,diff,alen,blen,errors;
	if ((argc <= 1)) {
		fprintf(stderr,"Format is: test_naturalcompare value1 value2 ...\n");
		exit(EXIT_FAILURE);
	}
	errors = 0;
	for (i = 1 ; i < (argc-1) ; i++) {
		for (j = i+1 ; j < argc ; j++) {
			alen = strlen(argv[i]);
			blen = strlen(argv[j]);
			diff = naturalcompare(argv[i],argv[j],alen,blen);
			if (diff < 0) {
				comp = '<';
			} else if (diff > 0) {
				comp = '>';
				errors++;
			} else {
				comp = '=';
			}
			fprintf(stdout,"%s %c %s (%d)\n",argv[i],comp,argv[j],diff);
			fflush(stdout);
		}
	}
	fprintf(stdout,"nr errors = %d\n",errors);
	exit(EXIT_SUCCESS);
}
