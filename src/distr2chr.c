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
#include "khash-dstring.h"

KHASH_MAP_INIT_DSTR(DSTRING, FILE*);

/*#include "dstring-khash.h"*/

int main(int argc, char *argv[]) {
	khash_t(DSTRING) *hashtable;
	FILE *o;
	DString *result = NULL, *buffer = NULL;
	DString *line = NULL,*chromosome = NULL;
	khiter_t k;
	int col = 0,header = 0,ret;
	if ((argc < 2)&&(argc > 4)) {
		fprintf(stderr,"Format is: distr2chr output_pre ?col? ?header?");
		exit(EXIT_FAILURE);
	}
	if (argc >= 3) {
		col = atoi(argv[2]);
	}
	line = DStringNew();
	if (argc == 4) {
		header = atoi(argv[2]);
		skip_header(stdin,line);
	}
	result = DStringArrayNew(col+1);
	hashtable = kh_init(DSTRING);
	while (!DStringGetTab(line,stdin,col,result,0)) {
		chromosome = result+col;
		k = kh_put(DSTRING,hashtable, chromosome, &ret);
		if (ret == 0) {
			/* key was already present in the hashtable */
			o = kh_value(hashtable, k);
		} else {
			buffer = DStringNew();
			DStringAppend(buffer,argv[1]);
			DStringAppendS(buffer,chromosome->string,chromosome->size);
			o = fopen64(buffer->string,"w");
			kh_value(hashtable, k) = o;
			DStringSetS(buffer,chromosome->string,chromosome->size);
			kh_key(hashtable, k) = buffer;
			buffer = NULL;
		}
		fprintf(o,"%s\n", line->string);
	}
	for (k = kh_begin(hashtable); k != kh_end(hashtable); ++k) {
		if (kh_exist(hashtable, k)) {
			DStringDestroy(kh_key(hashtable, k));
			fclose(kh_value(hashtable, k));
		}
	}
	if (line) {DStringDestroy(line);}
	if (result) {DStringArrayDestroy(result);}
	kh_destroy(DSTRING,hashtable);
	exit(EXIT_SUCCESS);
}
