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
#include "hash.h"

int main(int argc, char *argv[]) {
	FILE *f1;
	Hash_table *hashtable;
	Hash_bucket *bucket;
	Hash_iter iter;
	DStringArray *result1=NULL;
	DString *line1 = NULL;
	DString *key = NULL;
	unsigned int numfields1, pos, pos1;
	long int temp;
	int new;
	if ((argc != 2) && (argc != 3)) {
		fprintf(stderr,"Format is: tsv_histo histopos ?header?");
		exit(EXIT_FAILURE);
	}
	f1 = stdin;
	pos = atoi(argv[1]);
	hashtable = hash_init_size(500);
	line1 = DStringNew();
	/* The following allocation is not destroyed at end as it may point to something else */
	/* This will leak mem, but as the prog is finished anyway ... */
	result1 = DStringArrayNew(pos+2);
	if (argc == 3) {
		int header = atoi(argv[2]);
		if (header) {
			skip_header(f1,line1,&numfields1,&pos1);
		}
	}
	while (!DStringGetTab(line1,f1,pos,result1,1,&numfields1)) {
		key = result1->data+pos;
		bucket = hash_get(hashtable, (void *)(key), hash_DString_hash, hash_DString_compare, &new,1);
		if (new == 0) {
			/* key was already present in the hashtable */
			temp = (long int)hash_getvalue(bucket);
			temp++;
			hash_setvalue(bucket,temp);
		} else {
			DString *buffer = DStringNew();
			DStringSetS(buffer,key->string,key->size);
			hash_setkey(bucket,buffer);
			hash_setvalue(bucket,1);
			buffer = NULL;
		}
	}
	fclose(f1);
	bucket = hash_first(hashtable,&iter);
	while(bucket != NULL) {
		DString *ds = hash_getkey(bucket);
		temp = (long int)hash_getvalue(bucket);
		fprintf(stdout,"%s\t%ld\n",ds->string,temp);
		DStringDestroy(ds);
		bucket = hash_next(&iter);
	}
	hash_destroy(hashtable,NULL,NULL);
	if (line1) {DStringDestroy(line1);}
	if (result1) {DStringArrayDestroy(result1);}
	exit(EXIT_SUCCESS);
}
