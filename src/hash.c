#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "tools.h"
#include "hash.h"
#include "debug.h"

Hash_table *hash_init_size(int tablesize) {
	Hash_table *table;
	size_t size = tablesize*sizeof(Hash_tableitem);
	table = (Hash_table *)malloc(sizeof(Hash_table));
	table->table = (Hash_tableitem *)malloc(size);
	memset(table->table,0,size);
	table->max = tablesize-1;
	table->maxload = (int)((double)size*0.7);
	table->datasize = 0;
	table->items = NULL;
	return(table);
}

void hash_resize(Hash_table *table, hash_hash_func *hashfunc) {
	Hash_table *newtable;
	Hash_tableitem *olditem,*tableitem;
	Hash_bucket *oldbucket, *bucket;
	newtable = hash_init_size(2*(table->max+1));
	/* redo all stored buckets */
	olditem = table->items;
	while (olditem != NULL) {
		oldbucket = olditem->chain;
		while (oldbucket->key != NULL) {
			unsigned int hash;
			hash = hashfunc(oldbucket->key,newtable->max);
			/* can we find the bucket with the right key	*/
			/* will allways be new buckets, so do not check for existing here, no compar function needed */
			tableitem = newtable->table+hash;
			if (tableitem->chain == NULL) {
				/* not found, create new chain */
				bucket = (Hash_bucket *)malloc(2*sizeof(Hash_bucket));
				NODPRINT("new bucket")
				tableitem->chain = bucket;
				tableitem->next = newtable->items;
				newtable->items = tableitem;
			} else {
				int pos = 0;
				bucket = tableitem->chain;
				while(1) {
					bucket++;
					pos++;
					if (bucket->key == NULL) break;
				}
				NODPRINT("bucket at %d",pos)
				tableitem->chain = (Hash_bucket *)realloc(tableitem->chain,(pos+2)*sizeof(Hash_bucket));
				bucket = tableitem->chain + pos;
			}
			bucket->key = oldbucket->key;
			bucket->value = oldbucket->value;
			bucket[1].key = NULL;
			bucket[1].value = NULL;
			oldbucket++;
		}
		free(olditem->chain);
		olditem = olditem->next;
	}
	free(table->table);
	table->table = newtable->table;
	table->items = newtable->items;
	table->max = newtable->max;
	table->datasize = newtable->datasize;
	table->maxload = newtable->maxload;
	free(newtable);
}

Hash_bucket *hash_get(Hash_table *table, void *key, hash_hash_func *hashfunc, hash_key_compare_func *compare, int *new) {
	Hash_tableitem *tableitem;
	Hash_bucket *bucket;
	unsigned int hash;
	int comp;
	hash = hashfunc(key,table->max);
	NODPRINT("hash %s -> %d",((DString *)key)->string,hash);
	/* can we find the bucket with the right key	*/
	tableitem = table->table+hash;
	if (tableitem->chain == NULL) {
		/* not found, create new chain */
		bucket = (Hash_bucket *)malloc(2*sizeof(Hash_bucket));
		NODPRINT("new bucket")
		tableitem->chain = bucket;
		tableitem->next = table->items;
		table->items = tableitem;
	} else {
		int pos = 0;
		bucket = tableitem->chain;
		while(1) {
			comp = compare(key,bucket->key);
			if (comp == 0) {
				*new = 0;
				return bucket;
			}
			bucket++;
			pos++;
			if (bucket->key == NULL) break;
		}
		NODPRINT("bucket at %d",pos)
		tableitem->chain = (Hash_bucket *)realloc(tableitem->chain,(pos+2)*sizeof(Hash_bucket));
		bucket = tableitem->chain + pos;
	}
	table->datasize ++;
	bucket->key = key;
	bucket->value = NULL;
	bucket[1].key = NULL;
	bucket[1].value = NULL;
	*new = 1;
	/* do we need a resize */
	if (table->datasize >= table->maxload) {
		hash_resize(table,hashfunc);
	}
	return bucket;
}

Hash_bucket *hash_first(Hash_table *table,Hash_iter *iter) {
	iter->table = table;
	iter->item = table->items;
	if (table->items == NULL) {
		iter->bucket = NULL;
	} else {
		iter->bucket = iter->item->chain;
	}
	return(iter->bucket);
}

Hash_bucket *hash_next(Hash_iter *iter) {
	Hash_bucket *bucket;
	if (iter->item == NULL) {
		return NULL;
	}
	bucket = ++(iter->bucket);
	if (bucket->key != NULL) {
		return bucket;
	}
	iter->item = iter->item->next;
	if (iter->item == NULL) {
		return NULL;
	}
	iter->bucket = iter->item->chain;
	return iter->bucket;
}

void hash_destroy(Hash_table *table) {
	Hash_tableitem *item = table->items;
	while(1) {
		if (item == NULL) break;
		free((char *)item->chain);
		item = item->next;
	}
	free(table);
}

unsigned int hash_Dstring_hash (void *key,unsigned int hashtablemax)
{
	register unsigned int result;
	register int c;
	DString *ds = key;
	char *string = ds->string;
	int count = ds->size;
	if (count <= 0) {
		result = 0;
	} else {
		result = 0;
		while (--count) {
			c = *string;
			result += (result<<3) + c;
			string++;
		}
	}
	return result & hashtablemax;
}

int hash_Dstring_compare(const void* key1, const void* key2)
{
	DString *ds1 = (DString *)key1, *ds2 = (DString *)key2;
	return DStringCompare(ds1,ds2);
}
