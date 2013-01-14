typedef struct Hash_bucket {
	void *key;
	void *value;
} Hash_bucket;

typedef struct Hash_tableitem {
	Hash_bucket *chain;
	struct Hash_tableitem *next;
} Hash_tableitem;

typedef struct Hash_table {
	Hash_tableitem *table;
	Hash_tableitem *items;
	int max;
	int datasize;
	int maxload;
} Hash_table;

typedef struct Hash_iter {
	Hash_table *table;
	Hash_tableitem *item;
	Hash_bucket *bucket;
} Hash_iter;

typedef int hash_key_compare_func(const void* key1, const void* key2);
typedef int hash_hash_func(void* key, int hashtablesize);

Hash_table *hash_init_size(int tablesize);
void hash_destroy(Hash_table *table);
Hash_bucket *hash_get(Hash_table *table, void *key, hash_hash_func *hashfunc, hash_key_compare_func *compare, int *new);
Hash_bucket *hash_first(Hash_table *table,Hash_iter *bucket);
Hash_bucket *hash_next(Hash_iter *iter);

#define hash_getkey(bucket) (bucket->key);
#define hash_getvalue(bucket) (bucket->value);
#define hash_setkey(bucket,o) bucket->key = (void *)o;
#define hash_setvalue(bucket,o) bucket->value = (void *)o;
#define hash_init() (hash_init_size(64));

int hash_Dstring_hash (void *key,int hashtablesize);
int hash_Dstring_compare(const void* key1, const void* key2);
