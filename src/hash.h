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
typedef unsigned int hash_hash_func(void* key, unsigned int hashtablesize);
typedef void hash_free_func(void* key);

Hash_table *hash_init_size(int tablesize);
void hash_destroy(Hash_table *table,hash_free_func *freekeyfunc, hash_free_func *freevaluefunc);
Hash_bucket *hash_get(Hash_table *table, void *key, hash_hash_func *hashfunc, hash_key_compare_func *compare, int *new, int create);
Hash_bucket *hash_first(Hash_table *table,Hash_iter *bucket);
Hash_bucket *hash_next(Hash_iter *iter);

#define hash_getkey(bucket) (bucket->key);
#define hash_getvalue(bucket) (bucket->value);
#define hash_setkey(bucket,o) bucket->key = (void *)o;
#define hash_setvalue(bucket,o) bucket->value = (void *)o;
#define hash_init() (hash_init_size(64));

unsigned int hash_DString_hash(void *key,unsigned int hashtablesize);
int hash_DString_compare(const void* key1, const void* key2);

int char_hash_set(Hash_table *hashtable,char *key,void *value);
void *char_hash_get(Hash_table *hashtable,char *key);
int dstring_hash_set(Hash_table *hashtable,DString *key,void *value);
void *dstring_hash_get(Hash_table *hashtable,DString *key);
