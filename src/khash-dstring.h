#include "khash.h"

typedef DString *kh_dstr_t;

static kh_inline khint_t __kh_dstr_hash_string(DString *ds)
{
	char *s = ds->string;
	int count = ds->size;
	khint_t h = (khint_t) s;
	if (count <= 0) {
		h = 0;
	} else {
		s++;
		while (--count) {
			h = (h << 5) - h + (khint_t) * s;
			s++;
		}
	}
	return h;
}

#define kh_dstr_hash_func(key) __kh_dstr_hash_string(key)
#define kh_dstr_hash_equal(a, b) (DStringCompare(a, b) == 0)

#define KHASH_SET_INIT_DSTR(name, khval_t)								\
	KHASH_INIT(name, kh_dstr_t, khval_t, 0, kh_dstr_hash_func, kh_dstr_hash_equal)

#define KHASH_MAP_INIT_DSTR(name, khval_t)								\
	KHASH_INIT(name, kh_dstr_t, khval_t, 1, kh_dstr_hash_func, kh_dstr_hash_equal)

