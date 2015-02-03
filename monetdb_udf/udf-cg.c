/*
 * The contents of this file are subject to the MonetDB Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.monetdb.org/Legal/MonetDBLicense
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * The Original Code is the MonetDB Database System.
 *
 * The Initial Developer of the Original Code is CWI.
 * Portions created by CWI are Copyright (C) 1997-July 2008 CWI.
 * Copyright August 2008-2013 MonetDB B.V.
 * All Rights Reserved.
 */

/* monetdb_config.h must be the first include in each .c file */
#include "monetdb_config.h"
#include "udf-cg.h"
#include "float.h"

/* actual implementation */
/* all non-exported functions must be declared static */
static str
UDFlmin_(str *ret, str src)
{
	size_t len = 0;
	double min=INFINITY,value;
	char *cur,*next,*end,*result;
	int resultlen;
	str dst = NULL;
	/* assert calling sanity */
	assert(ret != NULL);
	/* handle NULL pointer and NULL value */
	if (src == NULL || strcmp(src, str_nil) == 0) {
		*ret = GDKstrdup(str_nil);
		if (*ret == NULL)
			throw(MAL, "udf.lmin",
			      "failed to create copy of str_nil");
		return MAL_SUCCEED;
	}
	cur = src;
	next = cur;
	result = cur;
	resultlen = 0;
	while (1) {
		if (*next == ',' || *next == ';' || *next == ' ' || *next == '\0') {
			value = strtod(cur,&end);
			if (end == next && end != cur) {
				if (resultlen == 0 || value < min) {
					min = value;
					result=cur; resultlen=next-cur;
				}
				cur=next+1; next = cur;
			}
		}
		if (*next == '\0') break;
		next++;
	}
	if (resultlen == 0) {
		*ret = GDKstrdup(str_nil);
		if (*ret == NULL)
			throw(MAL, "udf.lmin",
			      "failed to create copy of str_nil");
		return MAL_SUCCEED;
	}
	*ret = dst = GDKmalloc(resultlen+1);
	if (dst == NULL) {
		throw(MAL, "udf.lmin", "failed to allocate string ");
	}
	strncpy(dst, result, resultlen);
	dst[resultlen]='\0';
	return MAL_SUCCEED;
}

/* MAL wrapper */
str
UDFlmin(str *ret, str *arg)
{
	/* assert calling sanity */
	assert(ret != NULL && arg != NULL);
	return UDFlmin_ ( ret, *arg );
}

/* lmin a BAT of strings */
/*
 * Generic "type-oblivious" version,
 * using generic "type-oblivious" BAT access interface.
 */

#ifdef NEVER
/* actual implementation */
static str
UDFBATlmin_(BAT **ret, BAT *left)
{
	BATiter li;
	BAT *bn = NULL;
	BUN p = 0, q = 0;
	/* assert calling sanity */
	assert(ret != NULL);
	/* handle NULL pointer */
	if (left == NULL)
		throw(MAL, "batudf.lmin", RUNTIME_OBJECT_MISSING);
	/* check tail type */
	if (left->ttype != TYPE_str) {
		throw(MAL, "batudf.lmin",
		      "tail-type of input BAT must be TYPE_str");
	}
	/* allocate result BAT */
	bn = BATnew(left->htype, TYPE_dbl, BATcount(left));
	if (bn == NULL) {
		throw(MAL, "batudf.lmin", MAL_MALLOC_FAIL);
	}
	BATseqbase(bn, left->hseqbase);
	/* create BAT iterator */
	li = bat_iterator(left);
	/* advice on sequential scan */
	BATaccessBegin(left, USE_HEAD | USE_TAIL, MMAP_SEQUENTIAL);
	/* the core of the algorithm, expensive due to malloc/frees */
	BATloop(left, p, q) {
		str tr = NULL, err = NULL;
		/* get original head & tail value */
		ptr h = BUNhead(li, p);
		str t = (str) BUNtail(li, p);
		/* revert tail value */
		err = UDFlmin_(&tr, t);
		if (err != MAL_SUCCEED) {
			/* error -> bail out */
			BATaccessEnd(left, USE_HEAD | USE_TAIL,
				     MMAP_SEQUENTIAL);
			BBPreleaseref(bn->batCacheid);
			return err;
		}
		/* assert logical sanity */
		assert(tr != NULL);
		/* insert original head and lmin tail in result BAT */
		/* BUNins() takes care of all necessary administration */
		BUNins(bn, h, tr, FALSE);
		/* free memory allocated in UDFlmin_() */
		GDKfree(tr);
	}
	BATaccessEnd(left, USE_HEAD | USE_TAIL, MMAP_SEQUENTIAL);
	*ret = bn;
	return MAL_SUCCEED;
}

/* MAL wrapper */
str
UDFBATlmin(bat *ret, bat *bid)
{
	BAT *res = NULL, *left = NULL;
	str msg = NULL;
	/* assert calling sanity */
	assert(ret != NULL && bid != NULL);
	/* bat-id -> BAT-descriptor */
	if ((left = BATdescriptor(*bid)) == NULL)
		throw(MAL, "batudf.lmin", RUNTIME_OBJECT_MISSING);
	/* do the work */
	msg = UDFBATlmin_ ( &res, left );
	/* release input BAT-descriptor */
	BBPreleaseref(left->batCacheid);
	if (msg == MAL_SUCCEED) {
		/* register result BAT in buffer pool */
		BBPkeepref((*ret = res->batCacheid));
	}
	return msg;
}

#endif

/* actual implementation */
/* all non-exported functions must be declared static */
static str
UDFlmax_(str *ret, str src)
{
	size_t len = 0;
	double min=INFINITY,value;
	char *cur,*next,*end,*result;
	int resultlen;
	str dst = NULL;
	/* assert calling sanity */
	assert(ret != NULL);
	/* handle NULL pointer and NULL value */
	if (src == NULL || strcmp(src, str_nil) == 0) {
		*ret = GDKstrdup(str_nil);
		if (*ret == NULL)
			throw(MAL, "udf.lmax",
			      "failed to create copy of str_nil");
		return MAL_SUCCEED;
	}
	cur = src;
	next = cur;
	result = cur;
	resultlen = 0;
	while (1) {
		if (*next == ',' || *next == ';' || *next == ' ' || *next == '\0') {
			value = strtod(cur,&end);
			if (end == next && end != cur) {
				if (resultlen == 0 || value > min) {
					min = value;
					result=cur; resultlen=next-cur;
				}
				cur=next+1; next = cur;
			}
		}
		if (*next == '\0') break;
		next++;
	}
	if (resultlen == 0) {
		*ret = GDKstrdup(str_nil);
		if (*ret == NULL)
			throw(MAL, "udf.lmin",
			      "failed to create copy of str_nil");
		return MAL_SUCCEED;
	}
	*ret = dst = GDKmalloc(resultlen+1);
	if (dst == NULL) {
		throw(MAL, "udf.lmax", "failed to allocate string ");
	}
	strncpy(dst, result, resultlen);
	dst[resultlen]='\0';
	return MAL_SUCCEED;
}

/* MAL wrapper */
str
UDFlmax(str *ret, str *arg)
{
	/* assert calling sanity */
	assert(ret != NULL && arg != NULL);
	return UDFlmax_ ( ret, *arg );
}

