/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "debug.h"
#include "khash-dstring.h"

typedef struct dest {
	FILE *f;
	FILE *rf;
	uint64_t start;
	uint64_t lastpos;
} Dest;

KHASH_MAP_INIT_DSTR(DSTRING, Dest*);

/*#include "dstring-khash.h"*/

Dest *bcol_make_getout(khash_t(DSTRING) *hashtable,char *pre,DString *chromosome) {
	DString *buffer = NULL;
	Dest *o;
	khiter_t k;
	int ret;
	k = kh_put(DSTRING,hashtable, chromosome, &ret);
	if (ret == 0) {
		/* key was already present in the hashtable */
		o = kh_value(hashtable, k);
	} else {
		o = (Dest *)malloc(sizeof(Dest));
		o->start = 0;
		o->lastpos = -1;
		buffer = DStringNew();
		DStringAppend(buffer,pre);
		if (chromosome->size > 0) {
			DStringAppendS(buffer,chromosome->string,chromosome->size);
		}
		DStringAppend(buffer,".bcol");
		o->rf = fopen64(buffer->string,"w");
		if (o->rf == NULL) {
			fprintf(stderr,"Error opening file %s: %s.\n", buffer->string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		DStringAppend(buffer,".bin");
		o->f = fopen64(buffer->string,"w");
		if (o->f == NULL) {
			fprintf(stderr,"Error opening file %s: %s.\n", buffer->string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		kh_value(hashtable, k) = o;
		DStringSetS(buffer,chromosome->string,chromosome->size);
		kh_key(hashtable, k) = buffer;
		buffer = NULL;
	}
	return o;
}

/*
 *----------------------------------------------------------------------
 *
 * NeedReversing -- adapted from tclBinary.c
 * Copyright (c) 1997 by Sun Microsystems, Inc.
 * Copyright (c) 1998-1999 by Scriptics Corporation.
 * available under bsd license
 *
 *	This routine determines, if bytes of a number need to be re-ordered,
 *	and returns a numeric code indicating the re-ordering to be done.
 *	This depends on the endiannes of the machine and the desired format.
 *	It is in effect a table (whose contents depend on the endianness of
 *	the system) describing whether a value needs reversing or not. Anyone
 *	porting the code to a big-endian platform should take care to make
 *	sure that they define WORDS_BIGENDIAN though this is already done by
 *	configure for the Unix build; little-endian platforms (including
 *	Windows) don't need to do anything.
 *
 * Results:
 *	0	No re-ordering needed.
 *	1	Reverse the bytes:	01234567 <-> 76543210 (little to big)
 *	2	Apply this re-ordering: 01234567 <-> 45670123 (Nokia to little)
 *	3	Apply this re-ordering: 01234567 <-> 32107654 (Nokia to big)
 *
 * Side effects:
 *	None
 *
 *----------------------------------------------------------------------
 */

static int NeedReversing(int format) {
	switch (format) {
	/* native floats and doubles: never reverse */
	case 'd':
	case 'f':
	/* big endian ints: never reverse */
	case 'I':
	case 'S':
	case 'W':
#ifdef WORDS_BIGENDIAN
	/* native ints: reverse if we're little-endian */
	case 'n':
	case 't':
	case 'm':
	/* f: reverse if we're little-endian */
	case 'Q':
	case 'R':
#else /* !WORDS_BIGENDIAN */
	/* small endian floats: reverse if we're big-endian */
	case 'r':
#endif /* WORDS_BIGENDIAN */
	return 0;

#ifdef WORDS_BIGENDIAN
	/* small endian floats: reverse if we're big-endian */
	case 'q':
	case 'r':
#else /* !WORDS_BIGENDIAN */
	/* native ints: reverse if we're little-endian */
	case 'n':
	case 't':
	case 'm':
	/* f: reverse if we're little-endian */
	case 'R':
#endif /* WORDS_BIGENDIAN */
	/* small endian ints: always reverse */
	case 'i':
	case 's':
	case 'w':
		return 1;

#ifndef WORDS_BIGENDIAN
	/*
	 * The Q and q formats need special handling to account for the unusual
	 * byte ordering of 8-byte floats on Nokia 770 systems, which claim to be
	 * little-endian, but also reverse word order.
	 */

	case 'Q':
		return 1;
	case 'q':
		return 0;
#endif
	/* char: never reverse */
	case 'c':
		return 0;
	}

	fprintf(stderr,"unexpected fallthrough, do not know type \"%c\"\n",format);
	exit(EXIT_FAILURE);
}

/*
 *----------------------------------------------------------------------
 *
 * CopyNumber -- adapted from tclBinary.c
 * Copyright (c) 1997 by Sun Microsystems, Inc.
 * Copyright (c) 1998-1999 by Scriptics Corporation.
 * available under bsd license
 *
 *	This routine is called by FormatNumber and ScanNumber to copy a
 *	floating-point number. If required, bytes are reversed while copying.
 *	The behaviour is only fully defined when used with IEEE float and
 *	double values (guaranteed to be 4 and 8 bytes long, respectively.)
 *
 * Results:
 *	None
 *
 * Side effects:
 *	Copies length bytes
 *
 *----------------------------------------------------------------------
 */

static void CopyNumber(
	FILE *f,
	int reverse,
	const void *from,		  /* source */
	unsigned int length,	/* Number of bytes to copy */
	int type			 /* What type of thing are we copying? */
) {
	switch (reverse) {
	case 0: 
		fwrite(from,length,1,f);
		break;
	case 1: {
		const unsigned char *fromPtr = from;
		switch (length) {
		case 4:
			putc_unlocked((unsigned char) fromPtr[3],f);
			putc_unlocked((unsigned char) fromPtr[2],f);
			putc_unlocked((unsigned char) fromPtr[1],f);
			putc_unlocked((unsigned char) fromPtr[0],f);
			break;
		case 8:
			putc_unlocked((unsigned char) fromPtr[7],f);
			putc_unlocked((unsigned char) fromPtr[6],f);
			putc_unlocked((unsigned char) fromPtr[5],f);
			putc_unlocked((unsigned char) fromPtr[4],f);
			putc_unlocked((unsigned char) fromPtr[3],f);
			putc_unlocked((unsigned char) fromPtr[2],f);
			putc_unlocked((unsigned char) fromPtr[1],f);
			putc_unlocked((unsigned char) fromPtr[0],f);
			break;
		}
		break;
	}
	case 2: {
		const unsigned char *fromPtr = from;
		putc_unlocked((unsigned char) fromPtr[4],f);
		putc_unlocked((unsigned char) fromPtr[5],f);
		putc_unlocked((unsigned char) fromPtr[6],f);
		putc_unlocked((unsigned char) fromPtr[7],f);
		putc_unlocked((unsigned char) fromPtr[0],f);
		putc_unlocked((unsigned char) fromPtr[1],f);
		putc_unlocked((unsigned char) fromPtr[2],f);
		putc_unlocked((unsigned char) fromPtr[3],f);
		break;
	}
	case 3: {
		const unsigned char *fromPtr = from;
		putc_unlocked((unsigned char) fromPtr[3],f);
		putc_unlocked((unsigned char) fromPtr[2],f);
		putc_unlocked((unsigned char) fromPtr[1],f);
		putc_unlocked((unsigned char) fromPtr[0],f);
		putc_unlocked((unsigned char) fromPtr[7],f);
		putc_unlocked((unsigned char) fromPtr[6],f);
		putc_unlocked((unsigned char) fromPtr[5],f);
		putc_unlocked((unsigned char) fromPtr[4],f);
		break;
	}
	}
}

/*
 *----------------------------------------------------------------------
 *
 * FormatNumber -- adapted from tclBinary.c
 * Copyright (c) 1997 by Sun Microsystems, Inc.
 * Copyright (c) 1998-1999 by Scriptics Corporation.
 * available under bsd license
 *
 *	This routine is called by Tcl_BinaryObjCmd to format a number into a
 *	location pointed at by cursor.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	Moves the cursor to the next location to be written into.
 *
 *----------------------------------------------------------------------
 */

int FormatNumber(FILE *f,int reverse,int isunsigned,char *type,char *string) {
	long value;
	double dvalue;
	uint64_t wvalue;
	float fvalue;
	switch (type[0]) {
	case 'd':
	case 'q':
	case 'Q':
		errno = 0;
		dvalue = strtod(string,NULL);
		if (errno) {
			fprintf(stderr,"conversion error for type %s, value %s: %s.\n", type, string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		CopyNumber(f,reverse,&dvalue, sizeof(double), type[0]);
		return 1;
	case 'f':
	case 'r':
	case 'R':
		/*
		 * Single-precision floating point values. Tcl_GetDoubleFromObj
		 * returns TCL_ERROR for NaN, but we can check by comparing the
		 * object's type pointer.
		 */
		errno = 0;
		dvalue = strtof(string,NULL);
		if (errno) {
			fprintf(stderr,"conversion error for type %s, value %s: %s.\n", type, string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		/*
		 * Because some compilers will generate floating point exceptions on
		 * an overflow cast (e.g. Borland), we restrict the values to the
		 * valid range for float.
		 */
		fvalue = (float) dvalue;
		CopyNumber(f,reverse,&fvalue, sizeof(float), type[0]);
		return 1;
	/*
	 * 64-bit integer values.
	 */
	case 'w':
	case 'W':
	case 'm':
		errno = 0;
		if (isunsigned) {
			if (string[0] == '-') {
				fprintf(stderr,"conversion error for type %s, value %s: Numerical result out of range.\n", type, string);
				exit(EXIT_FAILURE);
			}
			wvalue = (uint64_t)strtoull(string,NULL,10);
		} else {
			wvalue = strtoll(string,NULL,10);
		}
		if (errno) {
			fprintf(stderr,"conversion error for type %s, value %s: %s.\n", type, string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		if (reverse) {
			putc_unlocked((unsigned char) wvalue,f);
			putc_unlocked((unsigned char) (wvalue >> 8),f);
			putc_unlocked((unsigned char) (wvalue >> 16),f);
			putc_unlocked((unsigned char) (wvalue >> 24),f);
			putc_unlocked((unsigned char) (wvalue >> 32),f);
			putc_unlocked((unsigned char) (wvalue >> 40),f);
			putc_unlocked((unsigned char) (wvalue >> 48),f);
			putc_unlocked((unsigned char) (wvalue >> 56),f);
		} else {
			putc_unlocked((unsigned char) (wvalue >> 56),f);
			putc_unlocked((unsigned char) (wvalue >> 48),f);
			putc_unlocked((unsigned char) (wvalue >> 40),f);
			putc_unlocked((unsigned char) (wvalue >> 32),f);
			putc_unlocked((unsigned char) (wvalue >> 24),f);
			putc_unlocked((unsigned char) (wvalue >> 16),f);
			putc_unlocked((unsigned char) (wvalue >> 8),f);
			putc_unlocked((unsigned char) wvalue,f);
		}
		return 1;
	/*
	 * 32-bit integer values.
	 */
	case 'i':
	case 'I':
	case 'n':
		errno = 0;
		if (isunsigned) {
			if (string[0] == '-') {
				fprintf(stderr,"conversion error for type %s, value %s: Numerical result out of range.\n", type, string);
				exit(EXIT_FAILURE);
			}
			value = (long)strtoul(string,NULL,10);
		} else {
			value = strtol(string,NULL,10);
		}
		if (errno) {
			fprintf(stderr,"conversion error for type %s, value %s: %s.\n", type, string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		if (reverse) {
			putc_unlocked((unsigned char) value,f);
			putc_unlocked((unsigned char) (value >> 8),f);
			putc_unlocked((unsigned char) (value >> 16),f);
			putc_unlocked((unsigned char) (value >> 24),f);
		} else {
			putc_unlocked((unsigned char) (value >> 24),f);
			putc_unlocked((unsigned char) (value >> 16),f);
			putc_unlocked((unsigned char) (value >> 8),f);
			putc_unlocked((unsigned char) value,f);
		}
		return 1;
	/*
	 * 16-bit integer values.
	 */
	case 's':
	case 'S':
	case 't':
		errno = 0;
		value = strtol(string,NULL,10);
		if (errno) {
			fprintf(stderr,"conversion error for type %s, value %s: %s.\n", type, string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		if (value > 65535 || (!isunsigned && value > 32767)) {
			fprintf(stderr,"conversion error for type %s: value %ld too large.\n",type, value);
			exit(EXIT_FAILURE);
		} else if ((isunsigned && value < 0) || value < -32768) {
			fprintf(stderr,"conversion error for type %s: value %ld too small.\n",type, value);
			exit(EXIT_FAILURE);
		}
		if (reverse) {
			putc_unlocked((unsigned char) value,f);
			putc_unlocked((unsigned char) (value >> 8),f);
		} else {
			putc_unlocked((unsigned char) (value >> 8),f);
			putc_unlocked((unsigned char) value,f);
		}
		return 1;
		/*
		 * 8-bit integer values.
		 */
	case 'c':
		errno = 0;
		value = strtol(string,NULL,10);
		if (errno) {
			fprintf(stderr,"conversion error for type %s, value %s: %s.\n", type, string, strerror(errno));
			exit(EXIT_FAILURE);
		}
		if (value > 255 || (!isunsigned && value > 127)) {
			fprintf(stderr,"conversion error for type %s: value %ld too large.\n",type, value);
			exit(EXIT_FAILURE);
		} else if ((isunsigned && value < 0) || value < -127) {
			fprintf(stderr,"conversion error for type %s: value %ld too small.\n",type, value);
			exit(EXIT_FAILURE);
		}
		putc_unlocked((unsigned char) value,f);
		return 1;
	default:
		fprintf(stderr,"unknown type %s\n",type);
		exit(EXIT_FAILURE);
	}
}


int main(int argc, char *argv[]) {
	khash_t(DSTRING) *hashtable;
	Dest *o,*po = NULL;
	DString *result = NULL;
	DString *line = NULL,*chromosome = NULL;
	khiter_t k;
	char *pre,*type,*defaultvalue;
	uint64_t offset, poffset = -1, size;
	int reverse = 0, isunsigned = 0;
	int col = 0,max = 0,offsetcol = -1,chrcol = -1;
	if ((argc < 2)||(argc > 8)) {
		fprintf(stderr,"Format is: bcol_make output_pre type ?col? ?chromosomecol? ?offsetcol? ?default?\n");
		exit(EXIT_FAILURE);
	}
	pre = argv[1];
	if (argc >= 3) {
		type = argv[2];
	}
	if (argc >= 4) {
		col = atoi(argv[3]);
		max = col;
	}
	if (argc >= 5) {
		chrcol = atoi(argv[4]);
		if (chrcol > max) {max = chrcol;}
	}
	if (argc >= 6) {
		offsetcol = atoi(argv[5]);
		if (offsetcol > max) {max = offsetcol;}
	}
	if (argc >= 7) {
		defaultvalue = argv[6];
	}
NODPRINT("bcol_make %s %s %d %d %d\n",pre,type,col,chrcol,offsetcol)
	line = DStringNew();
	reverse  = NeedReversing((int)type[0]);
	if (type[1] == 'u') {isunsigned = 1;}
	result = DStringArrayNew(max+1);
	hashtable = kh_init(DSTRING);
	if (chrcol == -1) {
		o = bcol_make_getout(hashtable,pre,DStringEmtpy());
		poffset = -1;
	}
	while (!DStringGetTab(line,stdin,max,result,0)) {
		if (chrcol != -1) {
			chromosome = result+chrcol;
			o = bcol_make_getout(hashtable,pre,chromosome);
			if (o != po) {
				po = o;
				poffset = -1;
			}
		}
		if (offsetcol != -1) {
			offset = atoll(result[offsetcol].string);
			if (poffset == -1) {
				o->start = offset;
			} else if (poffset != offset) {
				size = offset - poffset;
				if (size < 0) {
					fprintf(stderr, "error: cannot make position based bcol on unsorted data ($offset < $poffset sort on position first)\n");
					exit(EXIT_FAILURE);
				}
				while (poffset < offset) {
					FormatNumber(o->f,reverse,isunsigned,type,defaultvalue);
					poffset++;
				}
				o->lastpos = o->lastpos + size;
			}
			poffset = offset+1;
		}
		NODPRINT("s=%s\n",result[col].string)
		FormatNumber(o->f,reverse,isunsigned,type,result[col].string);
		o->lastpos ++;
	}
	for (k = kh_begin(hashtable); k != kh_end(hashtable); ++k) {
		if (kh_exist(hashtable, k)) {
			DStringDestroy(kh_key(hashtable, k));
			o = kh_value(hashtable, k);
			fclose(o->f);
			fprintf(o->rf,"# binary column\n");
			fprintf(o->rf,"# type %s\n",type);
			fprintf(o->rf,"# default %s\n","0");
			fprintf(o->rf,"begin\ttype\toffset\n");
			fprintf(o->rf,"%llu\t%s\t%d\n",o->start,type,0);
			fprintf(o->rf,"%llu\tend\t%d\n",o->start + o->lastpos,0);
			fclose(o->rf);
			free(o);
		}
	}
	if (line) {DStringDestroy(line);}
	if (result) {DStringArrayDestroy(result);}
	kh_destroy(DSTRING,hashtable);
	exit(EXIT_SUCCESS);
}

