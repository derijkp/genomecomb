/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#include "cg.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include "tools.h"
#include "tools_bcol.h"
#include "debug.h"

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

int bcol_NeedReversing(char format) {
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

void bcol_CopyNumber(
	const void *from,		  /* source */
	void *to,
	unsigned int length,	/* Number of bytes to copy */
	int reverse
) {
	switch (reverse) {
	case 0: 
		memcpy(to,from,length);
		break;
	case 1: {
		const unsigned char *fromPtr = from;
		unsigned char *toPtr = to;
		switch (length) {
		case 4:
			*(toPtr++) = fromPtr[3];
			*(toPtr++) = fromPtr[2];
			*(toPtr++) = fromPtr[1];
			*(toPtr++) = fromPtr[0];
			break;
		case 8:
			*(toPtr++) = fromPtr[7];
			*(toPtr++) = fromPtr[6];
			*(toPtr++) = fromPtr[5];
			*(toPtr++) = fromPtr[4];
			*(toPtr++) = fromPtr[3];
			*(toPtr++) = fromPtr[2];
			*(toPtr++) = fromPtr[1];
			*(toPtr++) = fromPtr[0];
			break;
		}
		break;
	}
	case 2: {
		const unsigned char *fromPtr = from;
		unsigned char *toPtr = to;
		*(toPtr++) = fromPtr[4];
		*(toPtr++) = fromPtr[5];
		*(toPtr++) = fromPtr[6];
		*(toPtr++) = fromPtr[7];
		*(toPtr++) = fromPtr[0];
		*(toPtr++) = fromPtr[1];
		*(toPtr++) = fromPtr[2];
		*(toPtr++) = fromPtr[3];
		break;
	}
	case 3: {
		const unsigned char *fromPtr = from;
		unsigned char *toPtr = to;
		*(toPtr++) = fromPtr[3];
		*(toPtr++) = fromPtr[2];
		*(toPtr++) = fromPtr[1];
		*(toPtr++) = fromPtr[0];
		*(toPtr++) = fromPtr[7];
		*(toPtr++) = fromPtr[6];
		*(toPtr++) = fromPtr[5];
		*(toPtr++) = fromPtr[4];
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

int bcol_printbin(FILE *f,int reverse,int isunsigned,char *type,char *string) {
	unsigned char buffer[8];
	long value;
	unsigned long uvalue;
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
		bcol_CopyNumber(&dvalue, buffer,sizeof(double), reverse);
		fwrite(buffer,sizeof(double),1,f);
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
		dvalue = strtod(string,NULL);
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
		bcol_CopyNumber(&fvalue, buffer,sizeof(float), reverse);
		fwrite(buffer,sizeof(float),1,f);
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
			uvalue = (long)strtoul(string,NULL,10);
			if (errno) {
				fprintf(stderr,"conversion error for type %s, value %s: %s.\n", type, string, strerror(errno));
				exit(EXIT_FAILURE);
			}
#if ULONG_MAX > 4294967295
			if (uvalue > 4294967295) {
				fprintf(stderr,"conversion error for type %s, value %s: Numerical result out of range.\n", type, string);
				exit(EXIT_FAILURE);
			}
#endif
		} else {
			value = strtol(string,NULL,10);
			if (errno) {
				fprintf(stderr,"conversion error for type %s, value %s: %s.\n", type, string, strerror(errno));
				exit(EXIT_FAILURE);
			}
#if LONG_MAX > 2147483647
			if (value > 2147483647 || value < -2147483648) {
				fprintf(stderr,"conversion error for type %s, value %s: Numerical result out of range.\n", type, string);
				exit(EXIT_FAILURE);
			}
#endif
		}
		if (isunsigned) {
			if (uvalue > UINT32_MAX) {errno = 1;}
			if (reverse) {
				putc_unlocked((unsigned char) uvalue,f);
				putc_unlocked((unsigned char) (uvalue >> 8),f);
				putc_unlocked((unsigned char) (uvalue >> 16),f);
				putc_unlocked((unsigned char) (uvalue >> 24),f);
			} else {
				putc_unlocked((unsigned char) (uvalue >> 24),f);
				putc_unlocked((unsigned char) (uvalue >> 16),f);
				putc_unlocked((unsigned char) (uvalue >> 8),f);
				putc_unlocked((unsigned char) uvalue,f);
			}
		} else {
			if (value < INT32_MIN || value > INT32_MAX) {errno = 1;}
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

char *clipzeros(char *sbuffer) {
	char *spointer, *end;
	spointer = sbuffer;
	while (*spointer != '\0' && *spointer != '.') {
		spointer++;
	}
	end = spointer;
	spointer++;
	while (*spointer != '\0') {
		if (*spointer != '0') {end=spointer+1;}
		spointer++;
	}
	*end = '\0';
	return sbuffer;
}

int bcol_printtext(FILE *f,int reverse,int isunsigned,char type,unsigned char *buffer,int precision) {
	long value;
	double dvalue;
	uint64_t wvalue;
	char *sbuffer[100];
	float fvalue;
	switch (type) {
	    case 'c':
		/*
		 * Characters need special handling. We want to produce a signed
		 * result, but on some platforms (such as AIX) chars are unsigned. To
		 * deal with this, check for a value that should be negative but
		 * isn't.
		 */
		value = buffer[0];
		if (!isunsigned) {
		    if (value & 0x80) {
			value |= -0x100;
		    }
		}
		goto returnNumericObject;
	    case 's':
	    case 'S':
	    case 't':
		/*
		 * 16-bit numeric values. We need the sign extension trick (see above)
		 * here as well.
		 */
		if (reverse) {
		    value = (long) (buffer[0] + (buffer[1] << 8));
		} else {
		    value = (long) (buffer[1] + (buffer[0] << 8));
		}
		if (!isunsigned) {
		    if (value & 0x8000) {
			value |= -0x10000;
		    }
		}
		goto returnNumericObject;
	    case 'i':
	    case 'I':
	    case 'n':
		/*
		 * 32-bit numeric values.
		 */
		if (reverse) {
		    value = (long) (buffer[0]
			    + (buffer[1] << 8)
			    + (buffer[2] << 16)
			    + (((long)buffer[3]) << 24));
		} else {
		    value = (long) (buffer[3]
			    + (buffer[2] << 8)
			    + (buffer[1] << 16)
			    + (((long)buffer[0]) << 24));
		}
	
		/*
		 * Check to see if the value was sign extended properly on systems
		 * where an int is more than 32-bits.
		 * We avoid caching unsigned integers as we cannot distinguish between
		 * 32bit signed and unsigned in the hash (short and char are ok).
		 */
	
		if (isunsigned) {
			fprintf(f,"%ld",value);
			return 1;
		}
		if ((value & (((unsigned int)1)<<31)) && (value > 0)) {
		    value -= (((unsigned int)1)<<31);
		    value -= (((unsigned int)1)<<31);
		}
	
	    returnNumericObject:
		fprintf(f,"%ld",value);
		return 1;
	    case 'w':
	    case 'W':
	    case 'm':
		/*
		 * Do not cache wide (64-bit) values; they are already too large to
		 * use as keys.
		 */
		if (reverse) {
		    wvalue = ((uint64_t) buffer[0])
			    | (((uint64_t) buffer[1]) << 8)
			    | (((uint64_t) buffer[2]) << 16)
			    | (((uint64_t) buffer[3]) << 24)
			    | (((uint64_t) buffer[4]) << 32)
			    | (((uint64_t) buffer[5]) << 40)
			    | (((uint64_t) buffer[6]) << 48)
			    | (((uint64_t) buffer[7]) << 56);
		} else {
		    wvalue = ((uint64_t) buffer[7])
			    | (((uint64_t) buffer[6]) << 8)
			    | (((uint64_t) buffer[5]) << 16)
			    | (((uint64_t) buffer[4]) << 24)
			    | (((uint64_t) buffer[3]) << 32)
			    | (((uint64_t) buffer[2]) << 40)
			    | (((uint64_t) buffer[1]) << 48)
			    | (((uint64_t) buffer[0]) << 56);
		}
		fprintf(f,"%" PRId64 "",wvalue);
/* need to handle unsigned properly ... (later)
		if (unsigned) {
			
		}
*/
		return 1;
	    case 'f':
	    case 'R':
	    case 'r':
		/*
		 * 32-bit IEEE single-precision floating point.
		 */
		bcol_CopyNumber(buffer, &fvalue, sizeof(float), reverse);
		if (precision == -1) {precision = 6;}
		sprintf((char *)sbuffer,"%.*f",precision,fvalue);
		fprintf(f,"%s",clipzeros((char *)sbuffer));
		return 1;
	    case 'd':
	    case 'Q':
	    case 'q':
		/*
		 * 64-bit IEEE double-precision floating point.
		 */
		bcol_CopyNumber(buffer, &dvalue, sizeof(double), reverse);
		if (precision == -1) {precision = 9;}
		sprintf((char *)sbuffer,"%.*f",precision,dvalue);
		fprintf(f,"%s",clipzeros((char *)sbuffer));
		return 1;
	    default:
		fprintf(stderr,"unknown type %s\n",&type);
		exit(EXIT_FAILURE);
	}
}

void bcol_close(BCol *bcol) {
	if (bcol->lz4 != NULL) {
		lz4_close(bcol->lz4);
	} else if (bcol->rz != NULL) {
		razf_close(bcol->rz);
	} else {
		fclose(bcol->f);
	}
	DStringDestroy(bcol->def);
	free(bcol->file);
	free(bcol->buffer);
	free(bcol);
}

int read_unlocked(FILE *f,unsigned char *buffer,int size) {
	int count = size,c;
	while (count--) {
		c = getc_unlocked(f);
		if (c == EOF) break;
		*buffer++ = c;
	}
	return(size-(count+1));
}

/* reads data from bcol, from position start to end (included): positions given in objects, so no need to * typesize */
int bcol_getbin(BCol *fbcol,uint64_t start,uint64_t end) {
	LZ4res *lz4 = fbcol->lz4;
	RAZF *rz = fbcol->rz;
	FILE *f = fbcol->f;
	int rsize,c;
	rsize = (end-start+1)*fbcol->typesize;
	if (fbcol->buffersize < rsize) {
		fbcol->buffer = realloc(fbcol->buffer,rsize);
		fbcol->buffersize = rsize;
	}
	start = (start-fbcol->start)*fbcol->typesize;
	if (start < 0) {
		c = 0;
	} else if (lz4 != NULL) {
		lz4_seek(lz4, start, SEEK_SET);
		c = lz4_read(lz4, fbcol->buffer, rsize);
	} else if (rz != NULL) {
		razf_seek(rz, start, SEEK_SET);
		c = razf_read(rz, fbcol->buffer, rsize);
	} else {
		fseeko(f, start, SEEK_SET);
		c = read_unlocked(f, fbcol->buffer, rsize);
	}
	if (c < rsize) {
		memset(fbcol->buffer+c,0,rsize-c);
	}
	return(c);
}

int bcol_getbinloc(BCol *bcol,DString *chromosome, uint64_t start,uint64_t end) {
	BCol_table *table = bcol->table;
	int tablesize = bcol->tablesize;
	int tablepos = bcol->tablepos;
	uint64_t binpos;
	DString *chromosome2 = table[tablepos].chr;
	int start2 = table[tablepos].begin;
	int end2 = table[tablepos].end;
	int comp;
	comp = DStringLocCompare(chromosome2,chromosome);
	while (tablepos < tablesize && ((comp < 0) || ((comp == 0) && ((end2 < start) || (end2 == start && start != end))))) {
		tablepos += 1;
		chromosome2 = table[tablepos].chr;
		start2 = table[tablepos].begin;
		end2 = table[tablepos].end;
		comp = DStringLocCompare(chromosome2,chromosome);
	}
	if (tablepos == tablesize) {
		tablepos = 0;
		while (tablepos < bcol->tablepos && ((comp < 0) || ((comp == 0) && ((end2 < start) || (end2 == start && start != end))))) {
			tablepos += 1;
			chromosome2 = table[tablepos].chr;
			start2 = table[tablepos].begin;
			end2 = table[tablepos].end;
			comp = DStringLocCompare(chromosome2,chromosome);
		}
		if (tablepos == bcol->tablepos) {return 0;}
	}
	if (comp> 0 || start < start2) {
		return 0;
	} else {
		binpos = table[tablepos].pos + start - start2;
		return bcol_getbin(bcol,binpos,binpos+end-start);
	}
}

int bcol_readbin(BCol *fbcol,int rsize,unsigned char *buffer) {
	LZ4res *lz4 = fbcol->lz4;
	RAZF *rz = fbcol->rz;
	FILE *f = fbcol->f;
	int c;
	if (lz4 != NULL) {
		c = lz4_read(lz4, buffer, rsize);
	} else if (rz != NULL) {
		c = razf_read(rz, buffer, rsize);
	} else {
		c = read_unlocked(f, buffer, rsize);
	}
	if (c < rsize) {
		memset(fbcol->buffer+c,0,rsize-c);
	}
	return(c);
}

int bcol_readdouble(BCol *fbcol,long double *result) {
	unsigned char buffer[9];
	long long value;
	float fvalue;
	double dvalue;
	int isunsigned = fbcol->isunsigned;
	int i, v;
	
	bcol_readbin(fbcol,fbcol->typesize,buffer);
	switch(fbcol->type) {
	case 'c':
		value = buffer[0];
		if (!isunsigned) {
			if (value & 0x80) {
				value |= -0x100;
			}
		}
		*result = (long double)value;
		break;
	case 's':
		if (fbcol->reverse) {
			value = (long long) (buffer[0] + (buffer[1] << 8));
		} else {
			value = (long long) (buffer[1] + (buffer[0] << 8));
		}
		if (!isunsigned) {
			if (value & 0x8000) {
				value |= -0x10000;
			}
		}
		*result = (long double)value;
		break;
	case 'i':
		if (fbcol->reverse) {
		value = (long long) (buffer[0]
			+ (buffer[1] << 8)
			+ (buffer[2] << 16)
			+ (((long long)buffer[3]) << 24));
		} else {
		value = (long long) (buffer[3]
			+ (buffer[2] << 8)
			+ (buffer[1] << 16)
			+ (((long long)buffer[0]) << 24));
		}
		if (!isunsigned) {
			if ((value & (((unsigned int)1)<<31)) && (value > 0)) {
			    value -= (((unsigned int)1)<<31);
			    value -= (((unsigned int)1)<<31);
			}
		}
		*result = (long double)value;
		break;
	case 'w':
		for(i = 0 ; i < 8 ; i++) {
			v = fgetc(stdin);
			if (v == EOF) {return 0;}
			buffer[i] = v;
		}
		if (fbcol->reverse) {
		value = ((long long) buffer[0])
			| (((long long) buffer[1]) << 8)
			| (((long long) buffer[2]) << 16)
			| (((long long) buffer[3]) << 24)
			| (((long long) buffer[4]) << 32)
			| (((long long) buffer[5]) << 40)
			| (((long long) buffer[6]) << 48)
			| (((long long) buffer[7]) << 56);
		} else {
		value = ((long long) buffer[7])
			| (((long long) buffer[6]) << 8)
			| (((long long) buffer[5]) << 16)
			| (((long long) buffer[4]) << 24)
			| (((long long) buffer[3]) << 32)
			| (((long long) buffer[2]) << 40)
			| (((long long) buffer[1]) << 48)
			| (((long long) buffer[0]) << 56);
		}
		*result = (long double)value;
		break;
	case 'f':
		for(i = 0 ; i < 4 ; i++) {
			v = fgetc(stdin);
			if (v == EOF) {return 0;}
			buffer[i] = v;
		}
		bcol_CopyNumber(buffer, &fvalue, sizeof(float), fbcol->reverse);
		*result = (long double)fvalue;
		break;
	case 'd':
		for(i = 0 ; i < 8 ; i++) {
			v = fgetc(stdin);
			if (v == EOF) {return 0;}
			buffer[i] = v;
		}
		bcol_CopyNumber(buffer, &dvalue, sizeof(double), fbcol->reverse);
		*result = (long double)dvalue;
		break;
	}
	return 1;
}

BCol *bcol_open(char *bcolfile) {
	FILE *f;
	BCol *result;
	DString *line = DStringNew();
	int typesize,len=strlen(bcolfile),i,cnt;
	result = (BCol *)malloc(sizeof(BCol));
	result->type = 'i';
	result->isunsigned = 1;
	result->table = NULL;
	result->multi = NULL;
	result->tablesize = 0;
	result->tablepos = 0;
	result->def = DStringNew();
	result->multi = DStringNew();
	result->precision = -1;
	f = fopen64_or_die(bcolfile, "r");
	while (1) {
		cnt = DStringGetLine(line,f);
		if (cnt == -1) break;
		if (line->string[0] != '#') break;
		if (line->size > 7 && strncmp(line->string,"# type ",7) == 0) {
			result->type = line->string[7];
			if (line->string[8] == 'u') {
				result->isunsigned = 1;
			} else {
				result->isunsigned = 0;
			}
		} else if (line->size > 10 && strncmp(line->string,"# default ",10) == 0) {
			DStringSetS(result->def,line->string+10,line->size-10);
		} else if (line->size > 8 && strncmp(line->string,"# multi ",8) == 0) {
			DStringSetS(result->multi,line->string+8,line->size-8);
		} else if (line->size > 12 && strncmp(line->string,"# precision ",12) == 0) {
			result->precision = atoi(line->string+12);
		}
	}
	if (naturalcompare(line->string,"chromosome\tbegin\tend",line->size,20) == 0) {
		result->version = 1;
		result->tablesize = 0;
		result->start = 0;
		if (result->table == NULL) {
			result->table = malloc(result->tablesize*sizeof(BCol_table));
		} else {
			result->table = realloc(result->table,result->tablesize*sizeof(BCol_table));
		}
		i = 0;
		while (1) {
			int strpos;
			cnt = DStringGetLine(line,f);
			if (cnt == -1) break;
			result->tablesize++;
			result->table = realloc(result->table,result->tablesize*sizeof(BCol_table));
			strpos = 0; while (strpos < line->size) {
				if (line->string[strpos] == '\t') break;
				strpos++;
			}
			result->table[i].chr = DStringNewFromCharS(line->string,strpos);
			if (sscanf(line->string+strpos+1,"%" PRId64 "\t%" PRId64 "",&(result->table[i].begin),&(result->table[i].end)) != 2) {
				fprintf(stderr,"error in bcol format");
				exit(1);
			}
			if (i == 0) {
				result->table[i].pos = 0;
			} else {
				result->table[i].pos = result->table[i-1].pos + result->table[i-1].end  - result->table[i-1].begin;
			}
			i++;
		}
	} else {
		result->version = 0;
		DStringGetLine(line,f);
		result->start = atoi(line->string);
		fclose(f);
	}
	DStringDestroy(line);
	result->file = (char *)malloc((len+9)*sizeof(char));
	strncpy(result->file,bcolfile,len);
	sprintf(result->file+len,".bin");
	result->f = fopen64(result->file, "r");
	if (result->f == NULL) {
		sprintf(result->file+len,".bin.lz4");
		result->lz4 = lz4_openfile(result->file,0);
		if (result->lz4 == NULL) {
			sprintf(result->file+len,".bin.rz");
			result->rz = razf_open(result->file, "r");
		} else {
			result->rz = NULL;
		}
	}
	switch (result->type) {
		case 'c':
			typesize = 1;
			break;
		case 's':
			typesize = 2;
			break;
		case 'i':
			typesize = 4;
			break;
		case 'w':
			typesize = 8;
			break;
		case 'f':
			typesize = 4;
			break;
		case 'd':
			typesize = 8;
			break;
		default:
			fprintf(stderr,"Unsupported bcol type: %c\n",result->type);
			exit(EXIT_FAILURE);
	}
	result->reverse = bcol_NeedReversing(result->type);
	result->typesize = typesize;
	result->buffer = (unsigned char *)malloc(typesize);
	result->buffersize = typesize;
	return(result);
}
