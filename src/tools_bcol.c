/*
 * Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
 * See the file "license.txt" for information on usage and redistribution of
 * this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
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

int bcol_NeedReversing(int format) {
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
		} else {
			value = strtol(string,NULL,10);
		}
		if (errno) {
			fprintf(stderr,"conversion error for type %s, value %s: %s.\n", type, string, strerror(errno));
			exit(EXIT_FAILURE);
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

int bcol_printtext(FILE *f,int reverse,int isunsigned,char *type,char *buffer) {
	long value;
	double dvalue;
	uint64_t wvalue;
	float fvalue;
	switch (type[0]) {
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
		fprintf(f,"%lld",wvalue);
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
		fprintf(f,"%f",fvalue);
		return 1;
	    case 'd':
	    case 'Q':
	    case 'q':
		/*
		 * 64-bit IEEE double-precision floating point.
		 */
		bcol_CopyNumber(buffer, &dvalue, sizeof(double), reverse);
		fprintf(f,"%f",fvalue);
		return 1;
	    default:
		fprintf(stderr,"unknown type %s\n",type);
		exit(EXIT_FAILURE);
	}
}

BCol *bcol_open(char *bcolfile,char *bcoltype) {
	BCol *result;
	int typesize,len=strlen(bcolfile);
	result = (BCol *)malloc(sizeof(BCol));
	if (len >= 3 && (strncmp(bcolfile+len-3,".rz",3)) == 0) {
		result->rz = razf_open(bcolfile, "r");
		result->f = NULL;
	} else {
		result->f = fopen(bcolfile, "r");
		result->rz = NULL;
	}
	result->type = bcoltype;
	switch (bcoltype[0]) {
		case 'c':
			typesize = 1;
			break;
		case 's':
			typesize = 2;
			break;
		case 'i':
			typesize = 4;
			break;
		default:
			fprintf(stderr,"Unsupported bcol type: %s",bcoltype);
			exit(EXIT_FAILURE);
	}
	result->reverse = bcol_NeedReversing((int)bcoltype[0]);
	if (bcoltype[1] == 'u') {
		result->isunsigned = 1;
	} else {
		result->isunsigned = 0;
	}
	result->file = bcolfile;
	result->type = bcoltype;
	result->typesize = typesize;
	result->buffer = (char *)malloc(typesize);
	result->buffersize = typesize;
	return(result);
}

void bcol_close(BCol *bcol) {
	if (bcol->rz != NULL) {
		razf_close(bcol->rz);
	} else {
		fclose(bcol->f);
	}
	free(bcol->buffer);
	free(bcol);
}

int bcol_get(BCol *fbcol,int start,int end) {
	RAZF *rz = fbcol->rz;
	int rsize,c;
	razf_seek(rz, start, SEEK_SET);
	rsize = (end-start+1)*fbcol->typesize;
	if (fbcol->buffersize < rsize) {
		fbcol->buffer = realloc(fbcol->buffer,rsize);
		fbcol->buffersize = rsize;
	}
	c = razf_read(rz, fbcol->buffer, rsize);
	return(c);
}
