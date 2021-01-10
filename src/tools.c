/* #define DEBUG 1 */

#ifndef TOOLS_LOADED
#define TOOLS_LOADED 1

#define _FILE_OFFSET_BITS 64

#define _GNU_SOURCE
#include "cg.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include "tools.h"
#include "debug.h"

void DStringInit(DString *dstring) {
	dstring->memsize = DSTRING_STATICLEN-1;
	dstring->size = 0;
	dstring->string = dstring->staticspace;
	dstring->string[0]='\0';
}

DString *DStringNew() {
	DString *dstring = (DString *)malloc(sizeof(DString));
	DStringInit(dstring);
	return dstring;
}

static DString empty_dstring;

DString *DStringEmtpy() {
	DStringInit(&empty_dstring);
	return &empty_dstring;
}

void DStringClear(DString *dstring) {
	if (dstring->string != dstring->staticspace && dstring->memsize != -1) {
		free(dstring->string);
		dstring->string = dstring->staticspace;
	}
	dstring->memsize = DSTRING_STATICLEN;
	dstring->size = 0;
	dstring->staticspace[0]='\0';
}

void DStringDestroy(DString *dstring) {
	if (dstring == NULL) return;
	if (dstring == &empty_dstring) return;
	DStringClear(dstring);
	free(dstring);
}

void DStringSetSize(DString *dstring, int size) {
	int ssize = dstring->size;
	if (ssize > size) {ssize = size;}
	size++;
	if (dstring->memsize < size) {
		if (dstring->string == dstring->staticspace) {
			dstring->string = (char *)malloc(size);
			strncpy(dstring->string,dstring->staticspace,ssize);
			dstring->string[ssize] = '\0';
		} else if (dstring->memsize == -1) {
			char *temp;
			temp = (char *)malloc(size);
			strncpy(temp,dstring->string,ssize);
			dstring->string = temp;
			dstring->string[ssize] = '\0';
		} else {
			dstring->string = (char *)realloc(dstring->string,size);
		}
		dstring->memsize = size;
	}
}

void DStringSet(DString *dstring, char *string) {
	int size = strlen(string);
	DStringSetSize(dstring,size);
	strncpy(dstring->string,string,size+1);
	dstring->size = size;
}

void DStringSetS(DString *dstring, char *string, int size) {
	DStringSetSize(dstring,size);
	strncpy(dstring->string,string,size+1);
	dstring->string[size] = '\0';
	dstring->size = size;
}

void DStringCopy(DString *dest, DString *src) {
	if (src == NULL) {
		dest->string[0] = '\0';
		dest->size = 0;
	} else {
		DStringSetSize(dest,src->size);
		strncpy(dest->string,src->string,src->size+1);
		dest->size = src->size;
	}
}

void DStringputs(DString *string,FILE *f) {
	char *cur;
	int size;
	if (string == NULL) return;
	cur = string->string;
	size = string->size;
	if (size > 0) {
		while(size--) {putc_unlocked(*cur++,f);}
	}
}

void DStringArrayPuts(DStringArray *array,char *join,FILE *f) {
	char *cur;
	int size=strlen(join),pos=0,count;
	if (array->size == 0) return;
	DStringputs(array->data+0,f);
	for(pos = 1; pos < array->size ; pos++) {
		cur = join; count = size;
		while(count--) {putc_unlocked(*cur++,f);}
		DStringputs(array->data+pos,f);
	}
}

void charputs(char *cur, int size,FILE *f) {
	if (size > 0) {
		while(size--) {putc_unlocked(*cur++,stdout);}
	}
}

void DStringPrintf(DString *dstring, char *format, ...) {
	va_list args;
	int size;
	va_start(args, format);
	size = vsnprintf(NULL,0, format, args);
	va_end(args);
	DStringSetSize(dstring,dstring->size+size);
	va_start(args, format);
	vsprintf(dstring->string+dstring->size, format, args);
	dstring->size += size;
	va_end(args);
}

DString *DStringNewFromChar(char *string) {
	DString *dstring = DStringNew();
	DStringSet(dstring,string);
	return dstring;
}

DString *DStringNewFromCharS(char *string, int size) {
	DString *dstring = DStringNew();
	DStringSetS(dstring,string,size);
	return dstring;
}

DString *DStringDup(DString *dstring) {
	DString *result = DStringNew();
	DStringSetS(result,dstring->string,dstring->size);
	return result;
}

DString *DStringNewFromInt(int i) {
	DString *dstring = DStringNew();
	int size=snprintf(NULL,0,"%d",i);
	DStringSetSize(dstring,size);
	sprintf(dstring->string,"%d",i);
	return dstring;
}

/* 
    "natural" sort order: deals with mixed alphabetic and numbers
 */

#define NM_NOTINNUM -1
#define NM_NOTNUM 0
#define NM_DIGIT 1
#define NM_MINUS 1
#define NM_DECIMAL 2
#define NM_PLUS 3
#define NM_E 4

#define LOC_UNKOWN -1
#define LOC_NONUM 0
#define LOC_SIMPLE 1
#define LOC_BEFORE 2
#define LOC_SIGN 3
#define LOC_SIGNNUM 4
#define LOC_DECIMAL 5
#define LOC_DECIMALNUM 6
#define LOC_E 7
#define LOC_ESIGN 8
#define LOC_ENUM 9

#define UCHAR(c) ((unsigned char) (c))
/*#define NATDIGIT(c) (isdigit(UCHAR(*(c))))*/
#define NATDIGIT(c) (*(c) > 47 && *(c)<58)
#define ISNUMBER(c) ((*(c) == '-')?NM_MINUS \
	:(*(c) == '.')?NM_DECIMAL \
	:(*(c) == '+')?NM_PLUS \
	:(*(c) == '.')?NM_DECIMAL \
	:(*(c) == 'e' || *(c) == 'E')?NM_E \
	:(*(c) > 47 && *(c)<58)?NM_DIGIT \
	:0)
#define BLANK(char) (*(char) == ' ' || *(char) == '\t' || *(char) == '\n')

/* 
	returns the start of the complete number, or NULL if there is someting (double e, etc.) preventing a complex number
*/
char *naturalcompare_numbercontext_before(char const *a, char *left, int *invert, int *context) {
	char *cur, *start;
	int digits;
	if (left == a) {return left;}
	cur = left - 1;
	digits = 0;
	while (cur >= a && NATDIGIT(cur)) {cur--; digits++;}
	/* if cur < a, number starts at the beginning with a digit -> no further context */
	*context = LOC_UNKOWN;
	if (cur < a) {
		if (digits) {
			*context = LOC_SIGNNUM;
		} else {
			*context = LOC_BEFORE;
		}
		return(cur+1);
	}
	start = cur+1;
	if (*cur == '-' || *cur == '+') {
		if (*cur == '-') {*invert = 1;}
		if (cur == a || BLANK(cur-1)) {
			if (digits) {
				*context = LOC_SIGNNUM;
			} else {
				*context = LOC_SIGN;
			}
			return cur;
		} else if (*(cur-1) == 'e' || *(cur-1) == 'E') {
			cur--;
			if (cur == a || BLANK(cur)) {
				*context = LOC_SIMPLE;
				return start;
			} else if (digits) {
				*context = LOC_ENUM;
			} else {
				*context = LOC_ESIGN;
			}
			cur--; digits = 0;
			while (cur >= a && NATDIGIT(cur)) {cur--; digits++;}
			if (!digits) {*context = LOC_SIMPLE; return start;}
			if (cur < a) {return (char *)a;}
		} else {
			*context = LOC_SIMPLE;
			return start;
		}
	} else if (*cur == 'e' || *cur == 'E') {
		if (cur == a || BLANK(cur)) {
			*context = LOC_SIMPLE;
			return start;
		} else if (digits) {
			*context = LOC_ENUM;
		} else {
			*context = LOC_E;
		}
		cur--; digits = 0;
		while (cur >= a && NATDIGIT(cur)) {cur--; digits++;}
		if (!digits) {*context = LOC_SIMPLE; return start;}
		if (cur < a) {return (char *)a;}
	}
	if (*cur == '.') {
		if (*context == LOC_UNKOWN) {
			if (digits) {
				*context = LOC_DECIMALNUM;
			} else {
				*context = LOC_DECIMAL;
			}
			if (cur == a) {return cur;}
		} else if (!digits) {
			*context = LOC_SIMPLE;
			return start;
		}
		cur--; digits = 0;
		while (cur >= a && NATDIGIT(cur)) {cur--; digits++;}
		if (!digits) {*context = LOC_SIMPLE; return start;}
		if (cur < a) {return cur+1;}
	}
	if (*cur == '-' || *cur == '+') {
		if (*cur == '-') {*invert = 1;}
		if (*context == LOC_UNKOWN) {
			if (digits) {
				*context = LOC_SIGNNUM;
			} else {
				*context = LOC_SIGN;
			}
			if (cur == a) {return cur;}
		} else if (!digits) {
			*context = LOC_SIMPLE;
			return start;
		}
		cur--;
	}
	if (*context == LOC_UNKOWN) {*context = LOC_SIGNNUM;}
	if (cur < a) {return cur+1;}
	if (BLANK(cur)) {return cur+1;}
	*context = LOC_SIMPLE;
	return start;
}

char *naturalcompare_numbercontext_after(char *cur, int curlen, char *start, int nmleft, int *context) {
	if (*context == LOC_BEFORE) {
		if (!curlen || BLANK(cur)) {
			*context = LOC_NONUM; return cur-1;
		}
		if (*cur == '-' || *cur == '+') {
			*context = LOC_SIGN;
			cur++; curlen--;
			if (!curlen || BLANK(cur) || !NATDIGIT(cur)) {
				*context = LOC_NONUM; return cur-1;
			}
		} else if (!NATDIGIT(cur)) {
			*context = LOC_NONUM; return cur-1;
		}
		*context = LOC_SIGNNUM;
	}
	if (*context == LOC_SIGN) {
		if (!curlen || BLANK(cur) || !NATDIGIT(cur)) {
			*context = LOC_NONUM; return cur-1;
		}
		*context = LOC_SIGNNUM;
	}
	if (*context == LOC_SIGNNUM) {
		while (curlen && NATDIGIT(cur)) {cur++; curlen--;}
		if (!curlen || BLANK(cur)) {return cur;}
		if (*cur == '.') {
			cur++; curlen--;
			if (!curlen || BLANK(cur) || !NATDIGIT(cur)) {
				return cur;
			}
			*context = LOC_DECIMALNUM;
		} else if (*cur == 'e' || *cur == 'E') {
			cur++; curlen--;
			*context = LOC_E;
		}
	}
	if (*context == LOC_DECIMAL) {
		if (!curlen || BLANK(cur) || !NATDIGIT(cur)) {
			*context = LOC_SIGNNUM; return cur;
		}
		*context = LOC_DECIMALNUM;
	}
	if (*context == LOC_DECIMALNUM) {
		while (curlen && NATDIGIT(cur)) {cur++; curlen--;}
		if (!curlen || BLANK(cur)) {return cur;}
		if (*cur == 'e' || *cur == 'E') {
			cur++; curlen--;
			*context = LOC_E;
		} else {
			return cur;
		}
	}
	if (*context == LOC_E) {
		if (*cur == '-' || *cur == '+') {
			cur++; curlen--;
		}
		if (curlen <= 0 || !NATDIGIT(cur)) {
			*context = LOC_DECIMALNUM; return cur;
		}
		*context = LOC_ENUM;
	}
	if (*context == LOC_ESIGN) {
		if (!curlen || BLANK(cur) || !NATDIGIT(cur)) {
			*context = LOC_DECIMALNUM; return cur;
		}
		*context = LOC_ENUM;
	}
	while (curlen && NATDIGIT(cur)) {cur++; curlen--;}
	if (!curlen || BLANK(cur)) {return cur;}
	return cur;
}

int naturalcompare_atof(char *a,char *b) {
	double da = atof(a);
	double db = atof(b);
	 DPRINT("atof %s vs %s: %f, %f",a,b,da,db);
	if (da < db) {
		return -1;
	} else if (da > db) {
		return 1;
	} else {
		return 0;
	}
}

/*
	set list [list { } ! \" \# $ % & \' ( ) * + , - . / 0 1 2 3 4 5 6 7 8 9 : \; < = > ? @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z \[ \\ \] ^ _ ` a b c d e f g h i j k l m n o p q r s t u v w x y z \{ | \} ~]
	set rlist [list { } ! \" \# $ % & \' ( ) , - . / 0 1 2 3 4 5 6 7 8 9 + : \; < = > ? @ A a B b C c D d E e F f G g H h I i J j K k L l M m N n O o P p Q q R r S s T t U u V v W w X x Y y Z z \[ \\ \] ^ _ ` \{ | \} ~ *]

	set cor [lmath_calc [list_cor $rlist $list] + 32]
	set temp [list_concat [list_fill 32 0 1] $cor]
	join $temp ,
*/

/*
	changes in order from ascii
	* sorts after everything except control+space, then -
	+ sorts after numbers
	letters sort as as A a B b ...
	chars between upper and lowercase letters moved to after letters
*/
const int char_reorder[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,126,56,42,43,44,45,46,47,48,49,50,51,52,53,54,55,57,58,59,60,61,62,63,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,117,118,119,120,121,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,101,103,105,107,109,111,113,115,122,123,124,125};

int naturalcompare_diff(int left,int right) {
	if (left == right) {
		return 0;
	} else {
		if (left < 127) {left = char_reorder[left];}
		if (right < 127) {right = char_reorder[right];}
		return left - right;
	}
}

/*
 * The code below compares the numbers in the two
 * strings without ever converting them to numbers.  It
 * does this by first comparing the number of digits after the diff
 * if they are the same, diff determines result.
 */
int naturalcompare_simplenum(char *left,int alen,char *right, int blen,int diff, int secondaryDiff, int invert) {
	 DPRINT("  compare_simplenum: %s vs %s    diff=%d invert=%d secondaryDiff=%d", left, right, diff,invert,secondaryDiff);
	while (alen && blen) {
		if (!NATDIGIT(left)) {
			if (!NATDIGIT(right)) {
				if (!invert) {return (diff<0)?-1:1;} else {return (diff<0)?1:-1;}
			} else {
				if (!invert) {return -1;} else {return 1;}
			}
		} else if (!NATDIGIT(right)) {
			if (!invert) {return 1;} else {return -1;}
		}
		if (diff == 0) {
			diff = naturalcompare_diff(*left,*right);
		}
		left++; right++; alen--; blen--;
	}
	if (!alen || !NATDIGIT(left)) {
		if (!blen || !NATDIGIT(right)) {
			if (diff == 0) {diff = secondaryDiff;}
			if (!invert) {return (diff<0)?-1:1;} else {return (diff<0)?1:-1;}
		} else {
			if (!invert) {return -1;} else {return 1;}
		}
	} else if (!blen || !NATDIGIT(right)) {
		/* a at end of digits already tested in previous */
		if (!invert) {return 1;} else {return -1;}
	}
	if (diff == 0) {diff = secondaryDiff;}
	if (!invert) {return (diff<0)?-1:1;} else {return (diff<0)?1:-1;}
	
}

/* 
    "natural" sort order: deals with mixed alphabetic and numbers
    returns 0 if equal, -1 if a < b and 1 if a > b
 */
int naturalcompare(char const *a, char const *b,int alen,int blen) {
	int diff, nmleft,nmright,result,simplenum,context,contextleft,contextright;
	int secondaryDiff = 0,invert=0;
	char *left=NULL,*right=NULL,*start=NULL,*rstart=NULL;
	if (a == NULL) {
		if (b == NULL) {
			return 0;
		} else {
			return 1;
		}
	} else if (b == NULL) {
		return -1;
	}
	left = (char *)a;
	right = (char *)b;
	 DPRINT("---------- naturalcompare {%s} vs {%s} ----------",a,b);
	/* find the first difference (not case) */
	while (1) {
		diff = naturalcompare_diff(*left,*right);
		 DPRINT("diff %c vs %c: %d", *left, *right, diff);
		if (!alen || *left == '\0') {
			if (!blen || *right == '\0') {return secondaryDiff;} else {break;}
		}
		if (!blen || *right == '\0') {break;}
		if (diff != 0) {
			/* only sort on case if no other diff -> keep secondaryDiff for case diff */
			/* remember: in naturalcompare_diff (reordered) upper and lower case are 1 apart */
			if (diff == -1 && islower(UCHAR(*right))) {
				secondaryDiff = -1;
			} else if (diff == 1 && islower(UCHAR(*left))) {
				secondaryDiff = 1;
			} else {
				break;
			}
		}
		left++; alen--;
		right++; blen--;
	}
	if (diff == 0) {return secondaryDiff;}
	nmleft = ISNUMBER(left);
	nmright = ISNUMBER(right);
	if (!nmleft && !nmright) {
		/* for sure no number involved, lexical sort */
		return (diff<0)?-1:1;
	}
	/* 
		(potential) number comparison
		-----------------------------
	 */
	 DPRINT("compare number: %s vs %s    digit %d vs %d    diff %d", left, right, nmleft, nmright, diff);
	/* move back to start of number to get context */
	
	start = left;
	simplenum = 0; /* becomes 1 if there is something (double e, etc.) preventing a special number interpretation */
	if (start == a || BLANK(start-1)) {
		/* diff at start of number, no need to parse before */
		/* The following must sort before numbers */
		/* empty or blank allways sorts first */
		if (!alen || BLANK(left)) {return -1;}
		if (!blen || BLANK(right)) {return 1;}
		context = LOC_BEFORE;
	} else {
		/* check for context before diff */
		start = naturalcompare_numbercontext_before(a,start,&invert,&context);
	 	DPRINT("checked before: %s context=%d invert=%d",start,context,invert);
		if (context <= LOC_SIMPLE) {simplenum = 1;}
	}
	rstart = (char *)b+(start-a);
	if (!simplenum) {
		/* check for special number after diff */
		contextleft = context;
		contextright = context;
		if (context == LOC_ENUM || context == LOC_ESIGN || context == LOC_SIGN || context == LOC_DECIMAL) {
			/* we can stop at these combination -> no furter number */
			if (nmleft != NM_DIGIT) {nmleft = 0;}
			if (nmright != NM_DIGIT) {nmright = 0;}
		}
		naturalcompare_numbercontext_after(left, alen, start, nmleft, &contextleft);
 		DPRINT("checked left after: %s left=%s contextleft=%d",start,left,contextleft);
		naturalcompare_numbercontext_after(right, blen, rstart, nmright, &contextright);
	 	DPRINT("checked right after: %s right=%s contextright=%d",rstart,right,contextright);
		/* sort numbers before the rest */
		if (contextleft > LOC_SIGN && contextright < LOC_SIGNNUM) {
			return -1;
		} else if (contextright > LOC_SIGN && contextleft < LOC_SIGNNUM) {
			return 1;
		}

		if (contextleft == LOC_NONUM || contextright == LOC_NONUM) {
			return (diff<0)?-1:1;
		}
		if ((nmleft || nmright) && !simplenum && (contextleft == LOC_ENUM || contextright == LOC_ENUM)) {
			/* diff is on the e of scientific notation -> compare via conversion to numbers using atof */
			result = naturalcompare_atof(start,rstart);
			if (result == 0) {
				int leftsc = (contextleft == LOC_ENUM);
				int rightsc = (contextright == LOC_ENUM);
				if (leftsc && !rightsc) {
					diff = 1;
				} else if (!leftsc && rightsc) {
					diff = -1;
				} else if (*start == '+' && *rstart != '+') {
					diff = 1;
				} else if (*start != '+' && *rstart == '+') {
					diff = -1;
				}
				if (!invert) {return (diff<0)?-1:1;} else {return (diff<0)?1:-1;}
			} else {
				return result;
			}
		}
		/* another few corner cases where the diff is not actually in the number anymore */
		if ((nmleft == NM_E || context == LOC_E || context == LOC_ESIGN) && contextleft != LOC_ENUM) {nmleft = 0;}
		if ((nmright == NM_E || context == LOC_E || context == LOC_ESIGN) && contextright != LOC_ENUM) {nmright = 0;}
		/* nmleft and nmright could have been changed to indicate diff is not really in a number */
		if (!nmleft && !nmright) {return (diff<0)?-1:1;}
	}
	if (simplenum) {
		 DPRINT("simpelnum - context %d diff=%d nmleft=%d nmright=%d", context,diff,nmleft,nmright);
		/* simplenum indicates we have embedded numbers */
		/* so we will just compare size of integer */
		/* we will not take into account invert, decimal */
		invert = 0;
		if ((left == a || !NATDIGIT(left-1)) && (nmleft != NM_DIGIT || nmright != NM_DIGIT)) {
			return (diff<0)?-1:1;
		} else {
			return naturalcompare_simplenum(left,alen,right,blen,diff,secondaryDiff,invert);
		}
	} else {
		 DPRINT("specialnum - invert %d - context %d", invert, context);
		/* special cases, +, -, . at diff */
		/* if not scientific notation (checked before) we can only have - or + at diff if at start */
		/* - always smaller */
		if (start == left) {
			 DPRINT("specialnum start: %s vs %s    diff %d", left, right, diff);
			/* difference is at start of number */
			if (*left == '-') {
				return -1;
			} else if (*right == '-') {
				return 1;
			} else if (*left == '+') {
				left++ ; alen--;
				if (alen == 0) {return 1;}
				result=naturalcompare(left,right,alen,blen);
				if (result == 0) {
					return 1;
				} else {
					return result;
				}
			} else if (*right == '+') {
				right++ ; blen--;
				if (blen == 0) {return -1;}
				result=naturalcompare(left,right,alen,blen);
				if (result == 0) {
					return -1;
				} else {
					return result;
				}
			} else if (*left == '.') {
				if (nmright == NM_DIGIT) {
					/* if the decimal point is in the dif and the other is a number, 
					   diff works, as . < number */
					if (!invert) {return (diff<0)?-1:1;} else {return (diff<0)?1:-1;}
				} else {
					/* if left is . its number is longer, and left > right */
					if (!invert) {return 1;} else {return -1;}
				}
			} else if (*right == '.') {
				if (nmleft == NM_DIGIT) {
					/* if the decimal point is in the dif and the other is a number, 
					   diff works, as . < number */
					if (!invert) {return (diff<0)?-1:1;} else {return (diff<0)?1:-1;}
				} else {
					/* if right is . its number is longer, and left < right */
					if (!invert) {return -1;} else {return 1;}
				}
			} else if (nmleft != NM_DIGIT || nmright != NM_DIGIT) {
				/* if one not a digit, diff decides (both not digit was tested earlier) */
				if (!invert) {return (diff<0)?-1:1;} else {return (diff<0)?1:-1;}
			} else {
				return naturalcompare_simplenum(left,alen,right,blen,diff,secondaryDiff,invert);
			}
		} else {
			 DPRINT("specialnum: %s vs %s  context=%d nmleft=%d nmright=%d diff %d", left, right, context, nmleft, nmright, diff);
			if (context < LOC_SIGNNUM) {
				if (nmleft != NM_DIGIT || nmright != NM_DIGIT) {
					/* like at start; if we are not yet in a number, use diff */
					return (diff<0)?-1:1;
				}
			} else {
				if (!nmleft) {
					if (!invert) {return -1;} else {return 1;}
				}
				if (!nmright) {
					if (!invert) {return 1;} else {return -1;}
				}
			}
			if (context == LOC_DECIMAL || context == LOC_DECIMALNUM) {
				/* after the decimal point, first diff decides */
				if (!invert) {return (diff<0)?-1:1;} else {return (diff<0)?1:-1;}
			}
			return naturalcompare_simplenum(left,alen,right,blen,diff,secondaryDiff,invert);
		}
	}
	fprintf(stderr,"internal error,we should not get here");
	exit(1);
}

int DStringCompare(DString *a, DString *b) {
	if (a == b) {return 0;}
	if (a == NULL) {
		return 1;
	} else if (b == NULL) {
		return -1;
	}
	return naturalcompare(a->string,b->string,a->size,b->size);
}

int loccompare (char const *a, char const *b,int alen,int blen) {
	if (alen >= 3) {
		if ((a[0] == 'C' || a[0] == 'c') && (a[1] == 'H' || a[1] == 'h') && (a[2] == 'R' || a[2] == 'r')) {
			a += 3; alen -= 3;
			if (alen && a[0] == '-') {
				a++; alen--;
			}
		}
	}
	if (blen >= 3) {
		if ((b[0] == 'C' || b[0] == 'c') && (b[1] == 'H' || b[1] == 'h') && (b[2] == 'R' || b[2] == 'r')) {
			b += 3; blen -= 3;
			if (blen && b[0] == '-') {
				b++; blen--;
			}
		}
	}
	return naturalcompare(a,b,alen,blen);
}

int chrmatch (char const *a, char const *b,int alen,int blen) {
	if (alen >= 3) {
		if ((a[0] == 'C' || a[0] == 'c') && (a[1] == 'H' || a[1] == 'h') && (a[2] == 'R' || a[2] == 'r')) {
			a += 3; alen -= 3;
			if (alen && a[0] == '-') {
				a++; alen--;
			}
		}
	}
	if (blen >= 3) {
		if ((b[0] == 'C' || b[0] == 'c') && (b[1] == 'H' || b[1] == 'h') && (b[2] == 'R' || b[2] == 'r')) {
			b += 3; blen -= 3;
			if (blen && b[0] == '-') {
				b++; blen--;
			}
		}
	}
	return naturalcompare(a,b,alen,blen);
}

int DStringLocCompare(DString *a, DString *b) {
	if (a == b) {return 0;}
	if (a == NULL) {
		return 1;
	} else if (b == NULL) {
		return -1;
	}
	return loccompare(a->string,b->string,a->size,b->size);
}

char *Loc_ChrString(DString *ds) {
	char *s;
	if (ds == NULL) {
		return (char *)"";
	} else {
		s = ds->string;
		if (ds->size > 3) {
			if ((s[0] == 'C' || s[0] == 'c') && (s[1] == 'H' || s[1] == 'h') && (s[2] == 'R' || s[2] == 'r')) {
				if (ds->size > 4 && s[3] == '-') {
					s += 4;
				} else {
					s += 3;
				}
			}
		}
		return s;
	}
}

void DStringAppend(DString *dstring, char *string) {
	int size = strlen(string);
	int nsize = dstring->size + size;
	DStringSetSize(dstring,nsize);
	strncpy(dstring->string+dstring->size,string,size+1);
	dstring->size = nsize;
}

void DStringAppendS(DString *dstring, char *string, int size) {
	int nsize = dstring->size + size;
	DStringSetSize(dstring,nsize);
	strncpy(dstring->string+dstring->size,string,size+1);
	dstring->string[nsize] = '\0';
	dstring->size = nsize;
}

int DStringGetLine(DString *linePtr,	FILE *f1) {
	register char *buf,*buf_end;
	int bsz;
	register char *cur;
	register int cnt;
	register int	c;
	if (linePtr->memsize != -1 && linePtr->string != linePtr->staticspace) {
		buf = linePtr->string;
		buf_end = linePtr->string+linePtr->memsize;
		bsz = linePtr->memsize;
	} else {
		buf=(char *)malloc(128);
		bsz=128;
		buf_end=buf+bsz;
	}
	cur=buf;
	cnt=0;
	while ((c=getc_unlocked(f1))!=EOF) {
		if (c == '\n') break;
		*cur++=c;
		cnt++;
		if (cur == buf_end) {
			buf=(char *)realloc(buf,bsz*2);
			cur=buf+bsz;
			bsz*=2;
			buf_end=buf+bsz;
		}
	}
	*cur='\0';
	linePtr->string = buf;
	linePtr->size = cnt;
	linePtr->memsize = buf_end-buf;
	if (c == EOF && cnt == 0) return -1;
	return cnt;
}

void SkipLine(FILE *f1) {
	register int c;
	while ((c=getc_unlocked(f1))!=EOF) {
		if (c == '\n') break;
	}
}

void InitBuffer(Buffer *buffer, int size) {
	buffer->memsize = size;
	buffer->data = (char *)malloc(size+1);
	buffer->pos = buffer->data;
	buffer->size = 0;
}

void DelBuffer(Buffer *buffer) {
	buffer->memsize = 0;
	free(buffer->data);
	buffer->data = NULL;
	buffer->pos = NULL;
	buffer->size = 0;
}

int DStringGetLine_b(DString *linePtr,	FILE *f1,Buffer *buffer) {
	char *cur = linePtr->string;
	int c,newdata=0;
	ssize_t size = 0;
NODPRINT("getline");
	while (1) {
		if (!buffer->size) {
			buffer->size = fread(buffer->data,sizeof(char),buffer->memsize,f1);
			buffer->pos = buffer->data;
NODPRINT("read buffer %d",buffer->size);
			if (!buffer->size) break;
		}
		newdata = 1;
		buffer->size--;
		c = *(buffer->pos)++;
		if  (c == '\n') {
			break;
		}
		if (linePtr->memsize <= size) {
			linePtr->size = size;
			DStringSetSize(linePtr,2*linePtr->memsize);
			cur = linePtr->string+size;
		}
		*cur++ = c;
		++size;
	}
	if (!newdata) {
		*cur = '\0';
		return -1;
	}
	linePtr->size = size;
	if (linePtr->memsize <= size) {
		DStringSetSize(linePtr,2*linePtr->memsize);
		cur = linePtr->string+size;
	}
	*cur = '\0';
	return size;
}

DStringArray *DStringArrayNew(int size) {
	DStringArray *dstringarray;
	int i;
	dstringarray = (DStringArray *)malloc(1*sizeof(DStringArray));
	dstringarray->data = (DString *)malloc(size*sizeof(DString));
	dstringarray->size = 0;
	dstringarray->memsize = size;
	dstringarray->datablock = NULL;
	for (i = 0; i < size ; i++) {
		DStringInit(dstringarray->data+i);
	}
	return dstringarray;
}

DStringArray *DStringArrayAppend(DStringArray *dstringarray,char *string,int size) {
	if (dstringarray->size == dstringarray->memsize) {
		DString *temp = NULL;
		int oldmemsize = dstringarray->memsize,i;
		dstringarray->memsize *= 2;
		/* 
			need to correct string pointer in all DStrings that point to the static space
		*/
		temp = (DString *)malloc(dstringarray->memsize*sizeof(DString));
		memcpy(temp,dstringarray->data,oldmemsize*sizeof(DString));
		for (i = 0 ; i < dstringarray->size ; i++) {
			if (dstringarray->data[i].string == dstringarray->data[i].staticspace) {
				temp[i].string = temp[i].staticspace;
			}
		}
		free(dstringarray->data);
		dstringarray->data = temp;
	}
	DStringInit(dstringarray->data+dstringarray->size);
	if (size < 0) {size = strlen(string);}
	DStringSetS(dstringarray->data+dstringarray->size, string, size);
	dstringarray->size++;
	return dstringarray;
}

DStringArray *DStringArrayFromChar(char *string,char sep) {
	DStringArray *result;
	char *cur,*prev;
	int count=1;
	cur=string;
	while(*cur) {
		if (*cur == sep) count++;
		cur++;
	}
	result = DStringArrayNew(count+1);
	cur = string;
	prev = cur;
	while(*cur) {
		if (*cur == sep) {
			DStringArrayAppend(result,prev,cur-prev);
			prev = cur+1;
		}
		cur++;
	}
	DStringArrayAppend(result,prev,cur-prev);
	return result;
}

DStringArray *DStringArrayFromCharM(char *string,char *seps) {
	DStringArray *result;
	char *cur,*prev;
	result = DStringArrayNew(5);
	cur = string;
	prev = cur;
	while(*cur) {
		if (strchr((const char *)seps,*cur)) {
			DStringArrayAppend(result,prev,cur-prev);
			do {cur++;} while(strchr((const char *)seps,*cur) && *cur);
			prev = cur;
			if (!*cur) break;
		} else {
			cur++;
		}
	}
	if (cur-prev > 0) {
		DStringArrayAppend(result,prev,cur-prev);
	}
	return result;
}

DStringArray *DStringArrayRange(DStringArray *dstringarray,int start, int end) {
	DStringArray *result;
	int i;
	result = DStringArrayNew(end-start+1);
	for (i = start; i <= end ; i++) {
		DStringCopy(result->data+(i-start), dstringarray->data+i);
	}
	result->size = end-start+1;
	return result;	
}

DStringArray *DStringArraySet(DStringArray *dstringarray,int pos,char *string,int size) {
	int i;
	if (pos >= dstringarray->memsize) {
		dstringarray->memsize = pos+1;
		dstringarray->data = (DString *)realloc(dstringarray->data,dstringarray->memsize*sizeof(DString));
		i = dstringarray->size;
		while (i <= pos) {
			DStringInit(dstringarray->data+i);
			i++;
		}
		dstringarray->size=pos+1;
	}
	if (size < 0) {size = strlen(string);}
	DStringSetS(dstringarray->data+pos, string, size);
	return dstringarray;
}

int DStringArraySearch(DStringArray *dstringarray,char *string,int size) {
	DString *astring;
	int i;
	for (i = 0; i < dstringarray->size ; i++) {
		astring = DStringArrayGet(dstringarray,i);
		if (size != -1) {
			if (astring->size == size && strncmp(astring->string,string,size) == 0) {
				return i;
			}
		} else {
			if (strcmp(astring->string,string) == 0) {
				return i;
			}
		}
	}
	return -1;
}

void DStringArrayDestroy(DStringArray *dstringarray) {
	if (dstringarray == NULL) return;
	int i=0;
	for (i =0; i < dstringarray->size ; i++) {
		if (dstringarray->data[i].string != NULL) {
			DStringClear(dstringarray->data+i);
		}
	}
	free(dstringarray->data);
	if (dstringarray->datablock) DStringDestroy(dstringarray->datablock);
	free(dstringarray);
}

void check_numfieldserror(int numfields,int numfields2,DString *line,char *filename,unsigned int *linenum) {
	if (numfields != numfields2) {
		fprintf(stderr,"ERROR in file %s",filename);
		if (linenum != NULL) {
			fprintf(stderr," at line %d",*linenum);
		}
		fprintf(stderr,": number of columns different from the header: %d instead of %d for:\n%s\n",numfields,numfields2,line->string);
		exit(1);
	}
}

int DStringGetTab(
	DString *linePtr,	FILE *f1, int maxtab, DStringArray *result, int setzero,unsigned int *numfields
) {
	register char *cur = linePtr->string;
	register int c,newdata=0;
	register unsigned int count=0,othertab=0;
	ssize_t size = 0;
NODPRINT("maxtab=%d result->memsize=%d",maxtab,result->memsize)
	maxtab += 1;
	if (maxtab > result->memsize) {
		fprintf(stderr,"cannot DStringGetTab, size of allocated DStringArray must be >= maxtab+2\n");
		exit(1);
	}
	result->data[count].string = cur;
	result->data[count].memsize = -1;
	/* fill current pos (size) in the array for now, convert to array element sizes later */
	while (1) {
		c = getc_unlocked(f1);
		if  (c == '\t') {
			if (count < maxtab) {
				result->data[count].size = size;
				/* fprintf(stdout,"count=%d size=%d\n",count, result->data[count].size); */
				count++;
			} else {
				othertab++;
			}
		} else if (c == '\n' || c == EOF) {
			if (count < maxtab) {
				result->data[count].size = size;
			}
			break;
		}
		newdata = 1;
		if (linePtr->memsize <= (size+1)) {
			linePtr->size = size;
			DStringSetSize(linePtr,2*linePtr->memsize);
			cur = linePtr->string+size;
		}
		*cur++ = c;
		++size;
	}
	*cur = '\0';
	if (numfields != NULL) {
		*numfields = count+othertab+1;
	}
	if (count < maxtab) {
		count++;
		result->size = count;
	} else {
		result->size = maxtab;
	}
	/* add \0 to line*/
	linePtr->size = size;
	if (linePtr->memsize <= size) {
		DStringSetSize(linePtr,linePtr->memsize+2);
		cur = linePtr->string+size;
	}
	*cur = '\0';
	if (!newdata) {
		return -1;
	}
	/* fill rest of tab separated elements in array */
	while (count <= maxtab) {
		result->data[count++].size = size;
	}
	/* make array */
	{
	register int prevsize;
	count = 0;
	prevsize = 0;
	while (count <= maxtab) {
		int pos = result->data[count].size;
		result->data[count].size = result->data[count].size-prevsize;
		result->data[count].string = linePtr->string + prevsize;
		if (setzero) {
			result->data[count].string[result->data[count].size] = '\0';
		}
		result->data[count].memsize = -1;
NODPRINT("final count=%d size=%d %s",count, result->data[count].size,result->data[count].string)
		if (pos < linePtr->size) {prevsize = pos+1;} else {prevsize = linePtr->size;}
		count++;
	}
	}
	NODPRINT("\n==== result ====\n%d %s",size,linePtr->string);
	NODPRINT("==================== getline finished ====================");
	return 0;
}

void DStringPrintTab(FILE *f, DString *linePtr) {
	int i;
	for(i = 0 ; i < linePtr->size ; i++) {
		if (linePtr->string[i] == '\0') {fputc('\t',stderr);} else {fputc(linePtr->string[i],stderr);}
	}
}

/* 
 * Creates a Tab from a DString
 * linePtr: DString to be converted
 * maxtab: maximum number that can be converted (size of DStringArray)
 * result: resulting DStringArray
 * setzero: if 1, the tabs in the string are changed to \0, 
 *          making the array strings compatible with C strings (but changes the original string)
 * numfields: returns actual number of fields found (NULL if not needed)
 */
int DStringSplitTab(
	DString *linePtr,	int maxtab, DStringArray *result, int setzero,unsigned int *numfields
) {
	register char *cur = linePtr->string;
	register int c;
	register unsigned int count=0,othertab=0;
	ssize_t size = 0;
NODPRINT("maxtab=%d result->memsize=%d",maxtab,result->memsize)
	maxtab += 1;
	if (maxtab > result->memsize) {
		fprintf(stderr,"cannot DStringGetTab, array memsize < maxtab+1\n");
		exit(1);
	}
	result->data[count].string = cur;
	result->data[count].memsize = -1;
	/* fill current pos (size) in the array for now, convert to array element sizes later */
	while (1) {
		c = *cur;
		if  (c == '\t') {
			if (count < maxtab) {
				result->data[count].size = size;
				/* fprintf(stdout,"count=%d size=%d\n",count, result->data[count].size); */
				count++;
			} else {
				othertab++;
			}
		} else if (size >= linePtr->size) {
			if (count < maxtab) {
				result->data[count].size = size;
			}
			break;
		}
		cur++;
		++size;
	}
	if (numfields != NULL) {
		*numfields = count+othertab+1;
	}
	if (count < maxtab) {
		count++;
		result->size = count;
	} else {
		result->size = maxtab;
	}
	/* fill rest of tab separated elements in array */
	while (count <= maxtab) {
		result->data[count++].size = size;
	}
	/* make array */
	{
	register int prevsize;
	count = 0;
	prevsize = 0;
	while (count <= maxtab) {
		int pos = result->data[count].size;
		result->data[count].size = result->data[count].size-prevsize;
		result->data[count].string = linePtr->string + prevsize;
		if (setzero) {
			result->data[count].string[result->data[count].size] = '\0';
		}
		result->data[count].memsize = -1;
NODPRINT("final count=%d size=%d %s",count, result->data[count].size,result->data[count].string)
		prevsize = pos+1;
		count++;
	}
	}
	NODPRINT("\n==== result ====\n%d %s",size,linePtr->string);
	NODPRINT("==================== getline finished ====================");
	return 0;
}

int parse_pos(char *arg, int **rresult, int *rnum) {
	char *pch;
	int *result;
	int num,memsize = 5;
	result = (int *)malloc(memsize*sizeof(int));
	pch = strtok (arg, " \t,-");
	num = 0;
	while (pch != NULL) {
		if (*pch == 'x') {
			result[num++] = -1;
		} else {
			result[num++] = atoi(pch);
		}
		if (num >= memsize) {
			memsize += memsize;
			result = (int *)realloc(result, memsize*sizeof(int));
		}
		pch = strtok (NULL, " \t,-");
	}
	*rnum = num;
	*rresult = result;
	return 0;
}

int get_region(FILE *f1, DString *linePtr, int chr1pos, int start1pos, int end1pos, int max1, DString **chromosome1, int *start1, int *end1) {
	char *linepos = NULL,*scanpos = NULL,*endpos;
	ssize_t read;
	int count;
	if (chromosome1 == NULL) {
		return 1;
	}
	if (f1 == NULL) {
		DStringDestroy(*chromosome1);
		*chromosome1 = NULL;
		return 1;
	}
	while ((read = DStringGetLine(linePtr, f1)) != -1) {
NODPRINT("%s",linePtr->string)
		linepos = linePtr->string;
		scanpos = linePtr->string;
		endpos = linePtr->string+linePtr->size;
		count = 0;
		while (1) {
			if (*linepos == '\t' || *linepos == '\0' || linepos == endpos) {
				if (count == chr1pos) {
					char *pos=scanpos;
					while (*pos != '\t' && *pos != '\0' && pos != endpos) {++pos;}
					DStringSetS(*chromosome1,scanpos,pos-scanpos);
				} else if (count == start1pos) {
					sscanf(scanpos,"%d",start1);
				} else if (count == end1pos) {
					sscanf(scanpos,"%d",end1);
				}
				if (count >= max1) break;
				count++;
				scanpos = linepos+1;
			}
			if (*linepos == '\0') break;
			if (linepos == endpos) break;
			linepos++;
		}
		if (count >= max1) break;
	}
	if (read == -1) {
		DStringDestroy(*chromosome1);
		*chromosome1 = NULL;
		return 1;
	} else {
		return 0;
	}
}

/*
 * this skips the header (file position will be just after the header line),
 * and returns the header line in linePtr.
 * if numfields != NUL, it will contain the number of fields in the header
 */
void skip_header(FILE *f1, DString *linePtr, unsigned int *numfields,unsigned int *pos) {
	ssize_t read;
	unsigned int curpos=0;
	if (pos != NULL) {curpos=*pos;}
	read = DStringGetLine(linePtr, f1);
	if (read == -1) return;
	if (strlen(linePtr->string) >= 16 && strncmp(linePtr->string,"##fileformat=VCF",16) == 0) {
		/* vcf style header */
		while (read != -1) {
			if (linePtr->string[0] == '\0') {
				DStringGetLine(linePtr, f1); curpos++;
				break;
			}
			if (linePtr->string[0] != '#' || linePtr->string[1] != '#') {
				break;
			}
			read = DStringGetLine(linePtr, f1); curpos++;
		}
	} else {
		while (read != -1) {
			if (linePtr->string[0] == '\0') {
				DStringGetLine(linePtr, f1); curpos++;
				break;
			}
			if (linePtr->string[0] != '#' && linePtr->string[0] != '>') break;
			read = DStringGetLine(linePtr, f1); curpos++;
		}
	}
	if (pos != NULL) {*pos = curpos;}
	if (numfields != NULL) {
		char *buffer;
		unsigned int count;
		buffer = linePtr->string;
		count = 1;
		while (*buffer != '\0') {
			if (*buffer == '\t') {count++;}
			buffer++;
		}
		*numfields = count;
	}
	NODPRINT("%s",linePtr->string)
}

FILE *fopen64_or_die(char *filename,char *mode) {
	FILE *f;
	f = fopen64(filename,mode);
	if (f == NULL) {
		fprintf(stderr,"Error opening file %s: %s.\n", filename, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return(f);
}

int checksort(DString *prevchromosome1,int *prevstart1,int *prevend1,DString *prevtype1,DString *prevalt1,DString *chromosome1,int start1,int end1,DString *type1,DString *alt1,char *filename,int *nextpos,int fillprev) {
	int comp,comptype,compalt, result = 2;
 	comp = DStringLocCompare(chromosome1,prevchromosome1);
	comptype = DStringCompare(type1,prevtype1);
	compalt = DStringCompare(alt1,prevalt1);
	if (comp < 0 || (comp == 0 && 
		(start1 < *prevstart1 || (start1 == *prevstart1 && 
		(end1 < *prevend1 || (end1 == *prevend1 &&
		(comptype < 0 || (comptype == 0 && compalt < 0)
	))))))) {
		fprintf(stderr,"File (%s) is not correctly sorted (sort correctly using \"cg select -s -\")\n",filename);
		fprintf(stderr,"%s:%d-%d:%s:%s came before %s:%d-%d:%s:%s\n",prevchromosome1->string,*prevstart1,*prevend1,prevtype1->string,prevalt1->string, chromosome1->string,start1,end1,type1->string,alt1->string);
		exit(1);
	} else if (comp > 0) {
		/* prevchromosome1 = chromosome1; */
		if (fillprev) {DStringCopy(prevchromosome1,chromosome1);}
		if (nextpos) {*nextpos = 0;}
	}
	if (comp == 0 && comptype == 0 && start1 == *prevstart1 && end1 == *prevend1) {
		if (compalt == 0) {
			result = 0;
		} else {
			result = 1;
		}
	}
	if (fillprev) {
		*prevstart1 = start1; *prevend1 = end1;
		if (comptype != 0) {DStringCopy(prevtype1,type1);}
		if (compalt != 0) {DStringCopy(prevalt1,alt1);}
	}
	return result;
}

/* checksortreg(prevchromosome1,prevstart1,prevend1,chromosome1,start1,end1,argv[1]); */
int checksortreg(DString *prevchromosome,int *prevstart,int *prevend,DString *chromosome,int start,int end,char *file) {
	int comp;
	if (!prevchromosome || !chromosome) {
		return 0;
	}
 	comp = DStringLocCompare(chromosome,prevchromosome);
	if (comp < 0 || (comp == 0 && (start < *prevstart || (start == *prevstart && end < *prevend)))) {
		fprintf(stderr,"File (%s) is not correctly sorted (sort correctly using \"cg select -s -\")\n",file);
		fprintf(stderr,"%*.*s:%d-%d came before %*.*s:%d-%d\n",
			prevchromosome->size,prevchromosome->size,prevchromosome->string,
			*prevstart,*prevend,
			chromosome->size,chromosome->size,chromosome->string,
			start,end
		);
		exit(1);
	}
	*prevstart = start; *prevend = end;
	if (comp > 0) {
		DStringCopy(prevchromosome,chromosome);
		return 1;
	} else {
		return 0;
	}
}

int fileexists(const char * filename) {
	/* try to open file to read */
	FILE *file = fopen(filename, "r");
	if (file) {
		FCLOSE(file);
		return 1;
	}
	return 0;
}

char *tempfilename() {
	int fd;
	char *tmpdir,*name;
	if ((tmpdir = getenv ("TEMP")) == NULL) {
		if ((tmpdir = getenv ("TMP")) == NULL) {
			if ((tmpdir = getenv ("TMPDIR")) == NULL) {
				tmpdir = "/tmp";
			}
		}
	}
	name = malloc(strlen(tmpdir)+15);
	sprintf(name,"%s/tmpfileXXXXXX",tmpdir);
	fd = mkstemp(name);
	close(fd);
	return name;
}



#endif

