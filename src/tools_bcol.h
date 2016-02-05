
#include "razf.h"

typedef struct BCol_table {
	DString *chr;
	unsigned long long begin;
	unsigned long long end;
	unsigned long long pos;
} BCol_table;

typedef struct BCol {
	char *file;
	char type;
	int reverse;
	int isunsigned;
	int typesize;
	int start;
	int version;
	BCol_table *table;
	int tablesize;
	unsigned char *buffer;
	int buffersize;
	RAZF *rz;
	FILE *f;
} BCol;

int bcol_NeedReversing(char format);
void bcol_CopyNumber(const void *from,	void *to,unsigned int length,	int reverse);

int bcol_getbin(BCol *fbcol,int start,int end);
int bcol_printbin(FILE *f,int reverse,int isunsigned,char *type,char *string);
int bcol_printtext(FILE *f,int reverse,int isunsigned,char type,unsigned char *bin);
int bcol_readdouble(BCol *fbcol,long double *result);

BCol *bcol_open(char *bcolfile);
void bcol_close(BCol *bcol);

