
#include "razf.h"

typedef struct BCol {
	char *file;
	char type;
	int reverse;
	int isunsigned;
	int typesize;
	int start;
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

BCol *bcol_open(char *bcolfile);
void bcol_close(BCol *bcol);

