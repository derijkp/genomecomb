
#include "razf.h"

typedef struct BCol {
	char *file;
	char *type;
	int reverse;
	int isunsigned;
	int typesize;
	char *buffer;
	int buffersize;
	RAZF *rz;
	FILE *f;
} BCol;

int bcol_NeedReversing(int format);
void bcol_CopyNumber(const void *from,	void *to,unsigned int length,	int reverse);

int bcol_get(BCol *fbcol,int start,int end);
int bcol_printbin(FILE *f,int reverse,int isunsigned,char *type,char *string);
int bcol_printtext(FILE *f,int reverse,int isunsigned,char *type,char *bin);

BCol *bcol_open(char *bcolfile,char *bcoltype);
void bcol_close(BCol *bcol);

