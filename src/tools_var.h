#ifndef TOOLS_VAR_H_LOADED
#define TOOLS_VAR_H_LOADED 1

#include "gztools.h"

typedef struct VariantPos {
	int chr;
	int start;
	int end;
	int type;
	int ref;
	int alt;
	int seq;
	int zyg;
	int a1;
	int a2;
	int max;
	int id;
} VariantPos;

typedef struct Variant {
	DString *chr;
	int start;
	int end;
	DString *type;
	DString *ref;
	DString *alt;
	DString *a1;
	DString *a2;
	int id;
} Variant;

typedef struct VarFile {
	char *file;
	GZFILE *f;
	DString *headerline;
	DString *prevline;
	DString *line;
	DStringArray *header;
	DStringArray *prevresult;
	DStringArray *result;
	VariantPos varpos;
	Variant *var;
	Variant *prevvar;
	int split;
	int max;
	int linenr;
	unsigned int numfields;
	unsigned int pos;
	unsigned int error;
} VarFile;

/* returns 0 if equal, 1 if equal not including alt, 2 if different */
int varchecksort(Variant *prev,Variant *var,char *filename,int *nextpos);
void varputs(Variant var,FILE *f);
void varputs_chr(DString *chr,FILE *f);
void result2var(DStringArray *result,VariantPos varpos, Variant *var);
int varcompare(Variant *var1, Variant *var2, int split);
int regcompare(Variant *var1, Variant *var2);
void varpos_init(VariantPos *varpos);
void var_init(Variant *var);
int varpos_max(VariantPos *varpos);
void varpos_fromheader(VariantPos *varpos,DStringArray *header);
VarFile *OpenVarfile(char *filename, int split);
void Varfile_checkbasicfields(VarFile *varfile);
void Varfile_checkbasicfields_basic(VarFile *varfile);
Variant *varfile_next(VarFile *varfile,int check);
void CloseVarfile(VarFile *varfile);

#endif
