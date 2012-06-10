#include "tcl.h"
#include "genomecomb.h"
#include <sys/types.h>
#include <time.h>
#include <math.h>

extern int genomecomb_tsv_select_ObjCmd _ANSI_ARGS_((ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]));

int
Genomecomb_Init(interp)
	Tcl_Interp *interp;		/* Interpreter to add extra commands */
{
#ifdef USE_TCL_STUBS
	if (Tcl_InitStubs(interp, "8.1", 0) == NULL) {
		return TCL_ERROR;
	}
#endif
	Tcl_CreateObjCommand(interp,"tsv_select",(Tcl_ObjCmdProc *)genomecomb_tsv_select_ObjCmd,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
	return TCL_OK;
}
