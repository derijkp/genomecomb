#include "tcl.h"
#include "genomecomb.h"
#include <sys/types.h>
#include <time.h>
#include <math.h>

extern int genomecomb_tsv_select_ObjCmd _ANSI_ARGS_((ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]));

extern int genomecomb_tsv_select_indexed_ObjCmd _ANSI_ARGS_((ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]));

extern int genomecomb_loc_compare_ObjCmd _ANSI_ARGS_((ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]));

extern int genomecomb_nat_compare_ObjCmd _ANSI_ARGS_((ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]));

extern int genomecomb_annotategene_findregc_ObjCmd _ANSI_ARGS_((ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]));

extern int genomecomb_bsortObjCmd _ANSI_ARGS_((ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]));

extern int genomecomb_alignedseqObjCmd _ANSI_ARGS_((ClientData clientData,
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
	Tcl_CreateObjCommand(interp,"tsv_selectc",(Tcl_ObjCmdProc *)genomecomb_tsv_select_ObjCmd,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateObjCommand(interp,"tsv_selectc_indexed",(Tcl_ObjCmdProc *)genomecomb_tsv_select_indexed_ObjCmd,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateObjCommand(interp,"loc_compare",(Tcl_ObjCmdProc *)genomecomb_loc_compare_ObjCmd,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateObjCommand(interp,"nat_compare",(Tcl_ObjCmdProc *)genomecomb_nat_compare_ObjCmd,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateObjCommand(interp,"annotategene_findregc",(Tcl_ObjCmdProc *)genomecomb_annotategene_findregc_ObjCmd,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateObjCommand(interp,"bsort",(Tcl_ObjCmdProc *)genomecomb_bsortObjCmd,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
	Tcl_CreateObjCommand(interp,"alignedseq",(Tcl_ObjCmdProc *)genomecomb_alignedseqObjCmd,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
	return TCL_OK;
}
