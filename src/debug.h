/* #define DEBUG 1 */

#ifdef DEBUG
#define DPRINT(...) fprintf(stderr, "-- %s:%d: ", __FUNCTION__, __LINE__);fprintf(stderr, __VA_ARGS__);fprintf(stderr, " --\n");
#define NODPRINT(fmt, ...) /**/
#else
#define DPRINT(fmt, ...) /**/
#define NODPRINT(fmt, ...) /**/
#endif
