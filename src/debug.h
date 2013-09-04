/* #define DEBUG 1 */

#ifdef DEBUG
#define DPRINT(...) fprintf(stderr, "-- ");fprintf(stderr, __VA_ARGS__);fprintf(stderr, "  (%s:%d)\n", __FUNCTION__, __LINE__);fflush(stderr);fflush(stdout);
#define NODPRINT(fmt, ...) /**/
#else
#define DPRINT(fmt, ...) /**/
#define NODPRINT(fmt, ...) /**/
#endif
