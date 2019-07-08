#include <stdio.h>

#include "easel.h"

#define RAISE_EWRITE ESL_EXCEPTION_SYS(eslEWRITE, "write failed")
#define XRAISE_EWRITE ESL_XEXCEPTION_SYS(eslEWRITE, "write failed")

#define FPRINTF(stream, fmt, ...)                                                           \
  if (fprintf(stream, fmt, __VA_ARGS__) < 0) RAISE_EWRITE
#define PRINTF(format, ...) FPRINTF(stdout, format, __VA_ARGS__)
#define FPUTS(str, stream)                                                             \
  if (fputs(str, stream) < 0) RAISE_EWRITE
#define PUTS(str) FPUTS(str, stdout)

#define XFPRINTF(stream, fmt, ...)                                                     \
  if (fprintf(stream, fmt, __VA_ARGS__) < 0) XRAISE_EWRITE
#define XPRINTF(format, ...) XFPRINTF(stdout, format, __VA_ARGS__)
#define XFPUTS(str, stream)                                                            \
  if (fputs(str, stream) < 0) XRAISE_EWRITE
#define XPUTS(str) XFPUTS(str, stdout)
