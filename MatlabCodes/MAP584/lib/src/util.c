#include <stdio.h>
#include "stdarg.h"
#include <stdlib.h>

#include "util.h"

#ifdef _DEBUG
static FILE *f_debug = NULL;

void debug_printf(char *format, ...)
{
  va_list ap;

  va_start(ap, format);

  if (f_debug) {
    vfprintf(f_debug, format, ap);
  }

  va_end(ap);
}

void debug_begin(const char *fname)
{
  if ((fname == NULL) || (*fname == '\0'))
    return;

  debug_end();
  
  f_debug = fopen(fname, "w");
}

void debug_end()
{
  if (f_debug) {
    fclose(f_debug); 
    f_debug= NULL;
  }
}

#else

void debug_begin(const char *s) {};
void debug_end() {};
void debug_printf(char *format, ...) {}

#endif
