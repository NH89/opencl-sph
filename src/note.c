#include <stdio.h>
#include <stdarg.h>

static unsigned int _loglevel = 1;

void set_log_level(unsigned int loglevel) {
    _loglevel = loglevel;
}

void note(unsigned int level, const char * fmt, ...) {
    if (level < _loglevel) return;

    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}
