/*  expr.c -- filter expression parsing and processing.

    Copyright (C) 2020 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef EXPR_H
#define EXPR_H

#include <htslib/kstring.h>

typedef struct {
    int is_str;   // Use .s vs .d
    kstring_t s;  // is_str and empty s permitted (eval as false)
    double d;     // otherwise this
} fexpr_t;

#define FEXPR_INIT {0, KS_INITIALIZE, 0}

typedef int (sym_func)(void *data, char *str, char **end, fexpr_t *res);
int evaluate_filter(void *data, sym_func *f, char *str, fexpr_t *res);

static inline void fexpr_free(fexpr_t *f) {
    ks_free(&f->s);
}

#endif /* EXPR_H */
