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

// TODO:
// - add maths functions.  pow, sqrt, log, min, max, ?
// - ?: operator for conditionals?
// - string literals (for variable comparison
// - variable lookup: supply a callback func to return value?

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>

#include "expr.h"

/*
 * This is designed to be C like with the same precedence rules.  It's not
 * full C (eg no bit-shifting), but good enough for our purposes.
 *
 * For now, all maths is strictly integer.  TODO: switch to double
 *
 * Supported syntax:
 *
 * Unary ops:     +, -, !, ~  eg -10 +10, !10 (0), ~5 (bitwise not)
 * Math ops:      +, -, *, /, %
 * Conditionals:  >, >=, <, <=, =, == (synonym for =)
 * Boolean:       &&, ||
 * Bit-wise:      &, |, ^ (XOR)
 * Grouping:      (, ),   eg "(1+2)*3"
 * Numerics:      integer only at present (no floats)
 * Variables:     any non-numeric
 */

// Skip to start of term
static char *ws(char *str) {
    while (*str && isspace(*str))
        str++;
    return str;
}

static double expression(void *data, sym_func *f,
			 char *str, char **end, int *err);

/*
 * simple_expr
 *     : identifier
 *     | constant
 * //  | string ?
 *     | '(' expression ')'
*/
static double simple_expr(void *data, sym_func *f,
			  char *str, char **end, int *err) {
    // FIXME: use double throughout so we
    // can handle int as well as floats?  Any rounding issues
    // to be concerned with?

    // Simple for now; ints only
    double ret = 0;

    // Main recursion step
    str = ws(str);
    if (*str == '(') {
        ret = expression(data, f, str+1, end, err);
        if (*err) return -1;
        str = ws(*end);
        if (*str != ')') {
            fprintf(stderr, "Missing ')'\n");
            *err = 1;
            return -1;
        }
        *end = str+1;

        return ret;
    }

    // Otherwise a basic element.
    // Ideally use hts_str2dbl, but it's internal only.
    // FIXME: update this if/when we migrate this from samtools to htslib.
    ret = strtod(str, end);
    if (str == *end) {
	// Not valid floating point syntax.
	// FIXME: add function call names in here; len(), sqrt(), pow(), etc
	if (f) {
	    // Look up variable
	    ret = f(data, str, end, err);
	    if (*err) return -1;
	} else {
	    *err = 1;
	    return -1;
	}
    }

    return ret;
}

/*
 * unary_expr
 *     : simple_expr
 *     | '+' simple_expr
 *     | '-' simple_expr
 *     | '!' unary_expr
 */
static double unary_expr(void *data, sym_func *f,
			 char *str, char **end, int *err) {
    double ret;
    str = ws(str);
    if (*str == '+') {
        ret = simple_expr(data, f, str+1, end, err);
    } else if (*str == '-') {
        ret = -simple_expr(data, f, str+1, end, err);
    } else if (*str == '!') {
        ret = !(int64_t)unary_expr(data, f, str+1, end, err);
    } else if (*str == '~') {
        ret = ~(int64_t)unary_expr(data, f, str+1, end, err);
    } else {
        ret = simple_expr(data, f, str, end, err);
    }
    return ret;
}


/*
 * mul_expr
 *     : unary_expr (
 *           unary_expr '*' unary_expr
 *         | unary_expr '/' unary_expr
 *         | unary_expr '%' unary_expr
 *       )*
 */
static double mul_expr(void *data, sym_func *f,
		       char *str, char **end, int *err) {
    double ret = unary_expr(data, f, str, end, err);
    if (*err) return -1;

    str = *end;
    while (*str && !*err) {
        str = ws(str);
        if (*str == '*')
            ret = ret * unary_expr(data, f, str+1, end, err);
        else if (*str == '/')
            ret = ret / unary_expr(data, f, str+1, end, err);
        else if (*str == '%')
            ret = (int64_t)ret % (int64_t)unary_expr(data, f, str+1, end, err);
        else
            break;

        str = *end;
    }

    return ret;
}

/*
 * add_expr
 *     : mul_expr (
 *           mul_expr '+' mul_expr
 *         | mul_expr '-' mul_expr
 *       )*
 */
static double add_expr(void *data, sym_func *f,
		       char *str, char **end, int *err) {
    double ret;

    ret = mul_expr(data, f, str, end, err);
    if (*err) return -1;

    str = *end;
    while (*str && !*err) {
        str = ws(str);
        if (*str == '+') {
            ret = ret + mul_expr(data, f, str+1, end, err);
            if (*err) return -1;
        } else if (*str == '-') {
            ret = ret - mul_expr(data, f, str+1, end, err);
            if (*err) return -1;
        } else {
            break;
        }
        str = *end;
    }

    return ret;
}

/*
 * cmp_expr
 *     : add_expr
 *     | cmp_expr '<=' add_expr
 *     | cmp_expr '<'  add_expr
 *     | cmp_expr '>=' add_expr
 *     | cmp_expr '>'  add_expr
 */
static double cmp_expr(void *data, sym_func *f,
		       char *str, char **end, int *err) {
    double ret = add_expr(data, f, str, end, err);
    if (*err) return -1;

    str = ws(*end);
    if (strncmp(str, ">=", 2) == 0)
        ret = ret >= cmp_expr(data, f, str+2, end, err);
    else if (*str == '>')
        ret = ret > cmp_expr(data, f, str+1, end, err);
    else if (strncmp(str, "<=", 2) == 0)
        ret = ret <= cmp_expr(data, f, str+2, end, err);
    else if (*str == '<')
        ret = ret < cmp_expr(data, f, str+1, end, err);

    if (*err) ret = -1;
    return ret;
}

/*
 * eq_expr
 *     : cmp_expr
 *     | eq_expr '==' cmp_expr
 *     | eq_expr '='  cmp_expr
 *     | eq_expr '!=' cmp_expr
 */
static double eq_expr(void *data, sym_func *f,
		      char *str, char **end, int *err) {
    double ret = cmp_expr(data, f, str, end, err);
    if (*err) return -1;

    str = ws(*end);
#if 1
    if (strncmp(str, "==", 2) == 0)
        ret = ret == eq_expr(data, f, str+2, end, err);
    else if (*str == '=') // synonym for ==
        ret = ret == eq_expr(data, f, str+1, end, err);
    else if (strncmp(str, "!=", 2) == 0)
        ret = ret != eq_expr(data, f, str+2, end, err);
#else
// For friendliness sake and the fact we're using doubles for representing
// integers, treat equality as within a small amount.  If we had a more
// expressive language we could do away with this perhaps.
// Eg 100.0/17*15*17/15 == 100
#define DBL_DELTA (DBL_EPSILON*100)  // approx 2e-14

    int n = 1;
    if (strncmp(str, "==", 2) == 0)
        ret -= eq_expr(data, f, str+2, end, err);
    else if (*str == '=') // synonym for ==
        ret -= eq_expr(data, f, str+1, end, err);
    else if (strncmp(str, "!=", 2) == 0)
        ret -= eq_expr(data, f, str+2, end, err), n=0;
    else
	return *err ? -1 : ret;

    ret = ret >= -DBL_DELTA && ret <= DBL_DELTA ? n : 1-n;
#endif

    if (*err) ret = -1;
    return ret;
}

/*
 * bitand_expr
 *     : eq_expr
 *     | bitand_expr '&' eq_expr
 */
static double bitand_expr(void *data, sym_func *f,
			  char *str, char **end, int *err) {
    double ret = eq_expr(data, f, str, end, err);
    if (*err) return -1;

    for (;;) {
        str = ws(*end);
        if (*str == '&' && str[1] != '&')
            ret = (int64_t)eq_expr(data, f, str+1, end, err) & (int64_t)ret;
        else
            break;
    }

    if (*err) ret = -1;
    return ret;
}

/*
 * bitxor_expr
 *     : bitand_expr
 *     | bitxor_expr '^' bitand_expr
 */
static double bitxor_expr(void *data, sym_func *f,
			  char *str, char **end, int *err) {
    double ret = bitand_expr(data, f, str, end, err);
    if (*err) return -1;

    for (;;) {
        str = ws(*end);
        if (*str == '^')
            ret = (int64_t)bitand_expr(data, f, str+1, end,err) ^ (int64_t)ret;
        else
            break;
    }

    if (*err) ret = -1;
    return ret;
}

/*
 * bitor_expr
 *     : xor_expr
 *     | bitor_expr '|' xor_expr
 */
static double bitor_expr(void *data, sym_func *f,
			 char *str, char **end, int *err) {
    double ret = bitxor_expr(data, f, str, end, err);
    if (*err) return -1;

    for (;;) {
        str = ws(*end);
        if (*str == '|' && str[1] != '|')
            ret = (int64_t)bitxor_expr(data, f, str+1, end,err) | (int64_t)ret;
        else
            break;
    }

    if (*err) ret = -1;
    return ret;
}

/*
 * and_expr
 *     : bitop_expr
 *     | and_expr 'and' bitop_expr
 *     | and_expr 'or'  bitop_expr
 */
static double and_expr(void *data, sym_func *f,
		       char *str, char **end, int *err) {
    double ret = bitor_expr(data, f, str, end, err);
    if (*err) return -1;

    for (;;) {
        str = ws(*end);
        if (strncmp(str, "&&", 2) == 0)
            ret = bitor_expr(data, f, str+2, end, err) && ret;
        else if (strncmp(str, "||", 2) == 0)
            ret = bitor_expr(data, f, str+2, end, err) || ret;
        else
            break;
    }

    if (*err) ret = -1;
    return ret;
}

static double expression(void *data, sym_func *f,
			 char *str, char **end, int *err) {
    return and_expr(data, f, str, end, err);
}

double evaluate_filter(void *data, sym_func *f, char *str, int *err) {
    char *end = NULL;

    double ret = expression(data, f, str, &end, err);
    if (end && *ws(end)) {
        fprintf(stderr, "Unable to parse expression at %s\n", str);
        *err = 1;
        return -1;
    }

    return ret;
}
