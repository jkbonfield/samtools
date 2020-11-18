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
#include <regex.h> // may need configure rule for this

#include "expr.h"

/*
 * This is designed to be mostly C like with mostly same the precedence rules.
 * It's not full C (eg no bit-shifting), but good enough for our purposes.
 *
 * Supported syntax, in order of precedence:
 *
 * Grouping:      (, ),   eg "(1+2)*3"
 * Values:        integers, floats, strings or variables
 * Unary ops:     +, -, !, ~  eg -10 +10, !10 (0), ~5 (bitwise not)
 * Math ops:      *, /, %  [TODO: add // for floor division?]
 * Math ops:      +, -
 * Conditionals:  >, >=, <, <=,
 * Equality:      ==, !=, =~ !~
 * Bit-wise:      &, |, ^  [NB as 3 precedence levels, in that order]
 * Boolean:       &&, ||
 */

// Skip to start of term
static char *ws(char *str) {
    while (*str && isspace(*str))
        str++;
    return str;
}

static int expression(void *data, sym_func *f,
		      char *str, char **end, fexpr_t *res);

/*
 * simple_expr
 *     : identifier
 *     | constant
 * //  | string ?
 *     | '(' expression ')'
*/
static int simple_expr(void *data, sym_func *f,
		       char *str, char **end, fexpr_t *res) {
    // Main recursion step
    str = ws(str);
    if (*str == '(') {
        if (expression(data, f, str+1, end, res)) return -1;
        str = ws(*end);
        if (*str != ')') {
            fprintf(stderr, "Missing ')'\n");
            return -1;
        }
        *end = str+1;

        return 0;
    }

    // Otherwise a basic element.
    // Ideally use hts_str2dbl, but it's internal only.
    // FIXME: update this if/when we migrate this from samtools to htslib.
    double d = strtod(str, end);
    if (str != *end) {
	res->s = NULL;
	res->d = d;
    } else {
	// Not valid floating point syntax.
	// FIXME: add function call names in here; len(), sqrt(), pow(), etc
	if (*str == '"') {
	    // string.  FIXME: cope with backslashes at some point.
	    res->s = str+1;
	    char *e = str+1;
	    while (*e && *e != '"')
		e++;

	    // FIXME: this modifies the query string.
	    // Alternatives are taking a copy and all the memory allocation
	    // bits that go with it, or storing string plus length.
	    //
	    // If we want to cope with escaping rules (\" etc) then
	    // kstring as one of the types is probably more robust.
	    if (*e == '"')
		*e++ = 0;

	    *end = e;
	} else if (f)
	    // Look up variable
	    return f(data, str, end, res);
	else
	    return -1;
    }

    return 0;
}

/*
 * unary_expr
 *     : simple_expr
 *     | '+' simple_expr
 *     | '-' simple_expr
 *     | '!' unary_expr // higher precedence
 *     | '~' unary_expr // higher precedence
 */
static int unary_expr(void *data, sym_func *f,
		      char *str, char **end, fexpr_t *res) {
    int err;
    str = ws(str);
    if (*str == '+') {
        err = simple_expr(data, f, str+1, end, res);
	err |= res->s != NULL;
    } else if (*str == '-') {
        err = simple_expr(data, f, str+1, end, res);
	err |= res->s != NULL;
	res->d = -res->d;
    } else if (*str == '!') {
        err = unary_expr(data, f, str+1, end, res);
	err |= res->s != NULL;
	res->d = !(int64_t)res->d;
    } else if (*str == '~') {
	err = unary_expr(data, f, str+1, end, res);
	err |= res->s != NULL;
        res->d = ~(int64_t)res->d;
    } else {
        err = simple_expr(data, f, str, end, res);
    }
    return err ? -1 : 0;
}


/*
 * mul_expr
 *     : unary_expr (
 *           unary_expr '*' unary_expr
 *         | unary_expr '/' unary_expr
 *         | unary_expr '%' unary_expr
 *       )*
 */
static int mul_expr(void *data, sym_func *f,
		    char *str, char **end, fexpr_t *res) {
    if (unary_expr(data, f, str, end, res))
	return -1;

    str = *end;
    while (*str) {
	fexpr_t val = {NULL, 0};
        str = ws(str);
	if (*str == '*' || *str == '/' || *str == '%') {
            if (unary_expr(data, f, str+1, end, &val)) return -1;
	    if (val.s || res->s) return -1; // arith on strings
	}

        if (*str == '*')
	    res->d *= val.d;
	else if (*str == '/')
	    res->d /= val.d;
	else if (*str == '%')
	    res->d = (int64_t)res->d % (int64_t)val.d;
	else
            break;

        str = *end;
    }

    return 0;
}

/*
 * add_expr
 *     : mul_expr (
 *           mul_expr '+' mul_expr
 *         | mul_expr '-' mul_expr
 *       )*
 */
static int add_expr(void *data, sym_func *f,
		    char *str, char **end, fexpr_t *res) {
    if (mul_expr(data, f, str, end, res))
	return -1;

    str = *end;
    while (*str) {
	fexpr_t val = {NULL, 0};
        str = ws(str);
	if (*str == '+' || *str == '-') {
            if (mul_expr(data, f, str+1, end, &val)) return -1;
	    if (val.s || res->s) return -1; // arith on strings
	}

        if (*str == '+')
	    res->d += val.d;
	else if (*str == '-')
	    res->d -= val.d;
	else
            break;

        str = *end;
    }

    return 0;
}

/*
 * cmp_expr
 *     : add_expr
 *     | cmp_expr '<=' add_expr
 *     | cmp_expr '<'  add_expr
 *     | cmp_expr '>=' add_expr
 *     | cmp_expr '>'  add_expr
 */
static int cmp_expr(void *data, sym_func *f,
		    char *str, char **end, fexpr_t *res) {
    if (add_expr(data, f, str, end, res)) return -1;

    str = ws(*end);
    fexpr_t val = {NULL, 0};
    int err = 0, numeric = 1;

    // Maybe consider > and < on strings.
    // eg "abba" < "acdc"?
    if (strncmp(str, ">=", 2) == 0) {
	err = cmp_expr(data, f, str+2, end, &val);
	res->d = res->d >= val.d;
    } else if (*str == '>') {
	err = cmp_expr(data, f, str+1, end, &val);
	res->d = res->d > val.d;
    } else if (strncmp(str, "<=", 2) == 0) {
	err = cmp_expr(data, f, str+2, end, &val);
	res->d = res->d <= val.d;
    } else if (*str == '<') {
	err = cmp_expr(data, f, str+1, end, &val);
	res->d = res->d < val.d;
    } else {
	numeric = 0;
    }

    if (numeric && (val.s || res->s))
	return -1;

    return err ? -1 : 0;
}

/*
 * eq_expr
 *     : cmp_expr
 *     | eq_expr '==' cmp_expr
 *     | eq_expr '!=' cmp_expr
 *     | eq_expr '=~' cmp_expr
 *     | eq_expr '!~' cmp_expr
 */
static int eq_expr(void *data, sym_func *f,
		   char *str, char **end, fexpr_t *res) {
    if (cmp_expr(data, f, str, end, res)) return -1;

    str = ws(*end);
#if 1
    int err = 0;
    fexpr_t val = {NULL, 0};
    if (strncmp(str, "==", 2) == 0) {
	// TODO: add =~ for strings
	err = eq_expr(data, f, str+2, end, &val);
	res->d = res->s && val.s
	    ? strcmp(res->s, val.s)==0
	    : (res->s || val.s) ? 0 : res->d == val.d;
	res->s = NULL;
    } else if (strncmp(str, "!=", 2) == 0) {
	err = eq_expr(data, f, str+2, end, &val);
	res->d = res->s && val.s
	    ? strcmp(res->s, val.s)!=0
	    : (res->s || val.s) ? 0 : res->d != val.d;
	res->s = NULL;
    } else if (strncmp(str, "=~", 2) == 0 || strncmp(str, "!~", 2) == 0) {
	err = eq_expr(data, f, str+2, end, &val);
	if (!val.s || !res->s) return -1;
	regex_t preg;
	int ec;
	// FIXME: cache compiled regexp
	if ((ec = regcomp(&preg, val.s, REG_EXTENDED | REG_NOSUB)) != 0) {
	    char errbuf[1024];
	    regerror(ec, &preg, errbuf, 1024);
	    fprintf(stderr, "Failed regex: %.1024s\n", errbuf);
	    return -1;
	}
	res->d = regexec(&preg, res->s, 0, NULL, 0) == 0
	    ? *str == '='  // matcn
	    : *str == '!'; // no-match
	res->s = NULL;
	regfree(&preg);
    }

#else
// For friendliness sake and the fact we're using fexpr_ts for representing
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

    return err ? -1 : 0;
}

/*
 * bitand_expr
 *     : eq_expr
 *     | bitand_expr '&' eq_expr
 */
static int bitand_expr(void *data, sym_func *f,
		       char *str, char **end, fexpr_t *res) {
    if (eq_expr(data, f, str, end, res)) return -1;

    for (;;) {
        str = ws(*end);
        if (*str == '&' && str[1] != '&') {
	    fexpr_t val = {NULL, 0};
	    if (eq_expr(data, f, str+1, end, &val)) return -1;
	    if (res->s || val.s) return -1;
	    res->d = (int64_t)res->d & (int64_t)val.d;
	} else {
            break;
	}
    }

    return 0;
}

/*
 * bitxor_expr
 *     : bitand_expr
 *     | bitxor_expr '^' bitand_expr
 */
static int bitxor_expr(void *data, sym_func *f,
		       char *str, char **end, fexpr_t *res) {
    if (bitand_expr(data, f, str, end, res)) return -1;

    for (;;) {
        str = ws(*end);
        if (*str == '^') {
	    fexpr_t val = {NULL, 0};
	    if (bitand_expr(data, f, str+1, end, &val)) return -1;
	    if (res->s || val.s) return -1;
	    res->d = (int64_t)res->d ^ (int64_t)val.d;
	} else {
            break;
	}
    }

    return 0;
}

/*
 * bitor_expr
 *     : xor_expr
 *     | bitor_expr '|' xor_expr
 */
static int bitor_expr(void *data, sym_func *f,
		      char *str, char **end, fexpr_t *res) {
    if (bitxor_expr(data, f, str, end, res)) return -1;

    for (;;) {
        str = ws(*end);
        if (*str == '|' && str[1] != '|') {
	    fexpr_t val = {NULL, 0};
	    if (bitxor_expr(data, f, str+1, end, &val)) return -1;
	    if (res->s || val.s) return -1;
	    res->d = (int64_t)res->d | (int64_t)val.d;
	} else {
            break;
	}
    }

    return 0;
}

/*
 * and_expr
 *     : bitop_expr
 *     | and_expr 'and' bitop_expr
 *     | and_expr 'or'  bitop_expr
 */
static int and_expr(void *data, sym_func *f,
		    char *str, char **end, fexpr_t *res) {
    if (bitor_expr(data, f, str, end, res)) return -1;

    for (;;) {
        str = ws(*end);
	fexpr_t val = {NULL, 0};
        if (strncmp(str, "&&", 2) == 0) {
	    if (bitor_expr(data, f, str+2, end, &val)) return -1;
	    res->d = (res->s || res->d) && (val.s || val.d);
	} else if (strncmp(str, "||", 2) == 0) {
	    if (bitor_expr(data, f, str+2, end, &val)) return -1;
	    res->d = res->s || res->d || val.s || val.d;
	} else {
            break;
	}
    }

    return 0;
}

static int expression(void *data, sym_func *f,
		      char *str, char **end, fexpr_t *res) {
    return and_expr(data, f, str, end, res);
}

int evaluate_filter(void *data, sym_func *f, char *str, fexpr_t *res) {
    char *end = NULL;

    if (expression(data, f, str, &end, res))
	return -1;

    if (end && *ws(end)) {
        fprintf(stderr, "Unable to parse expression at %s\n", str);
        return -1;
    }

    return 0;
}
