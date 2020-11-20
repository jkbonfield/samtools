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

// Could also cache fexpr_t stack here for kstring reuse?
#define MAX_REGEX 10
typedef struct sam_filter {
    char *str;
    int parsed;
    int curr_regex, max_regex;
    regex_t preg[MAX_REGEX];
} sam_filter_t;

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

static int expression(sam_filter_t *filt, void *data, sym_func *fn,
                      char *str, char **end, fexpr_t *res);

/*
 * simple_expr
 *     : identifier
 *     | constant
 * //  | string ?
 *     | '(' expression ')'
*/
static int simple_expr(sam_filter_t *filt, void *data, sym_func *fn,
                       char *str, char **end, fexpr_t *res) {
    // Main recursion step
    str = ws(str);
    if (*str == '(') {
        if (expression(filt, data, fn, str+1, end, res)) return -1;
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
        res->is_str = 0;
        res->d = d;
    } else {
        // Not valid floating point syntax.
        // FIXME: add function call names in here; len(), sqrt(), pow(), etc
        if (*str == '"') {
            // string.  FIXME: cope with backslashes at some point.
	    // Easiest with detection for \ here (rare event) and
	    // post ks_escape function to remove in-situ if needed.
            res->is_str = 1;
            char *e = str+1;
            while (*e && *e != '"')
                e++;

            kputsn(str+1, e-(str+1), ks_clear(&res->s));
            *end = e + (*e == '"');
        } else if (fn)
            // Look up variable.
            return fn(data, str, end, res);
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
static int unary_expr(sam_filter_t *filt, void *data, sym_func *fn,
                      char *str, char **end, fexpr_t *res) {
    int err;
    str = ws(str);
    if (*str == '+') {
        err = simple_expr(filt, data, fn, str+1, end, res);
        err |= res->is_str;
        res->is_true = res->d != 0;
    } else if (*str == '-') {
        err = simple_expr(filt, data, fn, str+1, end, res);
        err |= res->is_str;
        res->d = -res->d;
        res->is_true = res->d != 0;
    } else if (*str == '!') {
        err = unary_expr(filt, data, fn, str+1, end, res);
        if (res->is_str) {
            res->is_str = 0;
            res->d = 0;
        } else {
            res->d = !(int64_t)res->d;
        }
        res->is_true = !res->is_true;
    } else if (*str == '~') {
        err = unary_expr(filt, data, fn, str+1, end, res);
        err |= res->is_str;
        res->d = ~(int64_t)res->d;
        res->is_true = !res->is_true;
    } else {
        err = simple_expr(filt, data, fn, str, end, res);
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
static int mul_expr(sam_filter_t *filt, void *data, sym_func *fn,
                    char *str, char **end, fexpr_t *res) {
    if (unary_expr(filt, data, fn, str, end, res))
        return -1;

    str = *end;
    fexpr_t val = FEXPR_INIT;
    while (*str) {
        str = ws(str);
        if (*str == '*' || *str == '/' || *str == '%') {
            if (unary_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (val.is_str || res->is_str) {
                fexpr_free(&val);
                return -1; // arith on strings
            }
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
    fexpr_free(&val);

    return 0;
}

/*
 * add_expr
 *     : mul_expr (
 *           mul_expr '+' mul_expr
 *         | mul_expr '-' mul_expr
 *       )*
 */
static int add_expr(sam_filter_t *filt, void *data, sym_func *fn,
                    char *str, char **end, fexpr_t *res) {
    if (mul_expr(filt, data, fn, str, end, res))
        return -1;

    str = *end;
    fexpr_t val = FEXPR_INIT;
    while (*str) {
        str = ws(str);
        if (*str == '+' || *str == '-') {
            if (mul_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (val.is_str || res->is_str) {
                fexpr_free(&val);
                return -1; // arith on strings
            }
        }

        if (*str == '+')
            res->d += val.d;
        else if (*str == '-')
            res->d -= val.d;
        else
            break;

        str = *end;
    }
    fexpr_free(&val);

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
static int cmp_expr(sam_filter_t *filt, void *data, sym_func *fn,
                    char *str, char **end, fexpr_t *res) {
    if (add_expr(filt, data, fn, str, end, res)) return -1;

    str = ws(*end);
    fexpr_t val = FEXPR_INIT;
    int err = 0;

    if (strncmp(str, ">=", 2) == 0) {
        err = cmp_expr(filt, data, fn, str+2, end, &val);
        res->is_true=res->d = res->is_str && res->s.s && val.is_str && val.s.s
            ? strcmp(res->s.s, val.s.s) >= 0
            : !res->is_str && !val.is_str && res->d >= val.d;
        res->is_str = 0;
    } else if (*str == '>') {
        err = cmp_expr(filt, data, fn, str+1, end, &val);
        res->is_true=res->d = res->is_str && res->s.s && val.is_str && val.s.s
            ? strcmp(res->s.s, val.s.s) > 0
            : !res->is_str && !val.is_str && res->d > val.d;
        res->is_str = 0;
    } else if (strncmp(str, "<=", 2) == 0) {
        err = cmp_expr(filt, data, fn, str+2, end, &val);
        res->is_true=res->d = res->is_str && res->s.s && val.is_str && val.s.s
            ? strcmp(res->s.s, val.s.s) <= 0
            : !res->is_str && !val.is_str && res->d <= val.d;
        res->is_str = 0;
    } else if (*str == '<') {
        err = cmp_expr(filt, data, fn, str+1, end, &val);
        res->is_true=res->d = res->is_str && res->s.s && val.is_str && val.s.s
            ? strcmp(res->s.s, val.s.s) < 0
            : !res->is_str && !val.is_str && res->d < val.d;
        res->is_str = 0;
    }
    fexpr_free(&val);

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
static int eq_expr(sam_filter_t *filt, void *data, sym_func *fn,
                   char *str, char **end, fexpr_t *res) {
    if (cmp_expr(filt, data, fn, str, end, res)) return -1;

    str = ws(*end);

    int err = 0;
    fexpr_t val = FEXPR_INIT;

    // numeric vs numeric comparison is as expected
    // string vs string comparison is as expected
    // numeric vs string is false
    if (strncmp(str, "==", 2) == 0) {
        if ((err = eq_expr(filt, data, fn, str+2, end, &val))) {
            res->is_true = res->d = 0;
        } else {
            res->is_true = res->d = res->is_str
                ? (res->s.s && val.s.s ? strcmp(res->s.s, val.s.s)==0 : 0)
                : !res->is_str && !val.is_str && res->d == val.d;
        }
        res->is_str = 0;

    } else if (strncmp(str, "!=", 2) == 0) {
        if ((err = eq_expr(filt, data, fn, str+2, end, &val))) {
            res->is_true = res->d = 0;
        } else {
            res->is_true = res->d = res->is_str
                ? (res->s.s && val.s.s ? strcmp(res->s.s, val.s.s) != 0 : 1)
                : res->is_str != val.is_str || res->d != val.d;
        }
        res->is_str = 0;

    } else if (strncmp(str, "=~", 2) == 0 || strncmp(str, "!~", 2) == 0) {
        err = eq_expr(filt, data, fn, str+2, end, &val);
        if (!val.is_str || !res->is_str) {
            fexpr_free(&val);
            return -1;
        }
        if (val.s.s && res->s.s && val.is_true >= 0 && res->is_true >= 0) {
            regex_t preg_, *preg;
	    if (filt->curr_regex >= filt->max_regex) {
		// Compile regex if not seen before
		if (filt->curr_regex >= MAX_REGEX) {
		    preg = &preg_;
		} else {
		    preg = &filt->preg[filt->curr_regex];
		    filt->max_regex++;
		}

		int ec = regcomp(preg, val.s.s, REG_EXTENDED | REG_NOSUB);
		if (ec != 0) {
		    char errbuf[1024];
		    regerror(ec, preg, errbuf, 1024);
		    fprintf(stderr, "Failed regex: %.1024s\n", errbuf);
		    fexpr_free(&val);
		    return -1;
		}
            } else {
		preg = &filt->preg[filt->curr_regex];
	    }
            res->is_true = res->d = regexec(preg, res->s.s, 0, NULL, 0) == 0
                ? *str == '='  // matcn
                : *str == '!'; // no-match
	    if (preg == &preg_)
		regfree(preg);

	    filt->curr_regex++;
        } else {
            // nul regexp or input is considered false
            res->is_true = 0;
        }
        res->is_str = 0;
    }
    fexpr_free(&val);

    return err ? -1 : 0;
}

/*
 * bitand_expr
 *     : eq_expr
 *     | bitand_expr '&' eq_expr
 */
static int bitand_expr(sam_filter_t *filt, void *data, sym_func *fn,
                       char *str, char **end, fexpr_t *res) {
    if (eq_expr(filt, data, fn, str, end, res)) return -1;

    fexpr_t val = FEXPR_INIT;
    for (;;) {
        str = ws(*end);
        if (*str == '&' && str[1] != '&') {
            if (eq_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (res->is_str || val.is_str) {
                fexpr_free(&val);
                return -1;
            }
            res->is_true = res->d = (int64_t)res->d & (int64_t)val.d;
        } else {
            break;
        }
    }
    fexpr_free(&val);

    return 0;
}

/*
 * bitxor_expr
 *     : bitand_expr
 *     | bitxor_expr '^' bitand_expr
 */
static int bitxor_expr(sam_filter_t *filt, void *data, sym_func *fn,
                       char *str, char **end, fexpr_t *res) {
    if (bitand_expr(filt, data, fn, str, end, res)) return -1;

    fexpr_t val = FEXPR_INIT;
    for (;;) {
        str = ws(*end);
        if (*str == '^') {
            if (bitand_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (res->is_str || val.is_str) {
                fexpr_free(&val);
                return -1;
            }
            res->is_true = res->d = (int64_t)res->d ^ (int64_t)val.d;
        } else {
            break;
        }
    }
    fexpr_free(&val);

    return 0;
}

/*
 * bitor_expr
 *     : xor_expr
 *     | bitor_expr '|' xor_expr
 */
static int bitor_expr(sam_filter_t *filt, void *data, sym_func *fn,
                      char *str, char **end, fexpr_t *res) {
    if (bitxor_expr(filt, data, fn, str, end, res)) return -1;

    fexpr_t val = FEXPR_INIT;
    for (;;) {
        str = ws(*end);
        if (*str == '|' && str[1] != '|') {
            if (bitxor_expr(filt, data, fn, str+1, end, &val)) return -1;
            if (res->is_str || val.is_str) {
                fexpr_free(&val);
                return -1;
            }
            res->is_true = res->d = (int64_t)res->d | (int64_t)val.d;
        } else {
            break;
        }
    }
    fexpr_free(&val);

    return 0;
}

/*
 * and_expr
 *     : bitop_expr
 *     | and_expr 'and' bitop_expr
 *     | and_expr 'or'  bitop_expr
 */
static int and_expr(sam_filter_t *filt, void *data, sym_func *fn,
                    char *str, char **end, fexpr_t *res) {
    if (bitor_expr(filt, data, fn, str, end, res)) return -1;

    fexpr_t val = FEXPR_INIT;
    for (;;) {
        str = ws(*end);
        if (strncmp(str, "&&", 2) == 0) {
            if (bitor_expr(filt, data, fn, str+2, end, &val)) return -1;
            res->is_true = res->d =
                (res->is_true || (res->is_str && res->s.s) || res->d) &&
                (val.is_true  || (val.is_str && val.s.s) || val.d);
            res->is_str = 0;
        } else if (strncmp(str, "||", 2) == 0) {
            if (bitor_expr(filt, data, fn, str+2, end, &val)) return -1;
            res->is_true = res->d =
                res->is_true || (res->is_str && res->s.s) || res->d ||
                val.is_true  || (val.is_str  && val.s.s ) || val.d;
            res->is_str = 0;
        } else {
            break;
        }
    }
    fexpr_free(&val);

    return 0;
}

static int expression(sam_filter_t *filt, void *data, sym_func *fn,
                      char *str, char **end, fexpr_t *res) {
    return and_expr(filt, data, fn, str, end, res);
}

sam_filter_t *sam_filter_init(char *str) {
    sam_filter_t *f = calloc(1, sizeof(*f));
    if (!f) return NULL;

    f->str = str;
    return f;
}

void sam_filter_free(sam_filter_t *filt) {
    if (!filt)
	return;

    int i;
    for (i = 0; i < filt->max_regex; i++)
	regfree(&filt->preg[i]);

    free(filt);
}

int sam_filter_eval(sam_filter_t *filt, void *data, sym_func *fn,
		    fexpr_t *res) {
    char *end = NULL;

    memset(res, 0, sizeof(*res));

    filt->curr_regex = 0;
    if (expression(filt, data, fn, filt->str, &end, res))
        return -1;

    if (end && *ws(end)) {
        fprintf(stderr, "Unable to parse expression at %s\n", filt->str);
        return -1;
    }

    // Strings evaluate to true.  An empty string is also true, but an
    // absent (null) string is false.  An empty string has kstring length
    // of zero, but a pointer as it's nul-terminated.
    if (res->is_str)
        res->is_true = res->d = res->s.s != NULL;
    else
        res->is_true |= res->d != 0;

    return 0;
}
