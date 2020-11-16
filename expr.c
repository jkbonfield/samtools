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
// - add floating point (use double for either?)
// - ?: operator for conditionals?
// - string literals (for variable comparison
// - variable lookup: supply a callback func to return value?

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "expr.h"

/*

All rules depend on previous rules except for simple_expr -> expression
for the primary recursion.  This avoids any backtracking and keeps the
processing very simple.

simple_expr
    : identifier
    | constant
//  | string ?
    | '(' expression ')'

unary_expr
    : simple_expr
    | '+' simple_expr
    | '-' simple_expr
    | '!' unary_expr
    | '~' unary_expr

mul_expr
    : unary_expr (
          unary_expr '*' unary_expr
        | unary_expr '/' unary_expr
        | unary_expr '%' unary_expr
      )*

add_expr
    : mul_expr (
          mul_expr '+' mul_expr
        | mul_expr '-' mul_expr
      )*

cmp_expr
    : add_expr
    | cmp_expr '<=' add_expr
    | cmp_expr '<'  add_expr
    | cmp_expr '>=' add_expr
    | cmp_expr '>'  add_expr
    | cmp_expr '='  add_expr
    | cmp_expr '!=' add_expr

and_expr
    : cmp_expr
    | and_expr '&&' cmp_expr // logical
    | and_expr '||' cmp_expr // logical
    | and_expr '&'  cmp_expr // bit-wise, eg for FLAG
    | and_expr '|'  cmp_expr // bit-wise, eg for FLAG
    | and_expr '^'  cmp_expr // bit-wise, eg for FLAG

expression
    : and_expr
 */

// Skip to start of term
static char *ws(char *str) {
    while (*str && isspace(*str))
        str++;
    return str;
}

// // Length of a term.  Used for printing up error messages.
// static char term_len(char *str) {
//     char *end = str;
//     while (*end && !isspace(*end))
//      end++;
//     return end-str;
// }

static int expression(void *data, sym_func *f,
                      char *str, char **end, int *err);

/*
 * simple_expr
 *     : identifier
 *     | constant
 * //  | string ?
 *     | '(' expression ')'
*/
static int simple_expr(void *data, sym_func *f,
                       char *str, char **end, int *err) {
    // FIXME: use double throughout so we
    // can handle int as well as floats?  Any rounding issues
    // to be concerned with?

    // Simple for now; ints only
    int ret = 0;

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

    // Otherwise a basic element
    if (isdigit(*str)) {
        ret = *str++-'0';

        while (isdigit(*str))
            ret = ret*10 + *str++-'0';

        *end = str;
    } else if (f) {
        // Look up var
        ret = f(data, str, end, err);
        if (*err) return -1;
    } else {
        *err = 1;
        return -1;
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
static int unary_expr(void *data, sym_func *f,
                      char *str, char **end, int *err) {
    int ret;
    str = ws(str);
    if (*str == '+') {
        ret = simple_expr(data, f, str+1, end, err);
    } else if (*str == '-') {
        ret = -simple_expr(data, f, str+1, end, err);
    } else if (*str == '!') {
        ret = !unary_expr(data, f, str+1, end, err);
    } else if (*str == '~') {
        ret = ~unary_expr(data, f, str+1, end, err);
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
static int mul_expr(void *data, sym_func *f,
                    char *str, char **end, int *err) {
    int ret = unary_expr(data, f, str, end, err);
    if (*err) return -1;

    str = *end;
    while (*str && !*err) {
        str = ws(str);
        if (*str == '*')
            ret = ret * unary_expr(data, f, str+1, end, err);
        else if (*str == '/')
            ret = ret / unary_expr(data, f, str+1, end, err);
        else if (*str == '%')
            ret = ret % unary_expr(data, f, str+1, end, err);
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
static int add_expr(void *data, sym_func *f,
                    char *str, char **end, int *err) {
    int ret;

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
 *     | cmp_expr '='  add_expr
 *     | cmp_expr '!=' add_expr
 */
static int cmp_expr(void *data, sym_func *f,
                    char *str, char **end, int *err) {
    int ret = add_expr(data, f, str, end, err);
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
    else if (*str == '=')
        ret = ret == cmp_expr(data, f, str+1, end, err);
    else if (strncmp(str, "!=", 2) == 0)
        ret = ret != cmp_expr(data, f, str+2, end, err);

    if (*err) ret = -1;
    return ret;
}

/*
 * and_expr
 *     : cmp_expr
 *     | and_expr 'and' cmp_expr
 *     | and_expr 'or'  cmp_expr
 */
static int and_expr(void *data, sym_func *f,
                    char *str, char **end, int *err) {
    int ret = cmp_expr(data, f, str, end, err);
    if (*err) return -1;
    str = ws(*end);
    if (strncmp(str, "&&", 2) == 0)
        ret = and_expr(data, f, str+2, end, err) && ret;
    else if (strncmp(str, "||", 2) == 0)
        ret = and_expr(data, f, str+2, end, err) || ret;
    else if (*str == '&')
        ret = and_expr(data, f, str+1, end, err) & ret;
    else if (*str == '|')
        ret = and_expr(data, f, str+1, end, err) | ret;
    else if (*str == '^')
        ret = and_expr(data, f, str+1, end, err) ^ ret;

    if (*err) ret = -1;
    return ret;
}

static int expression(void *data, sym_func *f,
                      char *str, char **end, int *err) {
    return and_expr(data, f, str, end, err);
}

int evaluate_filter(void *data, sym_func *f, char *str, int *err) {
    char *end = NULL;

    int ret = expression(data, f, str, &end, err);
    if (end && *ws(end)) {
        fprintf(stderr, "Unable to parse expression at %s\n", str);
        *err = 1;
        return -1;
    }

    return ret;
}

#if defined(EXPR_TEST) || defined(EXPR_MAIN)
typedef struct {
    int   val;
    char *str;
} test_ev;

void test(void) {
    // These are all valid expressions that should work
    test_ev tests[] = {
        {  1, "1"},
        {  1, "+1"},
        { -1, "-1"},
        {  0, "!7"},
        {  1, "!0"},
        {  1, "!(!7)"},
        {  1, "!!7"},

        {  5, "2+3"},
        { -1, "2+-3"},
        {  6, "1+2+3"},
        {  1, "-2+3"},

        {  6, "2*3"},
        {  6, "1*2*3"},
        {  0, "2*0"},

        {  7, "(7)"},
        {  7, "((7))"},
        { 21, "(1+2)*(3+4)"},
        { 14, "(4*5)-(-2*-3)"},

        {  0, "1>2"},
        {  1, "1<2"},
        {  0, "3<3"},
        {  0, "3>3"},
        {  1, "9<=9"},
        {  1, "9>=9"},
        {  1, "2*4=8"},
        {  0, "2*4!=8"},
        {  1, "4+2<3+4"},
        {  0, "4*2<3+4"},
        {  8, "4*(2<3)+4"},  // boolean; 4*(1)+4

        {  1, "2 && 1"},
        {  0, "2 && 0"},
        {  0, "0 && 2"},
        {  1, "2 || 1"},
        {  1, "2 || 0"},
        {  1, "0 || 2"},
        {  1, "1 || 2 && 3"},
        {  1, "2 && 3 || 1"},
        {  1, "0 && 3 || 2"},
        {  0, "0 && 3 || 0"},

        {  1, "3 & 1"},
        {  1, "3 & 2"},
        {  3, "1 | 2"},
        {  3, "1 | 3"},
        {  7, "1 | 6"},
        {  2, "1 ^ 3"},

        {  0, " (2*3)&7  > 4"},
        {  1, "((2*3)&7) > 4"},
        {  1, "((2*3)&7) > 4 && 2*2 <= 4"},
    };

    int err = 0, r, i;
    for (i = 0; i < sizeof(tests) / sizeof(*tests); i++) {
        if ((r=evaluate_filter(NULL, NULL, tests[i].str, &err)) != tests[i].val) {
            fprintf(stderr, "Failed test: %s == %d, got %d\n",
                    tests[i].str, tests[i].val, r);
            exit(1);
        }
    }
}
#endif

#ifdef EXPR_MAIN
int lookup(void *data, char *str, char **end, int *err) {
    int foo = 17;
    int a = 1;
    int b = 2;
    int c = 3;
    if (strncmp(str, "foo", 3) == 0) {
        *end = str+3;
        return foo;
    } else if (*str == 'a') {
        *end = str+1;
        return a;
    } else if (*str == 'b') {
        *end = str+1;
        return b;
    } else if (*str == 'c') {
        *end = str+1;
        return c;
    }
    *err = 1;
    return -1;
}

int main(int argc, char **argv) {
    if (argc > 1) {
        int err = 0;
        printf("expr = %d\n", evaluate_filter(NULL, lookup, argv[1], &err));
        printf("err = %d\n", err);
        return 0;
    }
    test();
}
#endif
