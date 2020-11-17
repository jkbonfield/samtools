/*  test-expr.c -- Testing: filter expression parsing and processing.

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

#include <config.h>

#include <stdio.h>
#include <string.h>
#include "../expr.h"

int lookup(void *data, char *str, char **end, int *err) {
    int foo = 15551; // my favourite palindromic prime
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

typedef struct {
    int   val;
    char *str;
} test_ev;

int test(void) {
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

	{  1, "(1+2)*3=9"},
	{  1, "(1+2)*3==9"},
	{  1, "(1+2)*3!=8"},
	{  0, "(1+2)*3!=9"},
	{  0, "(1+2)*3==8"},

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

	{  1, "(1<2) == (3>2)"},
	{  1, "1<2 == 3>2"},

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
        {  2, "3 & 2"},
        {  3, "1 | 2"},
        {  3, "1 | 3"},
        {  7, "1 | 6"},
        {  2, "1 ^ 3"},

	{  1, "(1^0)&(4^3)"},
	{  2, "1 ^(0&4)^ 3"},
	{  2, "1 ^ 0&4 ^ 3"},  // precedence, & before ^

	{  6, "(1|0)^(4|3)"},
	{  7, "1 |(0^4)| 3"},
	{  7, "1 | 0^4 | 3"},  // precedence, ^ before |

	{  1, "4 & 2 || 1"},
	{  1, "(4 & 2) || 1"},
	{  0, "4 & (2 || 1)"},
	{  1, "1 || 4 & 2"},
	{  1, "1 || (4 & 2)"},
	{  0, "(1 || 4) & 2"},

        {  0, " (2*3)&7  > 4"},
        {  1, "((2*3)&7) > 4"},
        {  1, "((2*3)&7) > 4 && 2*2 <= 4"},
    };

    int err = 0, r, i;
    for (i = 0; i < sizeof(tests) / sizeof(*tests); i++) {
        if ((r=evaluate_filter(NULL, lookup, tests[i].str, &err))
	    != tests[i].val) {
            fprintf(stderr, "Failed test: %s == %d, got %d\n",
                    tests[i].str, tests[i].val, r);
            return 1;
        }
    }

    return 0;
}

int main(int argc, char **argv) {
    if (argc > 1) {
        int err = 0;
        printf("expr = %d\n", evaluate_filter(NULL, lookup, argv[1], &err));
        printf("err = %d\n", err);
        return 0;
    }

    return test();
}
