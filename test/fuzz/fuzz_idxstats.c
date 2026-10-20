/*  test/fuzz/fuzz_idxstats.c -- Fuzz driver for samtools idxstats

    Copyright (C) 2026 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#define _GNU_SOURCE
#include <string.h>
#include <setjmp.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

void *memmem(const void *haystack, size_t haystacklen,
             const void *needle, size_t needlelen);

int bam_idxstats(int argc, char *argv[]);

/*
Fuzzes samtools idxstats

Build with ./build_idxstats.sh

Run with: ./run_idx_stats.sh
export ASAN_SYMBOLIZER_PATH=/software/badger/opt/llvm/21.1.8/bin/llvm-symbolizer 
ASAN_OPTIONS=allow_addr2line=1:abort_on_error=1:detect_leaks=0 ./a.out  -rss_limit_mb=2560 -timeout=10 

We construct valid argument combinations in addition to the input data.
Both are controlled by the fuzzer.

Arguments may attempt to load files exclude,file, skip.bed, and RG.txt.
These should be constant that are viable given the initial test corpus.
Modifying the corpus will permit testing of invalid data relative to the
filtering files.

No validity of result is assumed.  We're simply checking for crashes or
illegal program behaviour.
*/

// Taken from bamtk.c which we exclude from our build line
const char *samtools_version(void) {
    return "fuzzer";
}

// Returns "data:,<base64-data>" as a malloced buffer
char *base64_encode(const uint8_t *data, size_t len) {
    static const char table[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    size_t out_len = 4 * ((len + 2) / 3);

    char *output = malloc(out_len + 20);
    if (!output)
        return NULL;

    size_t i = 0, j = 0;
    // data: URI
    strcpy(output, "data:;base64,");
    j = strlen(output);

    // base-64
    while (i < len) {
        unsigned int a =             data[i++];
        unsigned int b = (i < len) ? data[i++] : 0;
        unsigned int c = (i < len) ? data[i++] : 0;

        unsigned int triple = (a << 16) | (b << 8) | c;

        output[j++] = table[(triple >> 18) & 0x3F];
        output[j++] = table[(triple >> 12) & 0x3F];
        output[j++] = (i - 1 < len) ? table[(triple >> 6) & 0x3F] : '=';
        output[j++] = (i < len)     ? table[triple & 0x3F]        : '=';
    }

    output[j] = '\0';
    return output;
}

void dump(int argc, char **argv) {
    fprintf(stderr, "Calling: ");
    for (int i = 0; i < argc; i++)
	fprintf(stderr, " '%s'", argv[i]);
    fprintf(stderr, "\n");
}

// Mpileup calls exit() directly.  This is problamtic for the fuzzer, but
// we can compile it with -Dexit=fuzz_exit and write our own fuzz_exit
// function that unwinds back to LLVMFuzzerTestOneInput via setjmp
static jmp_buf fuzz_jmp;
void fuzz_exit(int status) {
    fprintf(stderr, "Caught exit of status %d\n", status);
    longjmp(fuzz_jmp, status + 1);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    int jret = setjmp(fuzz_jmp);
    if (jret != 0) {
	// if non-zero don't add it to the test corpus
	return jret ? -1 : 0;
    }

    if (size < 10)
	return 0;

    // Our input data is <bam/sam/cram data>#idx#<index data>.
    // Reject inputs that we cannot split in two
    uint8_t *idx;
    if (!(idx = memmem(data, size, "#idx#", 5)))
        return -1;

    // Split into two
    const uint8_t *data1 = data;
    size_t size1 = idx - data;
    const uint8_t *data2 = idx+5;
    size_t size2 = data+size - data2;

    char *argv[100] = {0}, argc = 0;

    argv[argc++] = "samtools";
    argv[argc++] = "idxstats";

    // Create a data file to read in.
    char *b64 = base64_encode(data1, size1);
    char *c64 = base64_encode(data2, size2);
    argv[argc++] = "-X";
    argv[argc++] = b64;
    argv[argc++] = c64;

    dump(argc, argv);

    optind = 1;
    opterr = 0;

    int ret = bam_idxstats(argc-1, argv+1);
    free(b64);
    free(c64);

    return 0;
}
