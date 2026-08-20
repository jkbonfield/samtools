#!/bin/sh -x
#
#    Copyright (C) 2026 Genome Research Ltd.
#
#    Author: James Bonfield <jkb@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

CC=${CC:-clang21}

# Build htslib too
# (cd ../../../htslib; make clean; make CC="$CC -fsanitize=address,undefined" CFLAGS="-g -O3 -DFUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION" -j8)

# cd ../..;make -j8 CC="$CC -fsanitize=address,undefined"  CFLAGS=-g

(cd ../..;$CC -I. -I../htslib -fsanitize=address,undefined,fuzzer -g -DFUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION -Dexit=fuzz_exit test/fuzz/fuzz_view.c -L./lz4 bam_aux.o bam_index.o bam_plcmd.o sam_view.c bam_fastq.o bam_cat.o bam_md.o bam_plbuf.o bam_reheader.o bam_sort.o bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o bam2bcf.o sample.o cut_target.o phase.o bam2depth.o coverage.o padding.o bedcov.o bamshuf.o faidx.o dict.o stats.o stats_isize.o bam_flags.o bam_split.o bam_tview.o bam_tview_curses.o bam_tview_html.o bam_lpileup.o bam_quickcheck.o bam_addrprg.o bam_markdup.o tmp_file.o bam_ampliconclip.o amplicon_stats.o bam_import.o bam_samples.o bam_consensus.o consensus_pileup.o reference.o reset.o cram_size.o bam_checksum.o ./lz4/lz4.o libst.a ../htslib/libhts.a -lpthread -lz -lm -lbz2 -llzma -ldeflate -lcurl -lcrypto -lncursesw -lm -lz  -lpthread -o test/fuzz/fuzz_view)
