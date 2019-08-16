# Makefile for samtools, utilities for the Sequence Alignment/Map format.
#
#    Copyright (C) 2008-2017 Genome Research Ltd.
#    Portions copyright (C) 2010-2012 Broad Institute.
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

CC       = gcc
AR       = ar
CPPFLAGS =
#CFLAGS   = -g -Wall -O2 -pedantic -std=c99 -D_XOPEN_SOURCE=600
CFLAGS   = -g -Wall -O2
LDFLAGS  =
LIBS     =

LZ4DIR   = ./lz4
LZ4_CPPFLAGS = -I$(LZ4DIR)
LZ4_LDFLAGS  = -L$(LZ4DIR)



LOBJS=      bam_aux.o bam.o sam.o \
            bam_plbuf.o
AOBJS=      bam_index.o bam_plcmd.o sam_view.o bam_fastq.o \
            bam_cat.o bam_md.o bam_reheader.o bam_sort.o bedidx.o \
            bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o \
            bamtk.o bam2bcf.o bam2bcf_indel.o sample.o \
            cut_target.o phase.o bam2depth.o coverage.o padding.o bedcov.o bamshuf.o \
            faidx.o dict.o stats.o stats_isize.o bam_flags.o bam_split.o \
            bam_tview.o bam_tview_curses.o bam_tview_html.o bam_lpileup.o \
            bam_quickcheck.o bam_addrprg.o bam_markdup.o tmp_file.o
LZ4OBJS  =  $(LZ4DIR)/lz4.o

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
datarootdir = $(prefix)/share
mandir      = $(datarootdir)/man
man1dir     = $(mandir)/man1

# Installation location for $(MISC_PROGRAMS) and $(MISC_SCRIPTS)
misc_bindir = $(bindir)

MKDIR_P = mkdir -p
INSTALL = install -p
INSTALL_DATA    = $(INSTALL) -m 644
INSTALL_DIR     = $(MKDIR_P) -m 755
INSTALL_MAN     = $(INSTALL_DATA)
INSTALL_PROGRAM = $(INSTALL)
INSTALL_SCRIPT  = $(INSTALL_PROGRAM)


PROGRAMS = samtools

MISC_PROGRAMS = \
	misc/ace2sam misc/maq2sam-long misc/maq2sam-short \
	misc/md5fa misc/md5sum-lite misc/wgsim

MISC_SCRIPTS = \
	misc/blast2sam.pl misc/bowtie2sam.pl misc/export2sam.pl \
	misc/interpolate_sam.pl misc/novo2sam.pl \
	misc/plot-bamstats misc/psl2sam.pl \
	misc/sam2vcf.pl misc/samtools.pl misc/seq_cache_populate.pl \
	misc/soap2sam.pl \
	misc/varfilter.py misc/wgsim_eval.pl misc/zoom2sam.pl

TEST_PROGRAMS = \
	test/merge/test_bam_translate \
	test/merge/test_rtrans_build \
	test/merge/test_trans_tbl_init \
	test/split/test_count_rg \
	test/split/test_expand_format_string \
	test/split/test_filter_header_rg \
	test/split/test_parse_args \
	test/vcf-miniview

all: $(PROGRAMS) $(MISC_PROGRAMS) $(TEST_PROGRAMS)

ALL_CPPFLAGS = -I. $(HTSLIB_CPPFLAGS) $(LZ4_CPPFLAGS) $(CPPFLAGS)
ALL_LDFLAGS  = $(HTSLIB_LDFLAGS) $(LZ4_LDFLAGS) $(LDFLAGS)
ALL_LIBS     = -lz $(LIBS)

# Usually config.mk and config.h are generated by running configure
# or config.status, but if those aren't used create defaults here.

config.mk:
	@sed -e '/^prefix/,/^LIBS/d;s/@Hsource@//;s/@Hinstall@/#/;s#@HTSDIR@#../htslib#g;s/@HTSLIB_CPPFLAGS@/-I$$(HTSDIR)/g;s/@CURSES_LIB@/-lcurses/g' config.mk.in > $@

config.h:
	echo '/* Basic config.h generated by Makefile */' > $@
	echo '#define HAVE_CURSES' >> $@
	echo '#define HAVE_CURSES_H' >> $@

include config.mk

# If not using GNU make, you need to copy the version number from version.sh
# into here.
PACKAGE_VERSION = $(shell ./version.sh)

# Force version.h to be remade if $(PACKAGE_VERSION) has changed.
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))

# If you don't have GNU Make but are building from a Git repository, you may
# wish to replace this with a rule that always rebuilds version.h:
# version.h: force
#	echo '#define SAMTOOLS_VERSION "`git describe --always --dirty`"' > $@
version.h:
	echo '#define SAMTOOLS_VERSION "$(PACKAGE_VERSION)"' > $@

print-version:
	@echo $(PACKAGE_VERSION)


.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) $(ALL_CPPFLAGS) -c -o $@ $<

LIBST_OBJS = sam_opts.o sam_utils.o


lib:libbam.a

libbam.a:$(LOBJS)
	$(AR) -csru $@ $(LOBJS)

samtools: $(AOBJS) $(LZ4OBJS) libbam.a libst.a $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ $(AOBJS) $(LZ4OBJS) libbam.a libst.a $(HTSLIB_LIB) $(CURSES_LIB) -lm $(ALL_LIBS) -lpthread

# For building samtools and its test suite only: NOT to be installed.
libst.a: $(LIBST_OBJS)
	@-rm -f $@
	$(AR) -rcs $@ $(LIBST_OBJS)


bam_h = bam.h $(htslib_bgzf_h) $(htslib_sam_h)
bam2bcf_h = bam2bcf.h $(htslib_hts_h) $(htslib_vcf_h)
bam_lpileup_h = bam_lpileup.h $(htslib_sam_h)
bam_plbuf_h = bam_plbuf.h $(htslib_sam_h)
bam_tview_h = bam_tview.h $(htslib_hts_h) $(htslib_sam_h) $(htslib_faidx_h) $(bam2bcf_h) $(htslib_khash_h) $(bam_lpileup_h)
bedidx_h = bedidx.h $(htslib_hts_h)
sam_h = sam.h $(htslib_sam_h) $(bam_h)
sam_opts_h = sam_opts.h $(htslib_hts_h)
sample_h = sample.h $(htslib_kstring_h)
samtools_h = samtools.h $(htslib_hts_defs_h) $(htslib_sam_h)
stats_isize_h = stats_isize.h $(htslib_khash_h)
tmp_file_h = tmp_file.h $(htslib_sam_h) $(LZ4DIR)/lz4.h

bam.o: bam.c config.h $(bam_h) $(htslib_kstring_h)
bam2bcf.o: bam2bcf.c config.h $(htslib_hts_h) $(htslib_sam_h) $(htslib_kstring_h) $(htslib_kfunc_h) $(bam2bcf_h)
bam2bcf_indel.o: bam2bcf_indel.c config.h $(htslib_hts_h) $(htslib_sam_h) $(bam2bcf_h) $(htslib_khash_h) $(htslib_ksort_h)
bam2depth.o: bam2depth.c config.h $(htslib_sam_h) $(samtools_h) $(sam_opts_h)
coverage.o: coverage.c config.h $(htslib_sam_h) $(htslib_hts_h) $(samtools_h) $(sam_opts_h)
bam_addrprg.o: bam_addrprg.c config.h $(htslib_sam_h) $(htslib_kstring_h) $(samtools_h) $(htslib_thread_pool_h) $(sam_opts_h)
bam_aux.o: bam_aux.c config.h $(bam_h)
bam_cat.o: bam_cat.c config.h $(htslib_bgzf_h) $(htslib_sam_h) $(htslib_cram_h) $(htslib_khash_h) $(samtools_h)
bam_color.o: bam_color.c config.h $(bam_h)
bam_fastq.o: bam_fastq.c config.h $(htslib_sam_h) $(htslib_klist_h) $(htslib_kstring_h) $(htslib_bgzf_h) $(htslib_thread_pool_h) $(samtools_h) $(sam_opts_h)
bam_index.o: bam_index.c config.h $(htslib_hts_h) $(htslib_sam_h) $(htslib_khash_h) $(samtools_h) $(sam_opts_h)
bam_lpileup.o: bam_lpileup.c config.h $(bam_plbuf_h) $(bam_lpileup_h) $(htslib_ksort_h)
bam_mate.o: bam_mate.c config.h $(htslib_thread_pool_h) $(sam_opts_h) $(htslib_kstring_h) $(htslib_sam_h) $(samtools_h)
bam_md.o: bam_md.c config.h $(htslib_faidx_h) $(htslib_sam_h) $(htslib_kstring_h) $(htslib_thread_pool_h) $(sam_opts_h) $(samtools_h)
bam_plbuf.o: bam_plbuf.c config.h $(htslib_hts_h) $(htslib_sam_h) $(bam_plbuf_h)
bam_plcmd.o: bam_plcmd.c config.h $(htslib_sam_h) $(htslib_faidx_h) $(htslib_kstring_h) $(htslib_khash_str2int_h) $(samtools_h) $(sam_opts_h) $(bam2bcf_h) $(sample_h)
bam_quickcheck.o: bam_quickcheck.c config.h $(htslib_hts_h) $(htslib_sam_h)
bam_reheader.o: bam_reheader.c config.h $(htslib_bgzf_h) $(htslib_sam_h) $(htslib_hfile_h) $(htslib_cram_h) $(samtools_h)
bam_rmdup.o: bam_rmdup.c config.h $(htslib_sam_h) $(sam_opts_h) $(samtools_h) $(bam_h) $(htslib_khash_h)
bam_rmdupse.o: bam_rmdupse.c config.h $(bam_h) $(htslib_sam_h) $(htslib_khash_h) $(htslib_klist_h) $(samtools_h)
bam_sort.o: bam_sort.c config.h $(htslib_ksort_h) $(htslib_hts_os_h) $(htslib_khash_h) $(htslib_klist_h) $(htslib_kstring_h) $(htslib_sam_h) $(sam_opts_h) $(samtools_h)
bam_split.o: bam_split.c config.h $(htslib_sam_h) $(htslib_khash_h) $(htslib_kstring_h) $(htslib_cram_h) $(htslib_thread_pool_h) $(sam_opts_h) $(samtools_h)
bam_stat.o: bam_stat.c config.h $(htslib_sam_h) $(samtools_h) $(sam_opts_h)
bam_tview.o: bam_tview.c config.h $(bam_tview_h) $(htslib_faidx_h) $(htslib_sam_h) $(htslib_bgzf_h) $(samtools_h) $(sam_opts_h)
bam_tview_curses.o: bam_tview_curses.c config.h $(bam_tview_h)
bam_tview_html.o: bam_tview_html.c config.h $(bam_tview_h)
bam_flags.o: bam_flags.c config.h $(htslib_sam_h)
bamshuf.o: bamshuf.c config.h $(htslib_sam_h) $(htslib_hts_h) $(htslib_ksort_h) $(samtools_h) $(htslib_thread_pool_h) $(sam_opts_h) $(htslib_khash_h)
bamtk.o: bamtk.c config.h $(htslib_hts_h) $(samtools_h) version.h
bedcov.o: bedcov.c config.h $(htslib_kstring_h) $(htslib_sam_h) $(htslib_thread_pool_h) $(samtools_h) $(sam_opts_h) $(htslib_kseq_h)
bedidx.o: bedidx.c config.h $(bedidx_h) $(htslib_ksort_h) $(htslib_kseq_h) $(htslib_khash_h)
cut_target.o: cut_target.c config.h $(htslib_hts_h) $(htslib_sam_h) $(htslib_faidx_h) $(samtools_h) $(sam_opts_h)
dict.o: dict.c config.h $(htslib_kseq_h) $(htslib_hts_h)
faidx.o: faidx.c config.h $(htslib_faidx_h) $(htslib_hts_h) $(htslib_hfile_h) $(htslib_kstring_h) $(samtools_h)
padding.o: padding.c config.h $(htslib_kstring_h) $(htslib_sam_h) $(htslib_faidx_h) $(sam_opts_h) $(samtools_h)
phase.o: phase.c config.h $(htslib_hts_h) $(htslib_sam_h) $(htslib_kstring_h) $(sam_opts_h) $(samtools_h) $(htslib_hts_os_h) $(htslib_kseq_h) $(htslib_khash_h) $(htslib_ksort_h)
sam.o: sam.c config.h $(htslib_faidx_h) $(sam_h)
sam_opts.o: sam_opts.c config.h $(sam_opts_h)
sam_utils.o: sam_utils.c config.h $(samtools_h)
sam_view.o: sam_view.c config.h $(htslib_sam_h) $(htslib_faidx_h) $(htslib_khash_h) $(htslib_thread_pool_h) $(samtools_h) $(sam_opts_h) $(bedidx_h)
sample.o: sample.c config.h $(sample_h) $(htslib_khash_h)
stats_isize.o: stats_isize.c config.h $(stats_isize_h) $(htslib_khash_h)
stats.o: stats.c config.h $(htslib_faidx_h) $(htslib_sam_h) $(htslib_hts_h) $(htslib_hts_defs_h) $(htslib_khash_str2int_h) $(samtools_h) $(htslib_khash_h) $(htslib_kstring_h) $(stats_isize_h) $(sam_opts_h) $(bedidx_h)
bam_markdup.o: bam_markdup.c config.h $(htslib_thread_pool_h) $(htslib_sam_h) $(sam_opts_h) $(samtools_h) $(htslib_khash_h) $(htslib_klist_h) $(htslib_kstring_h) $(tmp_file_h)
tmp_file.o: tmp_file.c config.h $(tmp_file_h) $(htslib_sam_h)

# Maintainer source code checks
# - copyright boilerplate presence
# - tab and trailing space detection
maintainer-check:
	test/maintainer/check_copyright.pl .
	test/maintainer/check_spaces.pl .

# test programs

# For tests that might use it, set $REF_PATH explicitly to use only reference
# areas within the test suite (or set it to ':' to use no reference areas).
# (regression.sh sets $REF_PATH to a subdirectory itself.)
#
# If using MSYS, avoid poor shell expansion via:
#    MSYS2_ARG_CONV_EXCL="*" make check
check test: samtools $(BGZIP) $(TEST_PROGRAMS)
	test/split/test_count_rg
	test/split/test_expand_format_string
	test/split/test_filter_header_rg
	test/split/test_parse_args
	REF_PATH=: test/test.pl --exec bgzip=$(BGZIP) $${TEST_OPTS:-}
	test/merge/test_bam_translate test/merge/test_bam_translate.tmp
	test/merge/test_rtrans_build
	test/merge/test_trans_tbl_init
	cd test/mpileup && ./regression.sh mpileup.reg
	cd test/mpileup && ./regression.sh depth.reg


test/merge/test_bam_translate: test/merge/test_bam_translate.o test/test.o libst.a $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ test/merge/test_bam_translate.o test/test.o libst.a $(HTSLIB_LIB) $(ALL_LIBS) -lpthread

test/merge/test_rtrans_build: test/merge/test_rtrans_build.o test/test.o libst.a $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ test/merge/test_rtrans_build.o test/test.o libst.a $(HTSLIB_LIB) $(ALL_LIBS) -lpthread

test/merge/test_trans_tbl_init: test/merge/test_trans_tbl_init.o test/test.o libst.a $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ test/merge/test_trans_tbl_init.o test/test.o libst.a $(HTSLIB_LIB) $(ALL_LIBS) -lpthread

test/split/test_count_rg: test/split/test_count_rg.o test/test.o libst.a $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ test/split/test_count_rg.o test/test.o libst.a $(HTSLIB_LIB) $(ALL_LIBS) -lpthread

test/split/test_expand_format_string: test/split/test_expand_format_string.o test/test.o libst.a $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ test/split/test_expand_format_string.o test/test.o libst.a $(HTSLIB_LIB) $(ALL_LIBS) -lpthread

test/split/test_filter_header_rg: test/split/test_filter_header_rg.o test/test.o libst.a $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ test/split/test_filter_header_rg.o test/test.o libst.a $(HTSLIB_LIB) $(ALL_LIBS) -lpthread

test/split/test_parse_args: test/split/test_parse_args.o test/test.o libst.a $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ test/split/test_parse_args.o test/test.o libst.a $(HTSLIB_LIB) $(ALL_LIBS) -lpthread

test/vcf-miniview: test/vcf-miniview.o $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ test/vcf-miniview.o $(HTSLIB_LIB) $(ALL_LIBS) -lpthread

test_test_h = test/test.h $(htslib_sam_h)

test/merge/test_bam_translate.o: test/merge/test_bam_translate.c config.h bam_sort.o $(test_test_h)
test/merge/test_rtrans_build.o: test/merge/test_rtrans_build.c config.h bam_sort.o
test/merge/test_trans_tbl_init.o: test/merge/test_trans_tbl_init.c config.h bam_sort.o
test/split/test_count_rg.o: test/split/test_count_rg.c config.h bam_split.o $(test_test_h)
test/split/test_expand_format_string.o: test/split/test_expand_format_string.c config.h bam_split.o $(test_test_h)
test/split/test_filter_header_rg.o: test/split/test_filter_header_rg.c config.h bam_split.o $(test_test_h)
test/split/test_parse_args.o: test/split/test_parse_args.c config.h bam_split.o $(test_test_h)
test/test.o: test/test.c config.h $(htslib_sam_h) $(test_test_h)
test/vcf-miniview.o: test/vcf-miniview.c config.h $(htslib_vcf_h)


# misc programs

misc/ace2sam: misc/ace2sam.o
	$(CC) $(LDFLAGS) -o $@ misc/ace2sam.o $(ALL_LIBS)

misc/maq2sam-short: misc/maq2sam-short.o
	$(CC) $(LDFLAGS) -o $@ misc/maq2sam-short.o $(ALL_LIBS)

misc/maq2sam-long: misc/maq2sam-long.o
	$(CC) $(LDFLAGS) -o $@ misc/maq2sam-long.o $(ALL_LIBS)

misc/md5fa: misc/md5fa.o $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ misc/md5fa.o $(HTSLIB_LIB) $(ALL_LIBS)

misc/md5sum-lite: misc/md5sum-lite.o $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ misc/md5sum-lite.o $(HTSLIB_LIB) $(ALL_LIBS)

misc/wgsim: misc/wgsim.o $(HTSLIB)
	$(CC) $(ALL_LDFLAGS) -o $@ misc/wgsim.o -lm $(HTSLIB_LIB) $(ALL_LIBS)

misc/ace2sam.o: misc/ace2sam.c config.h $(htslib_kstring_h) $(htslib_kseq_h)
misc/md5fa.o: misc/md5fa.c config.h $(htslib_kseq_h) $(htslib_hts_h)
misc/md5sum-lite.o: misc/md5sum-lite.c config.h $(htslib_hts_h)
misc/wgsim.o: misc/wgsim.c config.h version.h $(htslib_kseq_h) $(htslib_hts_os_h)

misc/maq2sam-short.o: misc/maq2sam.c config.h version.h
	$(CC) $(CFLAGS) $(ALL_CPPFLAGS) -c -o $@ misc/maq2sam.c

misc/maq2sam-long.o: misc/maq2sam.c config.h version.h
	$(CC) $(CFLAGS) -DMAQ_LONGREADS $(ALL_CPPFLAGS) -c -o $@ misc/maq2sam.c


install: $(PROGRAMS) $(MISC_PROGRAMS)
	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(misc_bindir) $(DESTDIR)$(man1dir)
	$(INSTALL_PROGRAM) $(PROGRAMS) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(MISC_PROGRAMS) $(DESTDIR)$(misc_bindir)
	$(INSTALL_SCRIPT) $(MISC_SCRIPTS) $(DESTDIR)$(misc_bindir)
	$(INSTALL_MAN) doc/samtools*.1 misc/wgsim.1 $(DESTDIR)$(man1dir)


testclean:
	-rm -f test/*.new test/*.tmp test/*/*.new test/*/*.tmp test/*/*.tmp.*
	-cd test/dat && rm -f test_input_*.bam.bai
	-cd test/mpileup && rm -f FAIL-*.out* PASS-*.out* anomalous.[bc]*am indels.[bc]*am mpileup.*.[cs]*am mpileup.*.crai overlap50.[bc]*am expected/1.out xx#depth*.bam*

mostlyclean: testclean
	-rm -f *.o misc/*.o test/*.o test/*/*.o version.h $(LZ4OBJS)

clean: mostlyclean
	-rm -f $(PROGRAMS) libbam.a libst.a $(MISC_PROGRAMS) $(TEST_PROGRAMS)

distclean: clean
	-rm -f config.cache config.h config.log config.mk config.status
	-rm -f TAGS
	-rm -rf autom4te.cache

clean-all: clean clean-htslib


tags:
	ctags -f TAGS *.[ch] misc/*.[ch]


force:


.PHONY: all check clean clean-all distclean force install
.PHONY: lib mostlyclean print-version tags test testclean
