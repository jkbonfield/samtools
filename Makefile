# The default version string in bam.h and bcftools/bcf.h can be overriden directly
#   make VERSION="-DVERSION='\\\"my-version\\\"'"
# or using the git-stamp rule
#   make git-stamp
VERSION=

CC=			gcc
CFLAGS=		-g -Wall $(VERSION) -O2
#LDFLAGS=		-Wl,-rpath,\$$ORIGIN/../lib
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=1 -DHAVE_LIBCURL -DSAMTOOLS=1
KNETFILE_O=	knetfile.o
LOBJS=		bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o	\
			bam_pileup.o bam_lpileup.o bam_md.o razf.o faidx.o bedidx.o \
			$(KNETFILE_O) bam_sort.o sam_header.o bam_reheader.o kprobaln.o bam_cat.o $(COBJS)
COBJS=  io_lib/cram_codecs.o \
	io_lib/cram_encode.o \
	io_lib/cram_io.o \
	io_lib/cram_decode.o \
	io_lib/cram_index.o \
	io_lib/cram_stats.o \
	io_lib/cram_samtools.o \
	io_lib/sam_header.o \
	io_lib/hash_table.o \
	io_lib/jenkins_lookup3.o \
	io_lib/vlen.o \
	io_lib/zfio.o \
	io_lib/mFILE.o \
	io_lib/md5.o \
	io_lib/open_trace_file.o \
	io_lib/pooled_alloc.o \
	io_lib/string_alloc.o \
	io_lib/dstring.o \
	io_lib/files.o \
	io_lib/error.o \
	io_lib/xalloc.o \
	io_lib/thread_pool.o
AOBJS=		bam_tview.o bam_plcmd.o sam_view.o \
			bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o \
			bamtk.o kaln.o bam2bcf.o bam2bcf_indel.o errmod.o sample.o \
			cut_target.o phase.o bam2depth.o padding.o bedcov.o bamshuf.o \
			bam_tview_curses.o bam_tview_html.o
PROG=		samtools
INCLUDES=	-I.
SUBDIRS=	. bcftools misc
LIBPATH=
LIBCURSES=	-lcurses # -lXCurses


.SUFFIXES:.c .o
.PHONY: all lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all-recur lib-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" $$target || exit 1; \
			cd $$wdir; \
		done;

all:$(PROG)

git-stamp:
		make VERSION="-DVERSION='\\\"`git describe --always --dirty`\\\"'"

.PHONY:all lib clean cleanlocal
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

samtools:lib-recur $(AOBJS)
		$(CC) $(CFLAGS) -o $@ $(AOBJS) $(LDFLAGS) libbam.a -Lbcftools -lbcf $(LIBPATH) $(LIBCURSES) -lm -lz -lpthread -lcurl

razip:razip.o razf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ $^ -lz

bgzip:bgzip.o bgzf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ $^ -lz -lpthread

bgzf.o:bgzf.c bgzf.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DBGZF_CACHE $(INCLUDES) bgzf.c -o $@

razip.o:razf.h
bam.o:bam.h razf.h bam_endian.h kstring.h sam_header.h
sam.o:sam.h bam.h
bam_import.o:bam.h kseq.h khash.h razf.h
bam_pileup.o:bam.h razf.h ksort.h
bam_plcmd.o:bam.h faidx.h bcftools/bcf.h bam2bcf.h
bam_index.o:bam.h khash.h ksort.h razf.h bam_endian.h
bam_lpileup.o:bam.h ksort.h
bam_tview.o:bam.h faidx.h bam_tview.h
bam_tview_curses.o:bam.h faidx.h bam_tview.h
bam_tview_html.o:bam.h faidx.h bam_tview.h
bam_sort.o:bam.h ksort.h razf.h
bam_md.o:bam.h faidx.h
sam_header.o:sam_header.h khash.h
bcf.o:bcftools/bcf.h
bam2bcf.o:bam2bcf.h errmod.h bcftools/bcf.h
bam2bcf_indel.o:bam2bcf.h
errmod.o:errmod.h
phase.o:bam.h khash.h ksort.h
bamtk.o:bam.h

faidx.o:faidx.h razf.h khash.h
faidx_main.o:faidx.h razf.h


libbam.1.dylib-local:$(LOBJS)
		libtool -dynamic $(LOBJS) -o libbam.1.dylib -lc -lz -lcurl

libbam.so.1-local:$(LOBJS)
		$(CC) -shared -Wl,-soname,libbam.so -o libbam.so.1 $(LOBJS) -lc -lz -lcurl

dylib:
		@$(MAKE) cleanlocal; \
		case `uname` in \
			Linux) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libbam.so.1-local;; \
			Darwin) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libbam.1.dylib-local;; \
			*) echo 'Unknown OS';; \
		esac


cleanlocal:
		rm -fr gmon.out *.o a.out *.exe *.dSYM razip bgzip $(PROG) *~ *.a *.so.* *.so *.dylib

clean:cleanlocal-recur
	-rm io_lib/*.o
