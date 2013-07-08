#ifndef _CRAM_SAMTOOLS_H_
#define _CRAM_SAMTOOLS_H_

/* Samtools compatible API */
#define bam_blk_size(b)  ((b)->data_len)
#define bam_set_blk_size(b,v) ((b)->data_len = (v))

#define bam_ref(b)       (b)->core.tid
#define bam_pos(b)       (b)->core.pos
#define bam_mate_pos(b)  (b)->core.mpos
#define bam_mate_ref(b)  (b)->core.mtid
#define bam_ins_size(b)  (b)->core.isize
#define bam_seq_len(b)   (b)->core.l_qseq
#define bam_cigar_len(b) (b)->core.n_cigar
#define bam_flag(b)      (b)->core.flag
#define bam_bin(b)       (b)->core.bin
#define bam_map_qual(b)  (b)->core.qual
#define bam_name_len(b)  (b)->core.l_qname
#define bam_name(b)      (b)->data

#define bam_qual(b)      bam1_qual((b))
#define bam_seq(b)       bam1_seq((b))
#define bam_cigar(b)     bam1_cigar((b))
#define bam_aux(b)       bam1_aux((b))
#define bam_seqi(s,i)    bam1_seqi((s),(i))

#define bam_dup(b)       bam_dup1((b))

enum cigar_op {
    BAM_CMATCH=0,
    BAM_CINS=1,
    BAM_CDEL=2,
    BAM_CREF_SKIP=3,
    BAM_CSOFT_CLIP=4,
    BAM_CHARD_CLIP=5,
    BAM_CPAD=6,
    BAM_CBASE_MATCH=7,
    BAM_CBASE_MISMATCH=8
};

#include "bam.h"

typedef bam1_t bam_seq_t;

// Sorry this is a bit loopy
#include "io_lib/cram.h"

int cram_get_bam1_seq(cram_fd *fd, bam1_t *b);
bam_header_t *cram_header_to_bam(SAM_hdr *h);
SAM_hdr *bam_header_to_cram(bam_header_t *h);

int bam_construct_seq(bam_seq_t **bp, size_t extra_len,
		      const char *qname, size_t qname_len,
		      int flag,
		      int rname,      // Ref ID
		      int pos,
		      int end,        // aligned start/end coords
		      int mapq,
		      uint32_t ncigar, const uint32_t *cigar,
		      int mrnm,       // Mate Ref ID
		      int mpos,
		      int isize,
		      int len,
		      const char *seq,
		      const char *qual);

#endif /* _CRAM_SAMTOOLS_H_ */
