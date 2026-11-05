# samtools sort2 — experimental pipelined, parallel sort

(This version is "try19", matched to samtools develop: the try8 design, plus per-phase `--stats`,
early deletion of temp files, merge ranges in whole rounds of workers,
`--write-index`, and every samtools sort order: `-t`, `-M`/`-K` and
`--template-coordinate`.  Later experiments, byte-sized write chunks
and in-memory output buffers, were slower on real long-read data and were
dropped.  `htslib-mem-write/` is an independent htslib bug fix, not used by
sort2.)

`samtools sort2` is a new subcommand for benchmarking side by side with
`samtools sort`.  It supports every `sort` order in samtools develop,
coordinate (default), `-n` (natural name), `-N` (byte name), `-t TAG`, `-M`
(with `-K`, `-R`, `-H`, `-I`, `-w`) and `--template-coordinate`, and
produces the same record order as `sort`; on a
73M-record collated BAM a coordinate sort was byte-identical to `sort`'s.

## The other sort orders

These follow samtools develop's `bam_sort.c` exactly, including how the options
combine (the last of `-n`, `-M` and `--template-coordinate` wins, and `-t`
then makes it a tag sort with coordinate or name as the secondary key) and
the `@HD` SO/GO/SS tags each order sets.  `--write-index` is ignored with a
warning for name, tag and template orders, as in `sort`.

- **`-t TAG`** (`bam1_cmp_by_tag()`): reads without the tag first; values of
  different types by their normalised type (A < B < H/Z < numbers, with
  integers and floats compared numerically); then the value; then coordinate
  order, or name order with `-n`/`-N`.  Encoded as a byte-comparable key.
- **`-M`, `-K INT`, `-R`, `-H`, `-I FILE`, `-w INT`** (`minhash()`,
  `build_minhash_index()`, `minhash_with_idx[_squash]()`,
  `worker_minhash()`, ported verbatim): unplaced reads are clustered by
  minimiser (both strands unless `-R`; homopolymers squashed with `-H`), or
  with `-I` by the reference position of their best indexed k-mer;
  reverse complemented (and 0x10 flipped) when the reverse strand wins.  The
  key is cached in pos/mpos/isize, so it survives temp files, and reset on
  output.  The index is read-only once built and shared by the threads.
  samtools' complement table is copied verbatim so output is identical,
  quirks included (see the samtools bug list).  The key is numeric, with a
  third field and one more stable radix pass.
- **`--template-coordinate`** (`template_coordinate_key()`): unclipped 5'
  ends (soft clips only; the mate's from its MC tag), strands, library (the
  read group's LB), CB, MI ignoring a trailing `/x`, name, and which end of
  the pair this is.
  Encoded as a byte-comparable key.  Reads that are paired with a mapped
  mate but have no MC tag are an error, as in `sort`.  `bam.h` provides
  `unclipped_start()` and friends, as for `bam_sort.c`.

Orders may differ from `sort` only in pathological cases: integer tags
beyond 2^53 tied with an equal float, NaN float tags (for which `sort`'s
comparator is not a consistent order), and non-ASCII MI tags (which `sort`
compares as signed `char`).

## Integrating into a samtools checkout

Copy `bam_sort2.c`, `sort2_engine.c` and `sort2_engine.h` into the samtools
source directory.  In **bamtk.c**, declare `int bam_sort2(int argc, char
*argv[]);` next to `bam_sort`, add
`else if (strcmp(argv[1], "sort2") == 0) ret = bam_sort2(argc-1, argv+1);`
to the dispatch chain, and a usage line under "File operations".  In the
**Makefile**, add `bam_sort2.o sort2_engine.o` to `AOBJS` with:

```make
bam_sort2.o: bam_sort2.c config.h $(htslib_sam_h) $(htslib_bgzf_h) $(htslib_hfile_h) $(htslib_hts_h) $(htslib_thread_pool_h) sam_opts.h samtools.h sort2_engine.h
sort2_engine.o: sort2_engine.c sort2_engine.h
```

## Options of note

- `-m` is per thread, as in `sort` (default 768M).
- `--serial-merge` merges on one thread through `sam_write1()` with htslib's
  threaded output, as `sort` does.  For comparison.
- `--stats` prints a per-phase breakdown, including user/system CPU and
  syscall read/write bytes per phase, temp read-back amplification, the time
  spent deleting temp files, and cleanup time after the output is closed;

## `--write-index`

As `samtools sort --write-index`: a CSI index (`min_shift` 14) named
`<out>.csi`, or `<out>.crai` for CRAM; coordinate order and an output file
are required (the `fn##idx##fnidx` form is not supported yet).  On the
serial paths htslib builds it (`sam_idx_init()` / `sam_idx_save()`).  In
the parallel BAM merge, records go to per-worker scratch files, so sort2
builds it the way htslib does for indexed BAM writes: `bgzf_flush_try()`
before each record, noting the scratch virtual offset, then the offset
after it; when a range's turn comes and its bytes reach the output, those
offsets are translated to output positions (the range's first byte in the
output plus the distance into the scratch file), and the worker holding the
turn adds the entries in output order with `hts_idx_amend_last()` and
`hts_idx_push()`.  `hts_idx_finish()` is given the end of the records and
the index is saved after the output is closed.  `--stats` reports the
output close time (which includes htslib's flush and sync) and which way
the index was built.

## Design

**Keys.** Coordinate order is (tid', ((pos+1) << 1) | rev), with unmapped
reads after every real tid, covering the full 64-bit position range SAM
allows.  A batch whose positions fit in 32 bits (all BAM data) packs this
into one 64-bit word and sorts it with a single parallel LSD radix sort
(constant digits skipped); otherwise two stable passes are used.  Name
order encodes the name so `memcmp` order equals `strnum_cmp` order, plus a
READ1/READ2 tie byte, and uses a parallel MSD radix sort past the batch's
common prefix.  Only 16-byte {key, ref} pairs move; records are gathered
once.  Ties keep input order.

**Pipeline.** Two batches, each half the memory: the main thread reads into
one while a sorter thread sorts the other and writes it as one sorted run.
Run writing is parallel: the sorted batch is cut into chunks of about
128 KB (a record count from the batch's average record size, at least 16
records) and groups of 16 chunks; workers claim groups in order,
serialise and compress them with `bgzf_compress()`, reserve file space under
a lock and `pwrite()` outside it.  Each chunk starts a BGZF block, giving a
sparse in-memory index of first keys.  The last batch stays in memory; if
it is the only one, nothing is written to disk.  If the input is already in
order the runs simply don't overlap.  Runs are merged early (serially, all
into one) once there are 512, or fewer if descriptors are scarce.

**Final merge.** The key space is cut into ranges at quantiles of the run
indexes; a cut point is a (key, run) pair so equal keys split without
changing tie order.  Ranges are sized to keep seek overhead per range low
and buffered output to about a quarter of the memory.  The range count is then
rounded up to a whole number of rounds of workers (at least four), so no
workers sit idle while a part-filled last round finishes.  Each worker seeks
every run to its range and merges it with its own loser tree.  Output:

- **Default, SAM and BAM:** each merge worker writes its ranges with
  `sam_write1()` into its own scratch `htsFile`, in the output's format, so
  htslib does all the encoding.  The scratch file is a `memfd` on Linux
  (else an unlinked temporary file); its header is written once and
  discarded.  Output order is a token: a worker may write only when
  `next_write` equals its range, and nobody else writes until it passes the
  turn on, so the copy needs no lock.  The scratch file is `mmap()`ed,
  written with one `hwrite()` and truncated.  The worker whose range is
  next streams its output every 16k records instead of buffering the range.
- **Serial:** CRAM output, `--serial-merge`, `-@` below 2, or too few
  descriptors.

Every merge worker opens every run; the worker count is limited to fit the
descriptor limit, which sort2 raises towards the hard limit if needed.
Temp files are unlinked as soon as nothing needs to open them by name (in
the parallel merge, once the last worker has opened them all), so their
space is released as the merge finishes rather than at exit.  `--stats`
reports the time spent deleting them and the cleanup time after the output
is closed.

## samtools sort bugs found along the way (develop)

- `reverse_complement()`'s table maps `Y`->`Y` and `y`->`y` (should be
  `R`/`r`), so `-M` corrupts Y bases in reverse-complemented reads; it also
  turns `=` into `N`.  sort2 copies the table verbatim to match.
- `minhash()`, `minhash_with_idx()` and `_squash()` use `(1L<<(2*kmer))-1`,
  undefined for `-K` >= 16 where `long` is 32 bits (Windows); sort2 uses a
  64-bit constant, identical on LP64.
- Minor: MI is compared as signed `char` while names/libraries use `strcmp`;
  NaN float tags make `bam1_cmp_by_tag()` inconsistent; `-M`'s 65535-pos
  tie-break saturates on reads over 64 kbp; other per-base tags (e.g. OQ)
  are not reversed by the minimiser's reverse complement.

## Testing

`test/run_tests.sh` builds against `test/mock/`, a functional stand-in for
the htslib calls used (real BGZF via zlib), and runs:

- engine unit tests (`test/test_engine.c`): name encoding against a verbatim
  `strnum_cmp`, radix and name sorts against a reference sort, loser tree;
- 278 end-to-end cases under ASan/UBSan (`test/test_e2e.py`, an independent
  oracle that also validates every BGZF block and decodes BAM): all orders, short and long reads (records up to 200 KB),
  BAM and text output (compressed and `-u`), parallel and serial merges,
  with and without forced head-of-line streaming, in-memory, many runs and
  intermediate merges, already-sorted input, heavy ties, 64-bit positions;
- `--write-index` cases that map every index entry through the output's real
  BGZF block layout onto its record;
- 51 cases for `-t`, `-M` (`-K -R -H -I -w`) and `--template-coordinate`
  (`test/test_orders.py`): whole records, including minimiser reverse
  complements, compared with an independent Python port of samtools
  develop's comparators, minimiser and index code, on every merge path,
  plus option interactions and `@HD` tags;
- ThreadSanitizer over the engine and the full pipeline.

On real htslib, check first: `hts_hopen()`/`hdopen()` scratch files,
`htell()`, `hseek()` and truncation between ranges, and that the scratch
EOF/header handling gives byte-identical output to `--serial-merge`.

## Benchmarking

```
SAMTOOLS=./samtools MEM_MB=768 TMPDIR=/local/tmp ./bench_sort.sh big.bam 8 16 32
```

Both tools get the same per-thread `-m`; the script reports wall, user and
system time and peak RSS, sort2's phase breakdown, and checks that record
order matches.  Comparing default and `--serial-merge` on
the same input shows what the parallel merge is worth.

## Parked / next steps

- Parallel BAM parsing on input (an htslib job).
- Tiered parallel intermediate merges, and skipping runs whose key span
  cannot overlap a range (both matter only with many runs).
- CRAM output through the htsFile-buffer merge (container concatenation,
  record counters, shared references).
