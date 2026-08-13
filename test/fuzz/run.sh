#!/bin/sh

# We deliberately run it many times over as the tools may leak memory or
# fail to close files, and we wish to restart again repeatedly to solve this.
#
# Also if it finds an error and aborts, we restart again too.
(for i in `seq 1 100000`;do ASAN_OPTIONS=allow_addr2line=1:abort_on_error=1:detect_leaks=1 ./fuzz_mpileup corpus_in/ -rss_limit_mb=12560 -timeout=1 -max_total_time=10;done) 2>&1; # | tee /tmp/fuzz.out
