Cram-Samtools
=============

This is a branch from samtools/samtools master upon which CRAM support has been added.
See http://www.ebi.ac.uk/ena/about/cram_toolkit for more information on CRAM.

The code here is primarily ripped out of Staden Package io_lib (aka libstaden-read),
which has a master copy in SVN at http://sourceforge.net/projects/staden/. As such it
contains vast tracks of code with duplicated functionality - my own hash tables, 
dynamic strings, malloc pools, etc.

In time we expect it to migrate to samtools khash and similar, but this is largely a
first pass to get things up and running quickly.

~ James Bonfield,
  Wellcome Trust Sanger Institute
