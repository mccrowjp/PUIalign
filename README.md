# PUIalign
#### Phylogenetically Unambiguous Indel Alignment

Multiple sequence alignment of amino acid sequences often involve gaps representing relative inserts or deletions (indels), which can be useful as phylogenetic characters, but are usually avoided due to alignment ambiguity. By evaluating the stability of indels as syapomorphic characters we can produce a list of only those indels that are suitable for phylogenetic analysis. Alignment of indels must be robust with respect to alignment model parameters and the neighborhood of suboptimal alignments. See below for more information about this work published in 2009.

Citation
--------
McCrow JP, Alignment of phylogenetically unambiguous indels in Shewanella, J Comput Biol. 2009 Nov;16(11):1517-28. doi: 10.1089/cmb.2009.0188.

http://online.liebertpub.com/doi/abs/10.1089/cmb.2009.0188

Installation
------------

To install, simply download the source code and run the make utility to compile the C++ code.
```bash
make
```
This will creat the main executable file PUIalign.

Usage
-----

| File | Description |
|------|-------------|
| list_potential_indels.pl | finds indels to be evaluated
| PUIalign | main executable |
| show_indels.pl | displays the results of phylogenetically unambiguous indels
|  |  |
| align_allindels.cc | source code for main executable |
| indelfinder.cpp | functions for evaluating indels 
| indelfinder.h | definitions for indelfinder.cpp
| seq_align.cpp | fuctions for sequence alignment
| seq_align.h | definitions for seq_align.cpp
| progress.h | progress bar for C
| progress.pm | progress bar for Perl
| BLOSUM62.csv | amino acid substitution matrix
| [example/](./example/) | example data of bacterial taxa

```
Phylogenetically Unambiguous Indel Alignment (PUIalign) v1.0
Created by John P. McCrow (9/5/2007)

Usage: ./PUIalign [Indel List File] [Parameters File]
```

Dependencies
------------

* Perl (https://www.perl.org/get.html)
* C++ (GNU g++)
