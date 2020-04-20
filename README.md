## Probe qPCR 
[![Build Status](https://travis-ci.com/pommevilla/probe_qpcr.svg?branch=master)](https://travis-ci.com/pommevilla/probe_qpcr)

A small script to perform *in silico* qPCR of TaqMan probes.

### Usage

Run the script on the command line by doing the following:

`python path/to/probe_qpcr.py taqman_primers.fasta targets.fasta`

The arguments are:

* `taqman_primers.fasta` - a fasta format of TaqMan primer sets with a `.F`, `.R`, or `.P` at the end of each sequence. For example:
```
>Primer_Set1.F
actgagtgctaa
>Primer_Set1.R
tgcacttagtaa
>Primer_Set1.P
ttttacgtagtaga
```

  
* `targets.fasta` - a file in fasta format

The results are output to stdout in a tab-delimited format. The fields are:

* `primer_name` - the name of the primer set resulting in the hit
* `target_name` - the name of the target sequence amplified
* `product_start` - the start location of the product
* `product_end` - the end location of the product
* `product_length` - the length of the product
* `product` - the product sequence


### To Do:

* Add support for ambiguous bases