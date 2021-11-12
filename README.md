# *In silico* randomized antibody variant generator 

Scripts to both parse the positions of the complementarity determining regions (CDRs) from light and heavy chain sequences (string format) and generate variants through in-place N-mutations across the CDRs.

## Overview
This code package provides the user the ability to efficiently parse out the complementarity-determining regions (CDRs) of an antibody sequence given heavy and light chains, and design *in silico* large sets of antibody variants for a given seed sequence by introducing random mutations in the CDRs. To help overcome the scarcity of labeled antibody data in the public domain required to help advance the area of antibody engineering, one can first use this generator code to quickly produce a large set of antibody variant designs. These variant designs can then be used in downstream high-throughput methods (e.g. AlphaSeq, [1]) to generate large scale quantitative datasets, such as antibody-antigen binding interactions. 

Various methods exist to number the residues and determine the CDR positions in an antibody sequence, and each method suffers some limitations due to the varying nature of the CDR lengths. The Martin scheme attempts to overcome shortcomings in several prior approaches [2] and this rule set is implemented in the generator code package to extract the CDRs from their approximate positions in each chain. A set of permutations and combinations will be used to perform in silico randomization to introduce user-defined k mutations in the CDRs. 

[1] 	D. e. a. Younger, "High-throughput characterization of protein-protein interactions by reprogramming yeast mating.," Proceedings of the National Academy of Sciences, pp. 114(46): 12166-12171, 2017. 

[2] 	A. Martin, "How to identify the CDRs by looking at a sequence," UCL, [Online]. Available: http://www.bioinf.org.uk/abs/info.html#martinnum. [Accessed 13 October 2021].

## Running the code: 

usage: run_AbGen [-h] [--output_dir OUTPUT_DIR] config_file

positional arguments:
  config_file           Path to configuration file to load inputs.

optional arguments:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR
                        Path to save output files. 

## Config File 
Please refer to sample_config.json for an example of filled in config file. 

| Field name | Description | 
| :--------: | :---------: | 
| Filename   | Name of FASTA file to save to | 
| heavychain | Sequence of amino acids | 
| lightchain | Sequence of amino acids | 
| heavy_mutations | List of lists. [[<num_mutations>,<sample_size>],...] |
| light_mutations | List of lists. [[<num_mutations>,<sample_size>],...] |

If no sample size provided, a default sample size will be used instead.
If on GPU, override default mutation number of 6 in scripts/run_library_gen.py , 
otherwise large mutation number (i.e. >6) will cause performance issues on CPU.

## License
This work is licensed under a
[BSD-2-Clause License][bsd-2-clause].

[bsd-2-clause]: https://opensource.org/licenses/bsd-license.php

## Disclaimer
DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.

This material is based upon work supported by the United States Air Force under Air Force Contract No. FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the United States Air Force.

© 2021 Massachusetts Institute of Technology.

The software/firmware is provided to you on an As-Is basis

SPDX-License-Identifier: BSD-2-Clause
