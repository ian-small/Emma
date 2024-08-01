# EMMA: Metazoan Mitochondrial Annotator
(Currently EMMA is optimised for vertebrates, especially fish, but the intention is to make it broadly useful for most metazoan mitochondrial genomes)

## Dependencies
EMMA relies on two amazing software packages from Sean Eddy's lab, [HMMER](http://hmmer.org) and [Infernal](http://eddylab.org/infernal/). Many programs from these packages must be accessible via your $PATH for Emma to function. Follow the installation instructions on the [HMMER](http://hmmer.org) and [Infernal](http://eddylab.org/infernal/) websites to install these packages.
EMMA is written in [Julia](https://julialang.org); follow the [download and installation instructions](https://julialang.org/downloads/).

## Installation
Clone this github repo:

`git clone https://github.com/ian-small/Emma.git`

Tell Julia to treat the repo as a Julia package:

`julia`

`]`

`dev '~/github/Emma'`

`instantiate`

(replace ~/github with whatever path you've cloned the Emma repo to)

This will install all the Julia packages that Emma needs, and the models it needs to annotate mitochondrial genomes.

You can now quit the Julia REPL with Ctrl-D

## Running EMMA
`julia --project=~/github/Emma ~/github/Emma/src/command.jl --help`                                             
Usage: Emma/src/command.jl [options] <FASTA file or directory>

Optional Arguments:

--gff file/dir for .gff output
            
--tbl file/dir for .tbl output (for GenBank submissions)
      
--fa file/dir for .fasta output; use this argument if you wish to rotate the genome to begin with tRNA-Phe
            
--svg file/dir for .svg output
      
Note: Use consistant inputs/outputs. If you wish to annotate a directory of fasta files, ensure that the output options are also directories.

## Examples

To annotate a single genome in .gff format:

`julia --project=~/github/Emma ~/github/Emma/src/command.jl --gff my_genome.gff my_genome.fasta`

To annotate many fasta files in the directory 'my_genomes' and save the generated .tbl files in the same directory:

`julia --project=~/github/Emma ~/github/Emma/src/command.jl --tbl my_genomes my_genomes`
      







