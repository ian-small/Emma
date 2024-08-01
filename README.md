# EMMA: Metazoan Mitochondrial Annotator
(Currently EMMA is optimised for vertebrates, especially fish, but the intention is to make it broadly useful for most metazoan mitochindrial genomes)

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

(replace ~/github with whatever path you've cloned the Emma repo to)

This will install all the Julia packages that Emma needs, and the models it needs to annotate mitochondrial genomes.

You can now quit the Julia REPL with Ctrl-D

## Running EMMA
`julia --project=~/github/Emma ~/github/Emma/src/command.jl --help`                                             
Usage: Emma/src/command.jl [options] <FASTA_file>

Note: Use consistant inputs/outputs. If you wish
to annotate a directory of fasta files, ensure that
the output parameters are also directories.

Positional Arguments:
FASTA_file
      file/dir for fasta input
      (Type: String, Required)

Option Arguments:
--gff
      file/dir for gff output
      (Type: String)
--gb
      file/dir for gb output
      (Type: String)
--fa
      file/dir for fasta output. Use this argument if you wish
      annotations to begin with tRNA-Phe
      (Type: String)
--svg
      file/dir for svg output
      (Type: String)







