dbgfm
=====

## Quick info

This is repo collects several different pieces of code to simplify the benchmarking of dbgfm:

1. The [original dbgfm package](https://github.com/jts/dbgfm)
2. [bwtdisk](https://people.unipmn.it/manzini/bwtdisk/), which is required for index construction
3. dbgfm2, another entrypoint to the data structure, developed by [Giulio Pibiri](https://github.com/jermp), which allows to test the impact of the sampling parameter for suffix array (used, e.g., in the [SSHash paper](https://doi.org/10.1093/bioinformatics/btac245))


## Usage

1. Clone and compile the repo:

	git clone https://github.com/karel-brinda/dbgfm && cd dbgfm && make -j

2. Build the FM index for your SPSS representation (unitigs / simplitigs; matchtigs probably not currently supported):

	./run_bwtdisk.sh your.fa


3. Run queries (note: the sampling tables are computed on the fly upon the execution)



## Old readme


An FM-index representation of a de Bruijn graph.

The code in this repository is a stand-alone version of the FM-index from SGA (github/jts/sga).
It is licensed under GPLv3.

### Compiling

The code has no dependencies and should build by just running:

	make

### Testing

To test the code is functioning correctly, you can run:

	make test

This will download human chromosome 20, index it with SGA then perform test queries using dbgfm.
You will need to modify the Makefile to point to your version of SGA.
This requires [bwtdisk](http://people.unipmn.it/manzini/bwtdisk/) is installed.

### API

A simple API for querying the structure of the de Bruijn graph is provided. See [dbg_query.h](/dbg_query.h/) and the [test driver](main.cpp).
