# pSCAN - Fast and Exact Structural Graph Clustering (with overlaps)

The paper: *["pSCAN: Fast and Exact Structural Graph Clustering"](https://www.cse.unsw.edu.au/~ljchang/pdf/icde16-pscan.pdf) by Lijun Chang, Wei Li, Xuemin Lin,Lu Qin and Wenjie Zhang, ICDE'16*. The [presentation slides](https://www.cse.unsw.edu.au/~ljchang/pdf/icde16s-pscan.pdf) are also available.

Author: Dr. [Lijun Chang](https://www.cse.unsw.edu.au/~ljchang/) <ljchang@cse.unsw.edu.au>

This is the refactored version of the pScan, the unsafe code and original data structures were NOT modified. The I/O formats are extended to be natively applicable for the [PyCaBeM](https://github.com/eXascaleInfolab/PyCABeM) clustering benchmark.  
Modified by Artem Lutov <artem@exascale.info>

# Usage
```
$ ./pscan -h
pSCAN 0.2

Clusters input network considering overlaps and building Exact Structural Graph

Usage: pSCAN [OPTIONS]... [input_network]...

input_network  - the input graph specified as either a file in the NSL format,
or a directory containing "b_degree.bin" and "b_adj.bin" binary files. If
the format is not specified explicitly then NSL file is expected and whether it
is NSA or NSE is identified by the file extension.

  -h, --help           Print help and exit
  -V, --version        Print version and exit
  -f, --format=ENUM    format of the input graph  (possible values="BIN",
                         "NSA", "NSE")
  -e, --epsilon=FLOAT  similarity threshold (typically E [1e-4, 0.5] )
                         (default=`0.2')
  -m, --mu=INT         size threshold  (default=`3')
  -l, --legacy         output clustering in the legacy pSCAN format instead of
                         the standard NSL  (default=off)
  -o, --output=STRING  output file if the resulting clustering should be saved
```
For example
```
./pscan -e=0.2 -m=3 -l -o test/clusters.txt -f=BIN test
```
The input network to be clustered is specified either in NSL (nsa/nse) format by 2 BINARY files:

1. NSL format (nsa - arcs, directed network; nse - edges, undirected network) specifies network links is each line of the file as node ids separated by the space delimiter:

	```
	# Example Network .nse (edges - undirected)
	# Nodes: 3  Edges: 3   Weighted: 0
	# Note that the number of links corresponds to the number of payload lines in the file
	0 1
	# Empty lines and comments are allowed
	0 2
	2 1
	```
2. Binary format:
  - Specification of the network properties and nodes (vertices) degrees, `b_degree.bin`:
  
	```
	4  // <id_len_in_bytes>
	13 // <nodes_number> (n)
	48 // <arcs_number> (= 2 * edges_num = 2*m)
	4  // <node1_degree>
	3  // <node2_degree>
	...
	3  // <noden_degree>
	```
  - Spase separated List of neigbours for each node in a new line, `b_adj.bin`:
  
	```
	1 2 3 4 // neighbours of vertex 0
	0 2 3   // neighbours of vertex 1
	...
	```
