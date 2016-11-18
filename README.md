# pSCAN - Fast and Exact Structural Graph Clustering (with overlaps)

The paper: *["pSCAN: Fast and Exact Structural Graph Clustering"](https://www.cse.unsw.edu.au/~ljchang/pdf/icde16-pscan.pdf) by Lijun Chang, Wei Li, Xuemin Lin,Lu Qin and Wenjie Zhang, ICDE'16*. The [presentation slides](https://www.cse.unsw.edu.au/~ljchang/pdf/icde16s-pscan.pdf) are also available.

Author: Dr. [Lijun Chang](https://www.cse.unsw.edu.au/~ljchang/) <ljchang@cse.unsw.edu.au>

This is the refactored version of the pScan, the unsafe code and most of the formatting were not modified. The I/O will be extended to be natively applicable for the [PyCaBeM](https://github.com/eXascaleInfolab/PyCABeM) clustering benchmark.  
Modified by Artem Lutov <artem@exascale.info>

# Usage
```
$ ./pscan -h
pSCAN 0.2

Clusters input network considering overlaps and building Exact Structural Graph

Usage: pSCAN [OPTIONS]... [input_network]...

input_network  - graph specified as either a file in the CNL format, or a
directory containing "b_degree.bin" and "b_adj.bin"

  -h, --help           Print help and exit
  -V, --version        Print version and exit
  -d, --dir            input network (graph) is specified by the
                         "b_degree.bin" and "b_adj.bin" files located in
                         the specified directory  (default=off)
  -e, --epsilon=FLOAT  epsilon  (default=`0.2')
  -m, --mu=INT         epsilon  (default=`3')
  -l, --legacy         output clustering in the legacy pSCAN format instead of
                         the standard CNL  (default=off)
  -o, --output=STRING  output file if the resulting clustering should be saved
```
For example
```
./pscan -d -e=0.2 -m=3 -l -o test/clusters.txt test
```
 `--output` should be specified to save the resulting clustering into the specified directory (${dir}/result-${epsilon}-${mu}.txt by default).

The input network to be clustered is specified by 2 BINARY files:
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
