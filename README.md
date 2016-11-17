# pSCAN - Fast and Exact Structural Graph Clustering (with overlaps)

The paper: *["pSCAN: Fast and Exact Structural Graph Clustering"](https://www.cse.unsw.edu.au/~ljchang/pdf/icde16-pscan.pdf) by Lijun Chang, Wei Li, Xuemin Lin,Lu Qin and Wenjie Zhang, ICDE'16*. The [presentation slides](https://www.cse.unsw.edu.au/~ljchang/pdf/icde16s-pscan.pdf) are also available.

Author: Dr. [Lijun Chang](https://www.cse.unsw.edu.au/~ljchang/) <ljchang@cse.unsw.edu.au>

This is the refactored version of the pScan, the unsafe code and most of the formatting were not modified. The I/O will be extended to be natively applicable for the [PyCaBeM](https://github.com/eXascaleInfolab/PyCABeM) clustering benchmark.  
Modified by Artem Lutov <artem@exascale.info>

# Usage
```
./pscan graph_directory epsilon mu [output]
```
For example
```
./pscan ./test 0.2 3 output
```
 `output` should be specified to save the resulting clustering into ${dir}/result-${epsilon}-${mu}.txt.

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
