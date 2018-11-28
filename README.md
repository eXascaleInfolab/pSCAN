# pSCAN - Fast and Exact Structural Graph Clustering (with overlaps)

The paper: *["pSCAN: Fast and Exact Structural Graph Clustering"](https://www.cse.unsw.edu.au/~ljchang/pdf/icde16-pscan.pdf) by Lijun Chang, Wei Li, Xuemin Lin,Lu Qin and Wenjie Zhang, ICDE'16*. The [presentation slides](https://www.cse.unsw.edu.au/~ljchang/pdf/icde16s-pscan.pdf) are also available.

Author: Dr. [Lijun Chang](https://www.cse.unsw.edu.au/~ljchang/) <ljchang@cse.unsw.edu.au>

This is the refactored version of the pScan, the unsafe code and original data structures were NOT modified. The I/O formats are extended to be natively applicable for the [PyCaBeM](https://github.com/eXascaleInfolab/PyCABeM) clustering benchmark.  
Extended by Artem Lutov <artem@exascale.info>

## Content
- [Deployment](#deployment)
	- [Dependencies](#dependencies)
	- [Compilation](#compilation)
- [Usage](#usage)
  - [Input](#input)
  - [Output](#output)
- [Related Projects](#related-projects)

# Deployment

## Dependencies
There no any dependencies for the execution or compilation.  
However, to extend the input options and automatically regenerate the input parsing,
[*gengetopt*](https://www.gnu.org/software/gengetopt) application should be installed: `$ sudo apt-get install gengetopt`.

## Compilation
Just execute `$ make [release | debug]`.  
To update/extend the input parameters just modify `args.ggo` and run `GenerateArgparser.sh` (calls `gengetopt`).

# Usage
```
$ ./pscan -h
pSCAN 0.2

Clusters unweighted undirected input network (graph) considering overlaps and
building Exact Structural Graph

Usage: pSCAN [OPTIONS]... [input_network]...

input_network  - the input graph specified as either a file in the NSL format,
or a directory containing "b_degree.bin" and "b_adj.bin" binary files. If
the format is not specified explicitly then NSL file is expected and whether it
is NSA or NSE is identified by the file extension or header.

  -h, --help           Print help and exit
  -V, --version        Print version and exit
  -f, --format=ENUM    format of the input graph  (possible values="BIN",
                         "NSA", "NSE")
  -e, --epsilon=FLOAT  similarity threshold (typically E [0.05, 0.95]), where a
                         higher value corresponds to the higher granularity and
                         more careful differentiation, a lower value yields
                         large clusters  (default=`0.35')
  -m, --mu=INT         size threshold  (default=`3')
  -l, --legacy         output clustering in the legacy pSCAN format instead of
                         the standard CNL (default=off)
  -o, --output=STRING  output file if the resulting clustering should be saved
```
For example
```
./pscan -e 0.7 -o graph-e7.cnl graph.nse

./pscan -e 0.2 -m 3 -l -o test/clusters.txt -f=BIN test
```

## Input
The undirected unweighted input network to be clustered is specified either in the NSL (nsa/nse) format or by the 2 BINARY files:

1. NSL format (nsa - arcs, directed network; nse - edges, undirected network) specifies network links in each line of the file as node ids separated by the space delimiter with optional `#` line comments and an optional header:

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
  - Space separated List of neighbors for each node in a new line, `b_adj.bin`:

	```
	1 2 3 4 // neighbors of vertex 0
	0 2 3   // neighbors of vertex 1
	...
	```

## Output
The CNL (clusters nodes list) output is a standard format, generalization of the Stanford SNAP ground-truth communities. It is an input format for various NMI-evaluation algorithms. Each line of the file corresponds to the single resulting cluster, where member nodes are specified separated by the space/tab with optional share. For example:
```
# Clusters: 2, Nodes: 5, Fuzzy: 0
0
1 3 2 4
```

# Related Projects
- [xmeasures](https://github.com/eXascaleInfolab/xmeasures)  - Extrinsic quality (accuracy) measures evaluation for the overlapping clustering on large datasets: family of mean F1-Score (including clusters labeling), Omega Index (fuzzy version of the Adjusted Rand Index) and standard NMI (for non-overlapping clusters).
- [GenConvNMI](https://github.com/eXascaleInfolab/GenConvNMI) - Overlapping NMI evaluation that is compatible with the original NMI (unlike the `onmi`).
- [OvpNMI](https://github.com/eXascaleInfolab/OvpNMI) - Another method of the NMI evaluation for the overlapping clusters (communities) that is not compatible with the standard NMI value unlike GenConvNMI, but it is much faster and yields exact results unlike probabilistic results with some variance in GenConvNMI.
- [Clubmark](https://github.com/eXascaleInfolab/clubmark) - A parallel isolation framework for benchmarking and profiling clustering (community detection) algorithms considering overlaps (covers).
- [ExecTime](https://bitbucket.org/lumais/exectime/)  - A lightweight resource consumption (RSS RAM, CPU, etc.) profiler.
