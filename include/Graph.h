#ifndef _GRAPH_H_
#define _GRAPH_H_

#ifdef _MSC_VER
	#define _CRT_SECURE_NO_WARNINGS
#endif

#include <string>
#include <vector>
#include <unordered_map>
#include "cmdline.h"


namespace pscan {

typedef unsigned  Id;
typedef int  Degree;  // ATTENTION: degrees are specified as int in the input file

using std::string;
using std::vector;
using std::unordered_map;
using std::pair;

class Graph {
	enum_format inpfmt;  // Input format
	string input;  // Input network (graph)
	Id n, m; //number of nodes and ARCS (2*edges) of the graph

	const float eps;
	int eps_a2, eps_b2, miu; // eps_a2/eps_b2 = eps^2

	Id *pstart; //offset of neighbors of nodes
	Id *edges; //adjacent ids of edges
	Id *reverse; //the position of reverse edge in edges
	int *min_cn; //minimum common neighbor: -2 means not similar; -1 means similar; 0 means not sure; > 0 means the minimum common neighbor

	Id *pa;
	int *rank; //pa and rank are used for the disjoint-set data structure

	Id *cid; //cluster ids

	Degree *degree;
	Degree *similar_degree; //number of adjacent edges with similarity no less than epsilon
	Degree *effective_degree; //number of adjacent edges not pruned by similarity

	vector<pair<int,int> > noncore_cluster;

	unordered_map<Id, Id>  ieids;  // Map from internal to external ids of nodes
public:
	Graph(float aeps, int amiu, const char *ainput, enum_format ainpfmt);
	~Graph();

	Graph(const Graph&)=delete;
	Graph(Graph&&)=delete;

	Graph& operator=(const Graph&)=delete;
	Graph& operator=(Graph&&)=delete;

	void load();
	void pSCAN();
		//eps_s and miu are the parameters (epsilon, miu) for the SCAN algorithm
	void cluster_noncore_vertices(int mu);
	void output(const char *outfile, bool lgcfmt=false);
protected:
	void loadBinary();
	void loadNSL();  // Load  network (graph) in NSE/A format

	void saveLegacy(const char *outfile);
	void saveCNL(const char *outfile);
private:
	Id binary_search(const Id *edges, Id b, Id e, Id val);
		//return the first pos, s.t. array[pos] >= val (may return e)
	int naive_similar_check(Id u, Id v);
	int similar_check(Id u, Id v);
	int similar_check_OP(Id u, Id idx);
	int check_common_neighbor(Id u, Id v, int c);
	int compute_common_neighbor_lowerbound(Id u,Id v);
	void prune_and_cross_link(int miu, int *cores, int &cores_e);

	Id find_root(Id u);
	void my_union(Id u, Id v);
};

}  // pscan

std::string to_string(enum_format fmt);  //!< Convert to string
void tolower(char* text);  //!< Convert string to lower case

#endif
