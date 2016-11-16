#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"

using namespace std;

class Graph {
private:
	string dir; //input graph directory
	Id n, m; //number of nodes and edges of the graph

	int eps_a2, eps_b2, miu; // eps_a2/eps_b2 = eps^2

	Id *pstart; //offset of neighbors of nodes
	Id *edges; //adjacent ids of edges
	Id *reverse; //the position of reverse edge in edges
	int *min_cn; //minimum common neighbor: -2 means not similar; -1 means similar; 0 means not sure; > 0 means the minimum common neighbor

	Id *pa;
	int *rank; //pa and rank are used for the disjoint-set data structure

	Id *cid; //cluster id

	Degree *degree;
	Degree *similar_degree; //number of adjacent edges with similarity no less than epsilon
	Degree *effective_degree; //number of adjacent edges not pruned by similarity

	vector<pair<int,int> > noncore_cluster;

public:
	Graph(const char *_dir=nullptr);
	~Graph();

	Graph(const Graph&)=delete;
	Graph(Graph&&)=delete;

	Graph& operator=(const Graph&)=delete;
	Graph& operator=(Graph&&)=delete;

	void read_graph(const char *filename=nullptr);
	void pSCAN(const char *eps_s, int miu);
		//eps_s and miu are the parameters (epsilon, miu) for the SCAN algorithm
	void cluster_noncore_vertices(int eps_a2, int eps_b2, int mu);
	void output(const char *eps_s, const char *miu);

private:
	Id binary_search(const Id *edges, Id b, Id e, Id val);
		//return the first pos, s.t. array[pos] >= val (may return e)
	int naive_similar_check(Id u, Id v, int eps_a2, int eps_b2);
	int similar_check(Id u, Id v, int eps_a2, int eps_b2);
	int similar_check_OP(Id u, Id idx, int eps_a, int eps_b);
	int check_common_neighbor(Id u, Id v, int c);
	int compute_common_neighbor_lowerbound(Id u,Id v,int eps_a2,int eps_b2);
	void prune_and_cross_link(int eps_a2, int eps_b2, int miu, int *cores, int &cores_e);

	Id find_root(Id u);
	void my_union(Id u, Id v);

	void get_eps(const char *eps_s);
};

#endif
