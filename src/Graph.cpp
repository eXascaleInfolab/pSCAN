#include <cstdlib>
#include <cstdio>
#include <cstring>  // strtok
#include <cmath>
#include <cstring>
#include <algorithm>
#include <sstream>
#include <cctype>  // tolower
//#include <iostream>
#include <set>
// #include <unordered_set>
#include <stdexcept>
//#include <system_error>
#include <fstream>
#include <cerrno>
#include <ctime>
#include <cassert>

#ifdef __unix__
#include <sys/time.h>  // gettimeofday
#endif

#include "Graph.h"

using std::set;
using std::invalid_argument;
using std::domain_error;
//using std::system_error;
using std::ios_base;
using std::ifstream;
using std::swap;
using std::to_string;
using namespace pscan;


// Accessory Functions ---------------------------------------------------------
string to_string(enum_format fmt)
{
	string res;
	switch(fmt) {
	case format_arg_BIN:
		res = "BINARY";
		break;
	case format_arg_NSA:
		res = "NSA";
		break;
	case format_arg_NSE:
		res = "NSE";
		break;
	case format__NULL:
	default:
		res = "UNDEFINED";
	}

	return res;
}

void tolower(char* text)
{
	if(!text)
		return;
	while(*text)
		*text++ = tolower(*text);
}

// Graph Methods ---------------------------------------------------------------
Graph::Graph(float aeps, int amiu, const char *ainput, enum_format ainpfmt)
: inpfmt(ainpfmt), input(ainput), n(0), m(0), eps(aeps), eps_a2(0), eps_b2(0), miu(amiu)
, pstart(nullptr), edges(nullptr), reverse(nullptr), min_cn(nullptr)
, pa(nullptr), rank(nullptr), cid(nullptr)
, degree(nullptr), similar_degree(nullptr), effective_degree(nullptr)
,noncore_cluster(), ieids()
{
	if(aeps < 0 || aeps > 1)
		throw invalid_argument("aeps should E [0, 1]");

	if(aeps == 0 || aeps == 1) {  // Note: strict == is fine here
		eps_b2 = 1;
		eps_a2 = aeps;
		return;
	}

	eps_b2 = 1e4;
	eps_a2 = round(aeps * eps_b2);  // Note: uint_32 is limited with 2^32, floor( log10(2^16) ) < 5  ==> 1e4 is max possible accuracy
	eps_a2 *= eps_a2;
	eps_b2 *= eps_b2;
}

Graph::~Graph() {
	if(pstart) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(edges) {
		delete[] edges;
		edges = nullptr;
	}
	if(reverse) {
		delete[] reverse;
		reverse = nullptr;
	}
	if(min_cn) {
		delete[] min_cn;
		min_cn = nullptr;
	}
	if(cid) {
		delete[] cid;
		cid = nullptr;
	}
	if(degree) {
		delete[] degree;
		degree = nullptr;
	}
	if(effective_degree) {
		delete[] effective_degree;
		effective_degree = nullptr;
	}
	if(similar_degree) {
		delete[] similar_degree;
		similar_degree = nullptr;
	}
	if(pa) {
		delete[] pa;
		pa = nullptr;
	}
	if(rank) {
		delete[] rank;
		rank = nullptr;
	}
}

void Graph::load() {
	if(inpfmt == format_arg_BIN)
		loadBinary();
	else loadNSL();

	// Verify loaded edges
	for(Id i = 0;i < n;i ++) {
		for(Id j = pstart[i];j < pstart[i+1];j ++) {
			if(edges[j] == i)
				throw domain_error("Self loops are not allowed\n");
			if(j > pstart[i] && edges[j] <= edges[j-1])
				throw domain_error("Edges not sorted in increasing id order!\nThe program may not run properly!\n");
		}
	}
}

void Graph::loadBinary() {
	string  fname = input + "/b_degree.bin";
	FILE *f = fopen(fname.c_str(), "rb");
	if(!f) {
		perror(fname.insert(0, "Error, on opening ").c_str());
		throw ios_base::failure(strerror(errno));
	}

	int tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != (int)sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	// printf("\tn = %u; m = %u\n", n, m/2);

	if(!degree)
		degree = new Degree[n];
	fread(degree, sizeof(unsigned), n, f);  // ATTENTION: degrees are specified as "unsigned int" in the input file

#ifdef _DEBUG_
	long long sum = 0;
	for(Id i = 0;i < n;i ++) sum += degree[i];
	assert(sum == m && "WA input graph");
#endif

	fclose(f);

	fname = input + "/b_adj.bin";
	f = fopen(fname.c_str(), "rb");
	if(!f) {
		perror(fname.insert(0, "Error, on opening ").c_str());
		throw ios_base::failure(strerror(errno));
	}

	if(!pstart)
		pstart = new Id[n+1];
	if(!edges)
		edges = new Id[m];
	if(!reverse)
		reverse = new Id[m];
	memset(reverse, 0, sizeof(Id)*m);
	if(!min_cn)
		min_cn = new int[m];
	memset(min_cn, 0, sizeof(int)*m);

	Degree *buf = new Degree[n];

	pstart[0] = 0;
	for(Id i = 0; i < n; i++) {
		//printf("%d %d\n", i, degree[i]);
		if(degree[i] > 0)
			fread(buf, sizeof(Degree), degree[i], f);

		for(Degree j = 0; j < degree[i]; j++)
			edges[pstart[i] + j] = buf[j];

		pstart[i+1] = pstart[i] + degree[i];

		++ degree[i];
	}

	delete[] buf;

	fclose(f);
}

void Graph::loadNSL() {
	assert(inpfmt != format_arg_BIN && "Inappropriate file format");

	if(inpfmt == format__NULL) {
		// Use extension to identify the file format
		size_t iext = input.rfind('.');
		if(iext != string::npos) {
			if(!strcmp(&input.c_str()[iext+1], "nse"))
				inpfmt = format_arg_NSE;
			else if(!strcmp(&input.c_str()[iext+1], "nsa"))
				inpfmt = format_arg_NSA;
		}
	}

	ifstream  finp;
	finp.exceptions(ifstream::badbit);  //  | ifstream::failbit  - raises exception on EOF
	finp.open(input);
	if(!finp.is_open()) {
		perror(("Error opening the file " + input).c_str());
		throw std::ios_base::failure(strerror(errno));
	}

	// Intermediate data structures to fill internal data structures
	unordered_map<Id, Id>  eiids;  // Map from external to internal ids
	vector<set<Id>>  nodes;  // Nodes with arcs specified with internal ids

	// Parse NSE/A file
	string line;

	// Parse the header
	// [Nodes: <nodes_num>[,]	<Links>: <links_num>[,] [Weighted: {0, 1}]]
	// Note: the comma is either always present as a delimiter or always absent
	while(getline(finp, line)) {
		// Skip empty lines
		if(line.empty())
			continue;
		// Consider only subsequent comments
		if(line[0] != '#')
			break;

		//// 1. Replace the staring comment mark '#' with space to allow "#clusters:"
		//line[0] = ' ';
		// 2. Replace ':' with space to allow "Clusters:<clsnum>"
		for(size_t pos = 0; pos != string::npos; pos = line.find(':', pos + 1))
			line[pos] = ' ';

		// Parse nodes num
        char *tok = strtok(const_cast<char*>(line.data()) + 1, " \t");  // Note: +1 to omit the starting '#'
        if(!tok)
			continue;
		tolower(tok);
		if(strcmp(tok, "nodes"))
			continue;
		// Read nodes num
		tok = strtok(nullptr, " \t");
		if(tok) {
			// Note: optional trailing ',' is allowed here
			n = strtoul(tok, nullptr, 10);
			// Read the number of links
			tok = strtok(nullptr, " \t");
			if(tok) {
				tolower(tok);
				if(inpfmt != format_arg_NSE && !strcmp(tok, "arcs"))
					inpfmt = format_arg_NSA;
				else if(inpfmt != format_arg_NSA && !strcmp(tok, "edges"))
					inpfmt = format_arg_NSE;
				else throw domain_error(string("The file format (").append(tok)
						.append(") is either inconsistent with the expected one(")
						+= to_string(inpfmt).append(")\n"));
				tok = strtok(nullptr, " \t");
				if(tok) {
					// Note: optional trailing ',' is allowed here
					m = strtoul(tok, nullptr, 10);
					// Read Weighted flag
					tok = strtok(nullptr, " \t");
					if(tok && (tolower(tok), !strcmp(tok, "weighted"))
					&& (tok = strtok(nullptr, " \t")) && strtoul(tok, nullptr, 10) != 0)
						fputs("WARNING, the network is weighted and this algorithm does not support weights"
							", so the weights are omitted.\n", stderr);
				}
			}
		}
	}
	assert((inpfmt == format_arg_NSE || inpfmt == format_arg_NSA) && "Unexpected inpfmt");

	// Preallocate containers if possible
	if(n) {
		nodes.reserve(n);
		eiids.reserve(n);
		ieids.reserve(n);
	}

	// Parse the body
	// Note: the processing is started from the read line
    size_t iline = 0;  // Payload line index, links counter
    do {
		// Skip empty lines and comments
		if(line.empty() || line[0] == '#')
			continue;

        char *tok = strtok(const_cast<char*>(line.data()), " \t");
        if(!tok)
			continue;
		Id sid = strtoul(tok, nullptr, 10);  // External source id
        tok = strtok(nullptr, " \t");
        if(!tok)
			throw domain_error(string(line).insert(0
				, "Destination link id is not specified in this line: ").c_str());
		Id did = strtoul(tok, nullptr, 10);  // External destination id
		// Make the mappings and fill the nodes
		auto ies = eiids.find(sid);  // Index of the external src id
		if(ies == eiids.end()) {
			ies = eiids.emplace(sid, eiids.size()).first;
			ieids.emplace(ies->second, sid);
			nodes.emplace_back();
		}
		auto ied = eiids.find(did);  // Index of the external dst id
		if(ied == eiids.end()) {
			ied = eiids.emplace(did, eiids.size()).first;
			ieids.emplace(ied->second, did);
			nodes.emplace_back();
		}
		nodes[ies->second].insert(ied->second);
		//fprintf(stderr, "+ arc: %u %u  [%u %u]\n", ies->second, ied->second, sid, did);

		// Insert back arc in case of edges
        if(inpfmt == format_arg_NSE)
			nodes[ied->second].insert(ies->second);
		++iline;
    } while(getline(finp, line));
    assert(eiids.size() == ieids.size() && nodes.size() == ieids.size()
		&& "Node mappings are not synchronized");

	// Initialize internal data structures using nodes
	assert((!n || n == nodes.size()) && "Nodes size validation failed");
	if(n != nodes.size())
		n = nodes.size();

	const size_t  arcsnum = iline * (1 + (inpfmt == format_arg_NSE));
	assert((!m || m == arcsnum) && "Arcs size validation failed");
	if(m != arcsnum)
		m = arcsnum;

	if(reverse)
		delete[] reverse;
	reverse = new Id[m];
	memset(reverse, 0, sizeof(Id)*m);
	if(min_cn)
		delete[] min_cn;
	min_cn = new int[m];
	memset(min_cn, 0, sizeof(int) * m);

	if(edges)
		delete[] edges;
	edges = new Id[m];
	if(pstart)
		delete[] pstart;
	pstart = new Id[n+1];
	if(degree)
		delete[] degree;
	degree = new Degree[n];

	size_t in = 0;
	size_t ie = 0;
	pstart[0] = 0;
	for(auto& nd: nodes) {
		//fprintf(stderr, "nodes[%lu] size: %lu\n", in, nd.size());
		degree[in++] = nd.size();
		for(auto did: nd)
			edges[ie++] = did;
		pstart[in] = ie;
		assert(ie == pstart[in-1] + degree[in-1] && "pstart item validation failed");
	}
	if(ie != m) {
		fprintf(stderr, "ie: %lu, m: %u\n", ie, m);
		assert(0 && "Arcs number validation failed");
	}
}

Id Graph::binary_search(const Id *edges, Id b, Id e, Id val) {
	assert(b <= e && "Invalid indices");
	if(b == e || edges[e-1] < val)
		return e;
	--e;
	while(b < e) {
		Id mid = b + (e-b)/2;
		if(edges[mid] >= val) e = mid;
		else b = mid+1;
	}
#ifdef _DEBUG_
	assert(edges[e] >= val && "Ordering validation failed");
#endif
	return e;
}

void Graph::cluster_noncore_vertices(int mu) {
	if(!cid)
		cid = new Id[n];
	for(Id i = 0; i < n; i++)
		cid[i] = n;

	for(Id i = 0; i < n; i++)
		if(similar_degree[i] >= mu) {
			Id x = find_root(i);
			if(i < cid[x])
				cid[x] = i;
		}

	noncore_cluster.clear();
	noncore_cluster.reserve(n);
	for(Id i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		for(Id j = pstart[i];j < pstart[i+1];j ++) if(similar_degree[edges[j]] < mu) {
			if(min_cn[j] >= 0) {
				min_cn[j] = similar_check_OP(i, j);
				if(reverse[reverse[j]] != j) printf("WA cluster_noncore\n");
				min_cn[reverse[j]] = min_cn[j];
				if(min_cn[j] == -1) {
					++ similar_degree[i];
					++ similar_degree[edges[j]];
				}
			}

			if(min_cn[j] == -1)
				noncore_cluster.emplace_back(cid[pa[i]], edges[j]);
		}
	}
}

void Graph::output(const char *outfile, bool lgcfmt) {
	if(lgcfmt)
		saveLegacy(outfile);
	else saveCNL(outfile);
}

void Graph::saveLegacy(const char *outfile) {
	string out_name = outfile ? string(outfile)
		: input + "/result-" + to_string(eps) + "-" + to_string(miu) + ".txt";
	FILE *fout = fopen(out_name.c_str(), "w");
	if(!fout) {
		perror(out_name.insert(0, "Error, on opening ").c_str());
		throw ios_base::failure(strerror(errno));
	}

	fputs("# vertex_id cluster_id\n", fout);
	if(!noncore_cluster.empty())
		fputs("# Core clusters ---\n", fout);
	// Map internal ids to the original external ids if required
	for(Id i = 0;i < n;i ++) if(similar_degree[i] >= miu)
		fprintf(fout, "%d %d\n", !ieids.empty() ? ieids[i] : i, cid[pa[i]]);


	if(!noncore_cluster.empty()) {
		fputs("\n# Non-core clusters ---\n", fout);
		sort(noncore_cluster.begin(), noncore_cluster.end());
		noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
		// Map internal ids to the original external ids if required
		for(Id i = 0;i < noncore_cluster.size();i ++)
			fprintf(fout, "%d %d\n", !ieids.empty()
				? ieids[noncore_cluster[i].second] : noncore_cluster[i].second
				, noncore_cluster[i].first);
	}

	fclose(fout);
}

void Graph::saveCNL(const char *outfile) {
	FILE *fout = fopen(outfile, "w");

	if(!fout) {
		perror(string("Error, on opening ").append(outfile).c_str());
		throw ios_base::failure(strerror(errno));
	}

	typedef vector<Id>  NodeIds;
	unordered_map<Id, NodeIds>  cls;
	cls.reserve(sqrt(n));  // Usually the number of clusters <= sqrt(nodes_num)
	for(Id i = 0; i < n; ++i)
		if(similar_degree[i] >= miu) {
			const auto key = cid[pa[i]];
			auto icl = cls.find(key);
			if(icl == cls.end());
				icl = cls.emplace(key, NodeIds()).first;
			// Map internal ids to the original external ids if required
			icl->second.push_back(!ieids.empty() ? ieids[i] : i);
		}

	// Save the optional Header
	fprintf(fout, "# Clusters: %lu, Nodes: %u, Fuzzy: 0\n"
		"# Each line contains member nodes of each resulting cluster\n", cls.size(), n);
	// Save the Body
	if(!noncore_cluster.empty())
		fputs("# Core clusters ---\n", fout);
	for(const auto& cl: cls) {
		for(const auto& nid: cl.second)
			fprintf(fout, "%u ", nid);
		fputs("\n", fout);
	}

	if(!noncore_cluster.empty()) {
		sort(noncore_cluster.begin(), noncore_cluster.end());
		noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
		cls.clear();
		for(Id i = 0; i < noncore_cluster.size(); ++i) {
			const auto key = noncore_cluster[i].first;
			auto icl = cls.find(key);
			if(icl == cls.end());
				icl = cls.emplace(key, NodeIds()).first;
			// Map internal ids to the original external ids if required
			icl->second.push_back(!ieids.empty() ? ieids[noncore_cluster[i].second] : noncore_cluster[i].second);
		}
		fputs("\n# Non-core clusters ---\n", fout);
		// Save the clusters
		for(const auto& cl: cls) {
			for(const auto& nid: cl.second)
				fprintf(fout, "%u ", nid);
			fputs("\n", fout);
		}
	}

	fclose(fout);
}

void Graph::pSCAN() {
	if(!similar_degree)
		similar_degree = new Degree[n];
	memset(similar_degree, 0, sizeof(Degree)*n);

	if(!effective_degree)
		effective_degree = new Degree[n];
	for(Id i = 0;i < n;i ++)
		effective_degree[i] = degree[i]-1;

	if(!pa)
		pa = new Id[n];
	if(!rank)
		rank = new int[n];
	for(Id i = 0;i < n;i ++) {
		pa[i] = i;
		rank[i] = 0;
	}

#ifdef __unix__
	struct timeval start;
	gettimeofday(&start, nullptr);
#else
	int start = clock();
#endif

	Id *edge_buf = new Id[n];
	int *cores = new int[n];
	int cores_n = 0;

	prune_and_cross_link(miu, cores, cores_n);
	//printf("\t*** Finished prune and cross link!\n");

#ifdef __unix__
	struct timeval end1;
	gettimeofday(&end1, nullptr);
#else
	int end1 = clock();
#endif

	int *bin_head = new int[n];
	int *bin_next = new int[n];
	for(Id i = 0;i < n;i ++)
		bin_head[i] = -1;

	int max_ed = 0;
	for(Id i = 0;i < n;i ++)
		if(effective_degree[i] >= miu) {
			Degree ed = effective_degree[i];
			if(ed > max_ed)
				max_ed = ed;
			bin_next[i] = bin_head[ed];
			bin_head[ed] = i;
		}

	while(true) {
		int u = -1;
		if(cores_n)
			u = cores[-- cores_n];
		else {
			while(max_ed >= miu&&u == -1) {
				for(int x = bin_head[max_ed];x != -1;) {
					int tmp = bin_next[x];
					Degree ed = effective_degree[x];
					if(ed == max_ed) {
						u = x;
						bin_head[max_ed] = bin_next[x];
						break;
					}
					else if(ed >= miu) {
						bin_next[x] = bin_head[ed];
						bin_head[ed] = x;
					}
					x = tmp;
				}
				if(u == -1) {
					bin_head[max_ed] = -1;
					-- max_ed;
				}
			}
		}
		if(u == -1)
			break;

		Id edge_buf_n = 0;
		for(Id j = pstart[u];j < pstart[u+1];j ++) {
			if(min_cn[j] == -2)
				continue;

			if(similar_degree[u] < miu||find_root(u) != find_root(edges[j]))
				edge_buf[edge_buf_n ++] = j;
		}

		Id i = 0;
		while(similar_degree[u] < miu&&effective_degree[u] >= miu&&i < edge_buf_n) {
			Id idx = edge_buf[i];
			if(min_cn[idx] != -1) {
#ifdef _DEBUG_
				assert(min_cn[idx] && "min_cn item is invalid");
#endif
				Id v = edges[idx];

				min_cn[idx] = min_cn[reverse[idx]] = similar_check_OP(u, idx);

				if(min_cn[idx] == -1) ++ similar_degree[u];
				else -- effective_degree[u];

				if(effective_degree[v] >= 0) {
					if(min_cn[idx] == -1) {
						++ similar_degree[v];

						if(similar_degree[v] == miu) cores[cores_n ++] = v;
					}
					else -- effective_degree[v];
				}
			}

			++ i;
		}

		effective_degree[u] = -1;

		if(similar_degree[u] < miu) continue;

		for(Id j = 0;j < edge_buf_n;j ++) {
			Id idx = edge_buf[j];
			if(min_cn[idx] == -1&&similar_degree[edges[idx]] >= miu)
				my_union(u, edges[idx]);
		}

		while(i < edge_buf_n) {
			Id idx = edge_buf[i];
			Id v = edges[idx];
			if(min_cn[idx] < 0||similar_degree[v] < miu||find_root(u) == find_root(v)) {
				++ i;
				continue;
			}

			min_cn[idx] = min_cn[reverse[idx]] = similar_check_OP(u, idx);

			if(effective_degree[v] >= 0) {
				if(min_cn[idx] == -1) {
					++ similar_degree[v];

					if(similar_degree[v] == miu) cores[cores_n ++] = v;
				}
				else -- effective_degree[v];
			}

			if(min_cn[idx] == -1) my_union(u,v);

			++ i;
		}
		//printf(")\n");
	}
	//printf("\t*** Finished clustering core vertices!\n");

	delete[] edge_buf; edge_buf = nullptr;
	delete[] cores; cores = nullptr;
	delete[] bin_head; bin_head = nullptr;
	delete[] bin_next; bin_next = nullptr;

#ifdef __unix__
	struct timeval end;
	gettimeofday(&end, nullptr);

//	long long seconds, useconds;
//	seconds = end.tv_sec - end1.tv_sec;
//	useconds = end.tv_usec - end1.tv_usec;
//
//	printf("Prune time: %lld\nRefine time: %lld\n", mtime1, seconds*1000000 + useconds);
#else
	int end = clock();

	printf("Prune time: %d\nSort time: %d\nRefine time: %d\n", end1-start,end2-2end1,end-end2);
#endif

	cluster_noncore_vertices(miu);
}

int Graph::check_common_neighbor(Id u, Id v, int c) {
	int cn = 2;

	if(degree[u] > degree[v])
		swap(u,v);

	Degree du = degree[u]+1, dv = degree[v]+1;
	Id i = pstart[u], j = pstart[v];
	while(i < pstart[u+1]&&j < pstart[v+1]&&cn < c&&du >= c&&dv >= c) {
		if(edges[i] < edges[j]) {
			-- du;
			++ i;
		}
		else if(edges[i] > edges[j]) {
			-- dv;
			++ j;
		}
		else {
			++ cn;
			++ i;
			++ j;
		}
	}

	if(cn >= c) return -1;
	return -2;
}

int Graph::similar_check_OP(Id u, Id idx) {
	Id v = edges[idx];

#ifdef _DEBUG_
	assert(min_cn[idx] >= 0 && "min_cn item validation failed");
#endif

	if(min_cn[idx] == 0) {
		Degree du = degree[u], dv = degree[v];
		int c = compute_common_neighbor_lowerbound(du,dv);

#ifdef _DEBUG_
		if(du < c||dv < c) return -2;
#endif

		if(c <= 2) return -1;

		min_cn[idx] = min_cn[reverse[idx]] = c;
	}

	return check_common_neighbor(u, v, min_cn[idx]);
}

int Graph::compute_common_neighbor_lowerbound(Id du,Id dv) {
	int c = (int)(sqrtl((((long double)du) * ((long double)dv)*eps_a2)/eps_b2));

#ifdef _DEBUG_
	assert(!(((long long)du)*dv*eps_a2 < 0 || ((long long)c)*c*eps_b2 < 0) && "Overflow");
#endif

	if(((long long)c)*((long long)c)*eps_b2 < ((long long)du)*((long long)dv)*eps_a2)
		++ c;

#ifdef _DEBUG_
	assert(((long long)c)*((long long)c)*eps_b2 >= ((long long)du)*((long long)dv)*eps_a2
		&& "Common neighbor validation failed");
#endif
	return c;
}

void Graph::prune_and_cross_link(int miu, int *cores, int &cores_e) {
	for(Id i = 0;i < n;i ++) { //must be iterating from 0 to n-1
		for(Id j = pstart[i];j < pstart[i+1];j ++) {
			if(edges[j] < i) {
				if(min_cn[j] == 0) min_cn[j] = -2;
				continue; //this edge has already been checked
			}

			Id v = edges[j];
			Degree a = degree[i], b = degree[v];
			if(a > b) swap(a, b);

			if(((long long)a)*eps_b2 < ((long long)b)*eps_a2) {
				min_cn[j] = -2;

				-- effective_degree[i];
				-- effective_degree[v];
			}
			else {
				int c = compute_common_neighbor_lowerbound(a, b);

#ifdef _DEBUG_
				assert(!(a < c || b < c) && "Invalid values");
#endif

				if(c <= 2) {
					min_cn[j] = -1;

					++ similar_degree[i];
					++ similar_degree[v];

					if(similar_degree[i] == miu) cores[cores_e ++] = i;
					if(similar_degree[v] == miu) cores[cores_e ++] = v;
				}
				else min_cn[j] = c;
			}

			if(min_cn[j] != -2) {
			//else {
				Id r_id = binary_search(edges, pstart[v], pstart[v+1], i);
				reverse[j] = r_id;
				reverse[r_id] = j;

				min_cn[r_id] = min_cn[j];
			}
		}
	}
}

Id Graph::find_root(Id u) {
	Id x = u;
	while(pa[x] != x)
		x = pa[x];

	while(pa[u] != x) {
		Id tmp = pa[u];
		pa[u] = x;
		u = tmp;
	}

	return x;
}

void Graph::my_union(Id u, Id v) {
	Id ru = find_root(u);
	Id rv = find_root(v);

	if(ru == rv) return ;

	if(rank[ru] < rank[rv])
		pa[ru] = rv;
	else if(rank[ru] > rank[rv])
		pa[rv] = ru;
	else {
		pa[rv] = ru;
		++ rank[ru];
	}
}
