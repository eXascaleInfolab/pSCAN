#include "Utility.h"
#include "Graph.h"

Graph::Graph(const char *_dir)
: dir(_dir ? _dir : ""), n(0), m(0), eps_a2(0), eps_b2(0), miu(0)
, pstart(nullptr), edges(nullptr), reverse(nullptr), min_cn(nullptr)
, pa(nullptr), rank(nullptr), cid(nullptr)
, degree(nullptr), similar_degree(nullptr), effective_degree(nullptr)
,noncore_cluster()  {}

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

void Graph::read_graph(const char *filename) {
	FILE *f = open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != (int)sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	// printf("\tn = %u; m = %u\n", n, m/2);

	degree = new Degree[n];
	fread(degree, sizeof(unsigned), n, f);  // ATTENTION: degrees are specified as "unsigned int" in the input file

#ifdef _DEBUG_
	long long sum = 0;
	for(Id i = 0;i < n;i ++) sum += degree[i];
	if(sum != m) printf("WA input graph\n");
#endif

	fclose(f);

	f = open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(!pstart)
		pstart = new Id[n+1];
	if(!edges)
		edges = new Id[m];
	if(!reverse)
		reverse = new Id[m];
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

	for(Id i = 0;i < n;i ++) {
		for(Id j = pstart[i];j < pstart[i+1];j ++) {
			if(edges[j] == i) {
				printf("Self loop\n");
				//exit(1);
			}
			if(j > pstart[i]&&edges[j] <= edges[j-1]) {
				printf("Edges not sorted in increasing id order!\nThe program may not run properly!\n");
				//exit(1);
			}
		}
	}
}

Id Graph::binary_search(const Id *edges, Id b, Id e, Id val) {
#ifdef _DEBUG_
	if(e < b)
		printf("??? WA1 in binary_search\n");
#endif
	if(b == e || edges[e-1] < val)
		return e;
	--e;
	while(b < e) {
		Id mid = b + (e-b)/2;
		if(edges[mid] >= val) e = mid;
		else b = mid+1;
	}
#ifdef _DEBUG_
	if(edges[e] < val)
		printf("??? WA2 in binary_search\n");
#endif
	return e;
}

void Graph::cluster_noncore_vertices(int eps_a2, int eps_b2, int mu) {
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
				min_cn[j] = similar_check_OP(i, j, eps_a2, eps_b2);
				if(reverse[reverse[j]] != j) printf("WA cluster_noncore\n");
				min_cn[reverse[j]] = min_cn[j];
				if(min_cn[j] == -1) {
					++ similar_degree[i];
					++ similar_degree[edges[j]];
				}
			}

			if(min_cn[j] == -1)
				noncore_cluster.push_back(make_pair(cid[pa[i]], edges[j]));
		}
	}
}

void Graph::output(const char *eps_s, const char *miu) {
	printf("\t*** Start write result into disk!\n");

	string out_name = dir + "/result-" + string(eps_s) + "-" + string(miu) + ".txt";
	FILE *fout = open_file(out_name.c_str(), "w");

	fprintf(fout, "c/n vertex_id cluster_id\n");

	int mu = atoi(miu);
	for(Id i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		fprintf(fout, "c %d %d\n", i, cid[pa[i]]);
	}

	sort(noncore_cluster.begin(), noncore_cluster.end());
	noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
	for(Id i = 0;i < noncore_cluster.size();i ++) {
		fprintf(fout, "n %d %d\n", noncore_cluster[i].second, noncore_cluster[i].first);
	}

	fclose(fout);
}

void Graph::pSCAN(const char *eps_s, int _miu) {
	get_eps(eps_s);
	miu = _miu;

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

	prune_and_cross_link(eps_a2, eps_b2, miu, cores, cores_n);
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
				if(min_cn[idx] == 0) printf("WA min_cn!\n");
#endif
				Id v = edges[idx];

				min_cn[idx] = min_cn[reverse[idx]] = similar_check_OP(u, idx, eps_a2, eps_b2);

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

			min_cn[idx] = min_cn[reverse[idx]] = similar_check_OP(u, idx, eps_a2, eps_b2);

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

	cluster_noncore_vertices(eps_a2, eps_b2, miu);
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

int Graph::similar_check_OP(Id u, Id idx, int eps_a2, int eps_b2) {
	Id v = edges[idx];

#ifdef _DEBUG_
	if(min_cn[idx] == -1||min_cn[idx] == -2) printf("??? WA in similar_check\n");
#endif

	if(min_cn[idx] == 0) {
		Degree du = degree[u], dv = degree[v];
		int c = compute_common_neighbor_lowerbound(du,dv,eps_a2,eps_b2);

#ifdef _DEBUG_
		if(du < c||dv < c) return -2;
#endif

		if(c <= 2) return -1;

		min_cn[idx] = min_cn[reverse[idx]] = c;
	}

	return check_common_neighbor(u, v, min_cn[idx]);
}

int Graph::compute_common_neighbor_lowerbound(Id du,Id dv,int eps_a2,int eps_b2) {
	int c = (int)(sqrtl((((long double)du)*((long double)dv)*eps_a2)/eps_b2));

#ifdef _DEBUG_
	if(((long long)du)*dv*eps_a2 < 0||((long long)c)*c*eps_b2 < 0)
		printf("??? Overflow in similar_check\n");
#endif

	if(((long long)c)*((long long)c)*eps_b2 < ((long long)du)*((long long)dv)*eps_a2) ++ c;

#ifdef _DEBUG_
	if(((long long)c)*((long long)c)*eps_b2 < ((long long)du)*((long long)dv)*eps_a2)
		printf("??? Wrong common neigbor computation in similar_check\n");
#endif
	return c;
}

void Graph::prune_and_cross_link(int eps_a2, int eps_b2, int miu, int *cores, int &cores_e) {
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
				int c = compute_common_neighbor_lowerbound(a, b, eps_a2, eps_b2);

#ifdef _DEBUG_
				if(a < c||b < c) printf("!!! HHH\n");
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

void Graph::get_eps(const char *eps_s) {
	int i = 0, eps_a = 0, eps_b = 1;
	while(eps_s[i] != '\0'&&eps_s[i] != '.') {
		eps_a = eps_a*10 + (eps_s[i]-'0');
		++ i;
	}

	if(eps_s[i] == '.') {
		++ i;
		while(eps_s[i] != '\0') {
			eps_a = eps_a*10 + (eps_s[i]-'0');
			eps_b *= 10;
			++ i;
		}
	}

	if(eps_a > eps_b||eps_b > 100||eps_a <= 0) {
		printf("??? Wrong eps format: %d/%d, %s\n", eps_a, eps_b, eps_s);
		exit(1);
	}

	eps_a2 = eps_a * eps_a;
	eps_b2 = eps_b * eps_b;
}
