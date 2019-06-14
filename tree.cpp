#include "tree.hpp"
#include "support_wt.hpp"

ll* bp;

ll* cum_bp;

ll* fcount_unary;

ll* cum_fcount;

ll* no_zero_node;

ll* leaf_lead_char;

uint8_t* pbwt;

int bp_size , no_zero_leaf_size , fcount_size;

wavelet_tree* wt_zero_dep;
wavelet_tree* wt_bwt;

int n, sigma_p;
int leaf_indx, rnk, no_of_bits , no_of_nodes;
int* prev_encoded;

int f(int b, int j){
	if(b >= n)return b;
	if(b <= j-1)return b;
	else return 0;
}

void P_suffixtree::encode(string s, bool* param){
	int last[256];
	memset(last,-1,sizeof(last));
	if(param == NULL){
		for(int i = 0;i < n-1;i++){
			if(last[s[i]] == -1)encod[i] = 0;
			else encod[i] = i - last[s[i]];
			last[s[i]] = i;
		}
		encod[n-1] = n;
	}
	else {
		map<char,int> tosort;
		for(int i = 0;i < n-1;i++){
			if(param[s[i]] == 0){
				tosort[s[i]] = 0;
			}
		}
		int j = 0;
		for(auto it = tosort.begin();it != tosort.end();it++,j++){
			it->second = n + j;
		}
		for(int i = 0;i < n-1;i++){
			if(param[s[i]] == 0){
				encod[i] = tosort[s[i]];
				continue;
			}
			if(last[s[i]] == -1)encod[i] = 0;
			else encod[i] = i - last[s[i]];
			last[s[i]] = i;
		}
		encod[n-1] = n+j;
	}
	
	// for(int i = 0;i < n;i++)cout<<s[i]<<" "<<encod[i]<<endl;

}

void P_suffixtree::construct_pbwt(string s, bool* param){
	pbwt = new uint8_t[n];
	int last[256];
	vector<int> v(prev_encoded, prev_encoded+2*n);
	wavelet_tree suprt(v);
	memset(last,-1,sizeof(last));
	for(int i = 2*n-1;i >= n;i--){
		last[s[i%n]] = i;
	}
	for(int i = n-1;i >= 0;i--){
		if(encod[i] >= n){
			pbwt[rank[(i+1)%n]] = encod[i] - n + sigma_p + 1;
		}
		else {
			pbwt[rank[(i+1)%n]] = suprt.twod(i+2, last[s[i]]+1, 0, i+1);
		}
		last[s[i]] = i;
	}
	// for(int i = 0;i < n;i++){
	// 	cout<<pbwt[i]<<" ";
	// }
	// cout<<endl;
}

// int P_suffixtree::plf(int i){
// 	if(pbwt[i] > sigma_p){
// 		return 1 + wt_bwt->twod(1,n,1,pbwt[i]-1) + wt_bwt->twod(1,i,pbwt[i],pbwt[i]);
// 	}
// 	else{
// 		node* z = zeronode(leaves[i]);   // i-th leaf can be found using bp encoding instead of leaves[i]
//		// node* v = z->parent;
// 		int l = z->lmost_leaf, r = z->rmost_leaf;  // l and r are assumed to be 1-based in the next line. l++,r++ if not
// 		int N1 = fsum(z);               // actually fsum(other node). find that other node    level ancestor of (i-1)th leaf
// 		int N2 = wt_bwt->twod(l,r,pbwt[i]+1,sigma_p) + wt_bwt->twod(l,i+1,pbwt[i],pbwt[i]);
// 		//int N3 = 0;
// 		int N4 = 0;
// 		if(leaf_lead_char[i]){       // leaf_lead_char is constructed in update_fcount function
// 			z = v->child(v->pcount); // how to find pcount ?
// 			N4 = wt_bwt->twod(r+1,z->rmost_leaf,pbwt[i],sigma_p);
// 		}
// 		return N1+N2+N4;
// 	}
// }

P_suffixtree::P_suffixtree(string s, bool* param, int sigma){
	sigma_p = sigma;
	s += '$';
	n = s.size();
	cout<<n<<"\n";
	root = new node;
	no_of_nodes++;
	encod = new int[n];
	sufix_aray = new int[n];
	rank = new int[n];
	leaves.resize(n);
	encode(s,param);
	oldhd = root;
	oldchild = NULL;

	int last[256];
	prev_encoded = new int[2*n];
	memset(last,-1,sizeof(last));
	for(int i = 0;i < 2*n;i++){
		if(param[s[i%n]] == 0){
			prev_encoded[i] = n+1;
			continue;
		}
		prev_encoded[i] = last[s[i%n]]+1;
		last[s[i%n]] = i;
	}
	vector<int> v(prev_encoded, prev_encoded+n);
	wt_zero_dep = new wavelet_tree(v);

	for(int i = 0;i < n;i++){
		insert(i);
	}
	delete wt_zero_dep;
	// suffix tree is built

	no_of_bits = no_of_nodes;
	int m = (2*no_of_nodes + 63)/64;
	bp = new ll[m];
	cum_bp = new ll[m];
	bp_size = 2*no_of_nodes;
	for(int i = 0;i < m;i++)cum_bp[i] = bp[i] = 0;
	for(int i = 1;i < m;i++)cum_bp[i] += cum_bp[i - 1];
	suffix_dfs(root);

	for(int i = 0;i < n;i++)rank[sufix_aray[i]] = i;

	m = (n+63)/64;
	no_zero_node = new ll[m];
	leaf_lead_char = new ll[m];
	no_zero_leaf_size = m;
	
	for(int i = 0;i < m;i++)no_zero_node[i] = leaf_lead_char[i] = 0;

	int nex[256];
	memset(nex, -1, sizeof(nex));
	for(int i = n-1;i >= 0;i--){
		if(nex[s[i]] != -1 and param[s[i]]){
			// cout<<i<<" "<<nex[s[i]]<<endl;
			update_fcount(i+1, nex[s[i]] - i);
			no_of_bits++;
		}
		nex[s[i]] = i;
	}
	m = (no_of_bits + 63)/64;
	fcount_unary = new ll[m];
	cum_fcount = new ll[m];
	fcount_size = m;

	no_of_bits = 0;
	for(int i = 0;i < m;i++)fcount_unary[i] = cum_fcount[i] = 0;
	rnk = 0;
	dfs(root);
	for(int i = 1;i < m;i++)cum_fcount[i] += cum_fcount[i-1];
	construct_pbwt(s, param);
	cout<<"Suffix array: ";
	for(int i = 0;i < n;i++)cout<<sufix_aray[i]<<" ";cout<<endl;
	cout<<"Burrows wheeler transform: ";
	for(int i = 0;i < n;i++)cout<<(unsigned)pbwt[i]<<" ";cout<<endl;
	int r = 1;
	pre_dfs(root,r);
	wt_zero_dep = new wavelet_tree(zero_depth_nodes);      // wt_zero_dep can be used for zero_node queries

	memset(nex,-1,sizeof(nex));
	for(int i = n-1;i >= 0;i--){
		if(nex[s[i]] == -1 or param[s[i]] == 0){
			int j = rank[(i+1)%n];
			no_zero_node[j/64] |= 1ll<<(63 - j%64);     // 1 if zeronode doesn't exist
		}
		nex[s[i]] = i;
	}
	delete[] prev_encoded;
}

// node* P_suffixtree::zeronode(int i){
// 	node* leaf = leaves[i];
// 	node* tmp = leaf->parent;
// 	if(no_zero_node[i/64] & (1ll<<(63 - i%64)))return NULL;
// 	int pre_rank = wt_zero_dep->predecessor(tmp->pre_ordr_rnk, pbwt[i]+1);
// 	tmp = pre_ordr_rank_to_node(pre_rank);                                  // bp encoding supports this
// 	tmp = lca(tmp, leaf);
// 	tmp = child(tmp,leaf);
// 	return tmp;
// }

int P_suffixtree::no_of_zeros(int i, int j, int k){ // no of zeros in prev(i) from index j to index k 
	return wt_zero_dep->twod(j+1,k+1,0,i);
}

void P_suffixtree::insert(int i){
	node* tmp;
	node* temphd;
	node* newhd;
	int j, edgeln, zeros;
	if(oldhd == root){
		j = 1;
		edgeln = 0;
		zeros = 0;
		tmp = oldhd;
		while(j <= n-i){
			if(edgeln == 0){
				temphd = tmp;
				if(temphd->edge.find(f(encod[i+j-1], j)) != temphd->edge.end()){
					
					// cout<<encod[i + j - 1]<<" "<<i<<" "<<j<<endl;
					
					tmp = temphd->edge[f(encod[i+j-1], j)];
					edgeln = tmp->r - tmp->l + 1;
					zeros = 0;
				}
				else{
					
					// cout<<encod[i + j - 1]<<" "<<i<<" "<<j<<endl;
					
					// create a new child here
					newhd = temphd;
					
					// if(i+j-1 > n-1)cout<<"nc"<<" "<<i+j-1<<" "<<n - 1<<endl;
					
					newhd->edge[f(encod[i+j-1], j)] = new node(n-i, i+j-1, n-1, newhd);
					node* new_child_created = newhd->edge[f(encod[i+j-1], j)];
					new_child_created->zero_depth = newhd->zero_depth + no_of_zeros(i,i+j-1, n-1);
					no_of_nodes++;
					oldhd = newhd;
					oldchild = temphd;
					return;
				}
			}
			else {
				int ind = tmp->r + 1 - edgeln;
				
				// cout<<encod[i + j - 1]<<" "<<i<<" "<<j<<endl;
				
				if(f(encod[ind], tmp->pathln - tmp->r + ind) != f(encod[i+j-1], j)){
					
					// if(tmp->l > ind-1)cout<<"nn"<<tmp->l<<" "<<ind-1<<" "<<temphd->r<<endl;
					// if(i+j-1 > n-1)cout<<"nc"<<i+j-1<<n-1<<endl;
					
					// create a new node and child here
					temphd->edge[f(encod[tmp->l], tmp->pathln - tmp->r + tmp->l)] = new node(j-1, tmp->l, ind-1, temphd);
					newhd = temphd->edge[f(encod[tmp->l], tmp->pathln - tmp->r + tmp->l)];
					newhd->zero_depth = temphd->zero_depth + zeros;

					newhd->edge[f(encod[i+j-1], j)] = new node(n-i, i+j-1, n-1, newhd);
					node* new_child_created = newhd->edge[f(encod[i+j-1], j)];
					new_child_created->zero_depth = newhd->zero_depth + no_of_zeros(i, i+j-1, n-1);

					newhd->edge[f(encod[ind], tmp->pathln - tmp->r + ind)] = tmp;
					tmp->parent = newhd;
					tmp->l = ind;
					// call update function and then
					update(tmp, newhd, oldchild, i);
					oldhd = newhd;
					oldchild = tmp;
					no_of_nodes += 2;
					return;
				}
			}
			if(f(encod[i+j-1], j) == 0)zeros++;
			edgeln--;
			j++;
		}
	}
	else{
		if(oldhd->sl == NULL){
			tmp = oldhd->parent->sl;
			j = tmp->pathln;
			
			// if(i == 18)cout<<"ide"<<tmp->l<<" "<<tmp->r<<" "<<tmp->pathln<<" "<<oldhd->pathln<<" "<<oldhd->r<<" ";
			
			while(j < oldhd->pathln - 1){
				tmp = tmp->edge[f(encod[i+j], j+1)];
				j = tmp->pathln;
			}
			
			// cout<<j<<endl;
			
			oldhd->sl = tmp;
		}
		tmp = oldhd->sl;
		temphd = tmp->parent;
		edgeln = tmp->pathln + 1 - oldhd->pathln;
		
		// if(tmp -> l==11)cout<<i<<"sl"<<edgeln<<" "<<tmp->r<<" "<<tmp->pathln<<" "<<temphd->pathln<<" "<<temphd->l<<" "<<temphd->r<<" "<<oldhd->pathln<<endl;
		
		j = oldhd->pathln;
		zeros = no_of_zeros(i, i + temphd->pathln, i+j-2);
		while(j <= n-i){
			if(edgeln == 0){
				temphd = tmp;
				if(temphd->edge.find(f(encod[i+j-1], j)) != temphd->edge.end()){

					// cout<<encod[i + j - 1]<<" "<<i<<" "<<j<<endl;
					
					tmp = temphd->edge[f(encod[i+j-1], j)];
					edgeln = tmp->r - tmp->l + 1;
					
					// if(tmp -> l == 11)cout<<i<<" "<<edgeln<<"k"<<tmp->r<<endl;
					
					zeros = 0;
				}
				else{
					
					// cout<<encod[i + j - 1]<<" "<<i<<" "<<j<<endl;
					// if(i + j - 1 > n - 1)cout<<"nc"<<i+j-1<<" "<<n-1<<endl;

					// create a new child here
					newhd = temphd;
					newhd->edge[f(encod[i+j-1], j)] = new node(n-i, i+j-1, n-1, newhd);
					node* new_child_created = newhd->edge[f(encod[i+j-1], j)];
					new_child_created->zero_depth = newhd->zero_depth + no_of_zeros(i,i+j-1, n-1);
					no_of_nodes++;
					break;
				}
			}
			else {
				int ind = tmp->r + 1 - edgeln;

				// if(tmp->l == 11)cout<<i<<" "<<edgeln<<"ok"<<tmp->r<<" "<<tmp->r - tmp->pathln + 1<<endl;
				// cout<<f(encod[i+j-1], j)<<" "<<encod[i + j - 1]<<" "<<i<<" "<<j<<endl;

				if(f(encod[ind], tmp->pathln - tmp->r + ind) != f(encod[i+j-1], j)){
					
					// if(tmp->l == 11)cout<<"nn"<<i<<" "<<tmp->l<<" "<<tmp->r<<" "<<ind-1<<" "<<edgeln<<" "<<tmp->pathln<<" "<<temphd->l<<" "<<temphd->r<<" "<<temphd->pathln<<endl;
					// if(i+j-1 > n-1)cout<<"nc"<<i+j-1<<" "<<n-1<<endl;

					// create a new node and child here
					temphd->edge[f(encod[tmp->l], tmp->pathln - tmp->r + tmp->l)] = new node(j-1, tmp->l, ind-1, temphd);
					newhd = temphd->edge[f(encod[tmp->l], tmp->pathln - tmp->r + tmp->l)];
					newhd->zero_depth = temphd->zero_depth + zeros;

					newhd->edge[f(encod[i+j-1], j)] = new node(n-i, i+j-1, n-1, newhd);
					node* new_child_created = newhd->edge[f(encod[i+j-1], j)];
					new_child_created->zero_depth = newhd->zero_depth + no_of_zeros(i, i+j-1, n-1);

					newhd->edge[f(encod[ind], tmp->pathln - tmp->r + ind)] = tmp;
					tmp->parent = newhd;
					tmp->l = ind;
					// call update function
					update(tmp, newhd, oldchild, i);
					if(newhd->pathln < oldhd->sl->pathln)oldhd->sl = newhd;
					no_of_nodes += 2;
					break;
				}
			}
			if(f(encod[i+j-1], j) == 0)zeros++;
			edgeln--;
			j++;
		}
		if(oldhd->sl->min == NULL or oldhd->sl->min->pathln > oldhd->pathln)oldhd->sl->min = oldhd;
		oldhd = newhd;
		oldchild = tmp;
	}
}

void P_suffixtree::update(node* g, node* newhd, node* oldchild, int i){
	if(g->min == NULL or g->min->pathln > 1 + newhd->pathln)return;
	node *tmp;
	tmp = g->min;
	while(tmp->pathln <= 1 + newhd->pathln and tmp->sl == g){
		tmp->sl = newhd;
		tmp = tmp->edge[0];
	}
	if(tmp->sl == g and tmp->edge.size() != 0)g->min = tmp;
}

void P_suffixtree::update_fcount(int i, int first_occ){
	node* tmp = root;
	node* temphd;
	int j = 0;
	while(j < first_occ){
		temphd = tmp;
		tmp = temphd->edge[f(encod[i+j], j+1)];
		
		// if(tmp == NULL){cout<<i<<" "<<first_occ<<" "<<temphd->pathln<<" "<<encod[i+j]<<" "<<j+1<<" "<<n<<"oh no\n";return;}
		
		j = tmp->pathln;
	}
	if(first_occ == temphd->pathln + 1){
		temphd->fcount++;
		i = rank[i];
		leaf_lead_char[i/64] |= 1ll<<(63 - i%64);  // 1 if fi = path(v) + 1
	}
	else tmp->fcount++;
}

void P_suffixtree::suffix_dfs(node* p){
	bp[rnk/64] |= (1ll << (63 - rnk%64));
	cout<<"(";
	cum_bp[rnk/64]++;
	rnk++;
	for(auto it = p->edge.begin();it != p->edge.end();it++){
		
		// cout<<rnk<<" "<<it->first<<endl;
		
		suffix_dfs(it->second);
	}
	if(p->edge.size() == 0){
		sufix_aray[leaf_indx] = p->r + 1 - p->pathln;
		leaves[leaf_indx++] = p;
	}
	rnk++;
	cout<<")";
}

void P_suffixtree::dfs(node* p){
	for(auto it = p->edge.begin();it != p->edge.end();it++){
		dfs(it->second);
	}
	p->pst_ordr_rnk = ++rnk;
	
	// cout<<p->r + 1 - p->pathln<<" "<<p->pst_ordr_rnk<<" "<<p->fcount<<" "<<p->zero_depth<<endl;

	no_of_bits += p->fcount;
	fcount_unary[(no_of_bits)/64] |= (1ll<<(63 - no_of_bits%64));
	cum_fcount[(no_of_bits)/64]++;
	no_of_bits++;
}

void P_suffixtree::pre_dfs(node* p, int &r){
	p->pre_ordr_rnk = r++;
	if(p->edge.size() != 0)zero_depth_nodes.push_back(p->zero_depth);
	for(auto it = p->edge.begin();it != p->edge.end();it++){
		pre_dfs(it->second,r);
	}
}

P_suffixtree::~P_suffixtree(){
	destruct(root);
	delete[] encod;
	delete[] sufix_aray;
	delete[] rank;
	delete[] pbwt;
	leaves.clear();
	zero_depth_nodes.clear();
}

void P_suffixtree::destruct(node* head){
	for(auto it = head->edge.begin();it != head->edge.end();it++){
		destruct(it->second);
	}
	delete head;
}

int bin_search(int num, ll* cum_fcount){

	// cout<<num<<"hi"<<endl;
	
	int m = (no_of_bits+63)/64;
	int l = 0,r = m-1;
	while(r-l > 1){
		int mid = (l+r)>>1;
		if(cum_fcount[mid] <= num)l = mid;
		else r = mid;
	}
	if(cum_fcount[l] >= num)r = l;
	return r;
}

int P_suffixtree::fsum(){
	return fsum(root);
}

int P_suffixtree::fsum(node* z){
	int ind = bin_search(z->pst_ordr_rnk, cum_fcount);
	// cout<<ind<<" hey "<<no_of_bits<<endl;
	int ans = 64*ind + select_int(fcount_unary[ind],z->pst_ordr_rnk - (ind > 0)*cum_fcount[ind-1]);
	return (ans - z->pst_ordr_rnk);
}