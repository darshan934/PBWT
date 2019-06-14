#include "tree.hpp"
#include "support_wt.hpp"
#include "succint_tree.hpp"
typedef unsigned long long ll;

// integer arrays for extra function , minimum and maximum values of Extra function.
int *e , *m , *M , *m_ind , *M_ind , *lf , *leaves;
int s = 16 , children = 3;
// each leaf node is a cluster of s bits from bp array and children is no of children of any parent node.
// this will be a complete k-ary tree (similar to heap)


//required functions - ith leaf , lmost and rmost leaf , fsum() , pcount , pre_order_rank_to_node , lca
//in zeronode() - int pre_rank = wt_zero_dep->predecessor(tmp->pre_ordr_rnk, pbwt[i]+1);  - what abt this and wt_bwt ?
//y = levelAncestor(Lx − 1, nodeDepth(lca(Lx;Lx − 1)) + 1) and fSum(x) = rankB(selectB(post-order(y), 0), 1)

//completed functions - parent

int len , cur = 1 , size;
int* tmp;

void build(int id , int x , int y){
	lf[id] = 0;
	if(x > y){
		lf[id] = -1;
		return;
	}
	if(x == y){
		e[id] = tmp[min(len - 1 , s * x + s - 1)];
		m[id] = INT_MAX;
		M[id] = INT_MIN;
		m_ind[id] = s*x;
		M_ind[id] = min(len , s*(x + 1)) - 1;
		for(int i = s*x;i < min(len , s*(x + 1));i++){
			m[id] = min(m[id] , tmp[i]);
			M[id] = max(M[id] , tmp[i]);
		}
		lf[id] = cur++;
		leaves[lf[id] - 1] = id;
		return;
	}
	int m1 = (2*x + y)/3 , m2 = (x + 2*y)/3;
	build(3*id - 1 , x , m1);
	build(3*id , m1 + 1 , m2);
	build(3*id + 1 , m2 + 1 , y);
	
	e[id] = e[3*id + 1];
	m[id] = min(m[3*id - 1] , min(m[3*id] , m[3*id + 1]));
	M[id] = max(M[3*id - 1] , min(M[3*id] , M[3*id + 1]));
	m_ind[id] = m_ind[3*id - 1];
	M_ind[id] = M_ind[3*id + 1];
}
// children of node i are 3*i - 1 , 3*i , 3 * i + 1
// parent of node i is (i + 1)/3

void construct(){
	len = bp_size;
	int k = (len + s - 1)/s; // no of leaf nodes
	size = 1;
	int cnt = 0;
	while(k){
		size *= 3;
		cnt++;
		k /= 3;
	}
	size /= 3;
	k = (len + s - 1)/s;
	if(k > size)size *= 3;
	size = (3 * size + 1)/2; // 3/2 * next nearest power of 3

	e = new int[size + 1];
	m = new int[size + 1];
	M = new int[size + 1];
	m_ind = new int[size + 1];
	M_ind = new int[size + 1];
	lf = new int[size + 1];
	leaves = new int[size + 1];

	tmp = new int[len];
	for(int i = 0;i < len;i++){
		k = (bp[i/64] >> (63 - i%64)) & 1;
		tmp[i] = 2*k - 1;
		if(i)tmp[i] += tmp[i - 1];
	}

	k = (len + s - 1)/s;
	build(1 , 0 , k - 1); // recursively builds the tree
	delete[] tmp;
}

int left_deep(int id , int d){
	if(lf[id]){
		// leaf node - iterate and return index
		int k = lf[id] - 1;
		int ans = -1 , sum = e[id] , tmp;
		for(int i = min(len , s*(k+1)) - 1;i >= s*k;i--){
			if(sum == d)ans = i;
			tmp = (bp[i/64] >> (63 - i%64)) & 1;
			sum += 1 - 2*tmp;
		}
		return ans;
	}
	for(int i = 3*id - 1;i < 3*id + 2;i++){
		if(lf[i] == -1)continue;
		if(m[i] <= d and M[i] >= d)return left_deep(i , d);
	}
}

int fwd_search(int ind , int d){
	int k = ind/s;
	//first check if it is in the same block - iterate - O(16)
	int sum = 0 , tmp;
	for(int i = ind;i < min(len , s*(k + 1));i++){
		tmp = (bp[i/64] >> (63 - i%64)) & 1;
		sum += 2*tmp - 1;
		if(sum == d)return i;
	}
	// if not in the same block check neighbouring blocks - big-o(children = 3)
	int id = leaves[k];
	if(k == (len - 1)/s)return -1;

	d += e[id] - sum;
	tmp = id;
	while(id > 1){
		id = (id + 1)/3;
		for(int i = tmp + 1;i < 3*id + 2;i++){
			if(lf[i] == -1)continue;
			if(m[i] <= d and M[i] >= d)return left_deep(i , d);
		}
		tmp = id;
	}
	return -1;
}

int right_deep(int id , int d){
	if(lf[id] > 0){
		// leaf node - iterate and return index
		int k = lf[id] - 1;
		int ans = -1 , sum = e[id] , tmp;
		for(int i = min(len , s*(k+1)) - 1;i >= s*k;i--){
			if(sum == d)return i + 1;
			tmp = (bp[i/64] >> (63 - i%64)) & 1;
			sum += 1 - 2*tmp;
		}
		if(sum == d)return s*k;
	}
	for(int i = 3*id + 1;i > 3*id - 2;i--){
		if(lf[i] == -1)continue;
		if(m[i] <= d and M[i] >= d)return right_deep(i , d);
	}
}

int bwd_search(int ind , int d){
	int k = ind/s;
	//first check if it is in the same block - iterate - O(16)
	int sum = 0 , tmp;
	for(int i = ind;i >= s*k;i--){
		tmp = (bp[i/64] >> (63 - i%64)) & 1;
		sum += 2*tmp - 1;
		if(sum == d)return i;
	}
	// if not in the same block check neighbouring blocks - O(children = 3)
	int id = leaves[k];
	if(k == 0)return -2;

	d = e[leaves[k - 1]] + sum - d; // k > 0 since it passed the previous if condition => id - 1 is also a leaf node
	tmp = id;
	while(id > 1){
		id = (id + 1)/3;
		for(int i = tmp - 1;i > 3*id - 2;i--){
			if(lf[i] == -1)continue;
			if(m[i] <= d and M[i] >= d)return right_deep(i , d);
		}
		tmp = id;
	}
	if(d == 0)return 0;	
	return -1;
}

int sum_deep(int id , int j){
	int sum = 0,tmp;
	if(lf[id]){
		for(int i = s*id;i <= j;i++){
			tmp = (bp[i/64] >> (63 - i%64)) & 1;
			sum += 2*tmp - 1;
		}
		return sum;
	}
	for(int i = 3*id - 1;i < 3*id + 2;i++){
		if(M_ind[i] < j){
			if(m_ind[i] == 0){
				sum += e[i];
			}
			else sum += e[i] - e[i - 1];
		}
		else return (sum += sum_deep(i , j));
	}
}

int sum(int i , int j){
	if(i > j)swap(i , j);
	int ind1 = i/s , ind2 = j/s;
	if(ind1 == ind2){ // same leaf node, iterate and return value
		int sum = 0 , tmp;
		for(int l = i;l <= j;l++){
			tmp = (bp[l/64] >> (63 - l%64)) & 1;
			sum += 2*tmp - 1;
		}
		return sum;
	}
	int sum = 0 , tmp;
	for(int l = i;l < min(len , s * (ind1 + 1));l++){
		tmp = (bp[l/64] >> (63 - l%64)) & 1;
		sum += 2*tmp - 1;
	}
	int id = ind1;
	tmp = id;
	while(id > 1){
		id = (id + 1)/3;
		for(int l = tmp + 1;l < 3 * id + 2;l++){
			if(M_ind[l] < j)sum += e[l] - e[l - 1];
			else return (sum += sum_deep(l , j));
		}
	}
}

int enclose(int i){
	int tmp = (bp[i/64] >> (63 - i%64)) & 1;
	//bwd_search(i , 2) if i is '(', bwd_search(i , 1) if i is ')'
	if(tmp)return fwd_search(i , 0);
	else return bwd_search(i , 0);
}

int parent(int i){
	//parent in the suffix tree represented by bp, not in the current tree(i.e. (i + 1)/3)
	int tmp = (bp[i/64] >> (63 - i%64)) & 1;
	//bwd_search(i , 2) if i is '(', bwd_search(i , 1) if i is ')'
	if(i == 0 or i == bp_size - 1)return -1;
	return bwd_search(i , 1 + tmp);
}

int rank_bit(int i){
	int j = i/64;
	int ans = 0;
	if(j)ans = cum_bp[j - 1];
	i %= 64;
	ans += rank_int(bp[j] , i + 1);
	return ans;
}

int rank_(int i , int j){
	if(i)return rank_bit(j);
	else return (j + 1 - rank_bit(j));
}

int select1(int rnk){
	int l = 0 , r = (len + 63)/64 - 1, mid;
	while(r - l > 1){
		mid = (l + r) >> 1;
		if(cum_bp[mid] < rnk)l = mid;
		else r = mid;
	}
	if(cum_bp[r] < rnk)return -1;
	if(cum_bp[l] >= rnk)r = l;
	if(r)rnk -= cum_bp[r - 1];
	return (select_int(bp[r] , rnk) - 1);
}

int select0(int rnk){
	int l = 0 , r = (len + 63)/64 - 1, mid;
	while(r - l > 1){
		mid = (l + r) >> 1;
		if(64 * (mid + 1) - cum_bp[mid] < rnk)l = mid;
		else r = mid;
	}

	mid = (len + 63)/64 - 1;
	if(r == mid){
		if(len - cum_bp[r] < rnk)return -1;
	}
	else {
		if(64 * (r + 1) - cum_bp[r] < rnk)return -1;
	}

	if(l == r){
		return (select_int(~bp[r] , rnk) - 1);
	}
	else {
		if(64 * (l + 1) - cum_bp[l] >= rnk)r = l;
		if(r)rnk -= 64 * r - cum_bp[r - 1];
		return (select_int(~bp[r] , rnk) - 1);
	}
}

// need rank and select for these 2 functions
int lmost_leaf(int i){//left most leaf of node at i
	int tmp = (bp[i/64] >> (63 - i%64)) & 1;
	if(tmp == 0)i = bwd_search(i , 0);

	tmp = (bp[(i + 1)/64] >> (63 - (i + 1)%64)) & 1;
	if(tmp == 0)return -1;//it is a leaf => no children

	int r = rank_(0 , i);
	// cout<<r<<" "<<rank_(1,i)<<endl;
	return select0(r + 1) - 1;
}

int rmost_leaf(int i){//right most leaf of node at i
	int tmp = (bp[i/64] >> (63 - i%64)) & 1;
	if(tmp == 1)i = fwd_search(i , 0);

	tmp = (bp[(i - 1)/64] >> (63 - (i - 1)%64)) & 1;
	if(tmp == 1)return -2;//it is a leaf => no children

	int r = rank_(1 , i);
	return select1(r);
}

void destroy(){
	delete[] e;
	delete[] m;
	delete[] M;
	delete[] m_ind;
	delete[] M_ind;
	delete[] lf;
	delete[] leaves;
}
