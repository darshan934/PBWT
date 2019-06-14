#include "tree.hpp"
#include "support_wt.hpp"

int rank_int(ll v,int pos){
	ll r;       // Resulting rank of bit at pos goes here.
	// Shift out bits after given position.
	r = v >> (sizeof(v) * CHAR_BIT - pos);
	// Count set bits in parallel.
	// r = (r & 0x5555...) + ((r >> 1) & 0x5555...);
	r = r - ((r >> 1) & ~0UL/3);
	// r = (r & 0x3333...) + ((r >> 2) & 0x3333...);
	r = (r & ~0UL/5) + ((r >> 2) & ~0UL/5);
	// r = (r & 0x0f0f...) + ((r >> 4) & 0x0f0f...);
	r = (r + (r >> 4)) & ~0UL/17;
	// r = r % 255;
	r = (r * (~0UL/255)) >> ((sizeof(v) - 1) * CHAR_BIT);
	return r;
}

int select_int(ll v,int r){
	if(r == 0)return 0;
    // Input value to find position with rank r.
    // Input: bit's desired rank [1-64].
    int s;      // Output: Resulting position of bit with rank r [1-64]
    ll a, b, c, d; // Intermediate temporaries for bit count.
    int t;      // Bit count temporary.

    // Do a normal parallel bit count for a 64-bit integer,
    // but store all intermediate steps.
    // a = (v & 0x5555...) + ((v >> 1) & 0x5555...);
    a =  v - ((v >> 1) & ~0UL/3);
    // b = (a & 0x3333...) + ((a >> 2) & 0x3333...);
    b = (a & ~0UL/5) + ((a >> 2) & ~0UL/5);
    // c = (b & 0x0f0f...) + ((b >> 4) & 0x0f0f...);
    c = (b + (b >> 4)) & ~0UL/0x11;
    // d = (c & 0x00ff...) + ((c >> 8) & 0x00ff...);
    d = (c + (c >> 8)) & ~0UL/0x101;
    t = (d >> 32) + (d >> 48);
    // Now do branchless select!
    s  = 64;
    // if (r > t) {s -= 32; r -= t;}
    s -= ((t - r) & 256) >> 3; r -= (t & ((t - r) >> 8));
    t  = (d >> (s - 16)) & 0xff;
    // if (r > t) {s -= 16; r -= t;}
    s -= ((t - r) & 256) >> 4; r -= (t & ((t - r) >> 8));
    t  = (c >> (s - 8)) & 0xf;
    // if (r > t) {s -= 8; r -= t;}
    s -= ((t - r) & 256) >> 5; r -= (t & ((t - r) >> 8));
    t  = (b >> (s - 4)) & 0x7;
    // if (r > t) {s -= 4; r -= t;}
    s -= ((t - r) & 256) >> 6; r -= (t & ((t - r) >> 8));
    t  = (a >> (s - 2)) & 0x3;
    // if (r > t) {s -= 2; r -= t;}
    s -= ((t - r) & 256) >> 7; r -= (t & ((t - r) >> 8));
    t  = (v >> (s - 1)) & 0x1;
    // if (r > t) s--;
    s -= ((t - r) & 256) >> 8;
    s = 65 - s;
    return s;
}

wavelet_tree::wavelet_tree(vector<int> &s){
	set<int> charset;
	ll n = s.size();
	for(int i = 0;i < n;i++){
		charset.insert(s[i]);
	}
	n = charset.size();
	arr.resize(n);
	int j = 0;
	for(auto it = charset.begin();it != charset.end();it++){
		arr[j++] = *it;
	}
	head = new node;
	head->par = NULL;
	construct(s,0,n-1,head);
}

void wavelet_tree::construct(vector<int> &s, int l, int r, node* head){
	head->n = s.size();
	int m = (head->n + 63)/64;
	head->bits = new ll[m];
	head->cumulative = new ll[m];
	for(ll i = 0;i < m;i++){
		head->cumulative[i] = 0;
		head->bits[i] = 0;
	}
	vector<int> sl, sr;
	ll mid = (l+r)/2;
	head->div = arr[mid];
	for(ll i = 0;i < head->n;i++){
		if(s[i] <= head->div){
			head->bits[i/64] &= ~(1ll<<(63-i%64));
			sl.push_back(s[i]);
		}
		else{
			head->bits[i/64] |= (1ll<<(63-i%64));
			head->cumulative[i/64]++;
			sr.push_back(s[i]);
		}
	}
	for(ll i = 1;i < m;i++){
		head->cumulative[i] += head->cumulative[i-1];
	}
	// if(l == r)cout<<head->div<<" "<<l<<" "<<r<<" "<<head->n<<endl;
	head->left = head->right = NULL;
	if(l == r)return;

	head->left = new node;
	head->right = new node;
	head->left->par = head;
	head->right->par = head;
	construct(sl,l,mid,head->left);
	construct(sr,mid+1,r,head->right);
}

int wavelet_tree::twod(int l,int r,int c,int d){
	if(r > head->n) r = head->n;
	int ans = twod(0,arr.size()-1,l,r,c,d,head);
	return ans;
}

int wavelet_tree::twod(int x,int y,int l,int r,int c,int d,node* head){
	int left = 0,right = 0;
	int ans = 0;
	if(l > r)return 0;
	if(c <= arr[x] and d >= arr[y])return r-l+1;
	if(d < arr[x] or c > arr[y])return 0;
	if(head->left == NULL){
		if(head->div <= d and head->div >= c)return r-l+1;
		else return 0;
	}
	int mid = (x+y)/2;
	if(c <= head->div){
		if(l == 1){
			left = 1;
		}
		else{
			left = l - 2 - (l-2)%64;
			left += (l-2)%64 + 1 - rank_int(head->bits[(l-2)/64],(l-2)%64 + 1);
			left++;
			if((l-2)/64)left -= head->cumulative[(l-2)/64 - 1];
		}
		right = r - 1 - (r-1)%64;
		right += (r-1)%64 + 1 - rank_int(head->bits[(r-1)/64],(r-1)%64 + 1);
		if((r-1)/64)right -= head->cumulative[(r-1)/64 - 1];
		ans += twod(x,mid,left,right,c,d,head->left);
	}
	left = right = 0;
	if(d > head->div){
		if((r-1)/64)right = head->cumulative[(r-1)/64 - 1];
		right += rank_int(head->bits[(r-1)/64],(r-1)%64 + 1);
		if(l == 1)left = 1;
		else{
			if((l-2)/64)left = head->cumulative[(l-2)/64 - 1];
			left += 1 + rank_int(head->bits[(l-2)/64],(l-2)%64 + 1);
		}
		ans += twod(mid+1,y,left,right,c,d,head->right);
	}
	return ans;
}

wavelet_tree::~wavelet_tree(){
	destruct(head);
}

void wavelet_tree::destruct(node* head){
	if(head->left != NULL)destruct(head->left);
	if(head->right != NULL)destruct(head->right);
	delete head;
}

int wavelet_tree::bin_search0(int num,node *head){
	ll n = head->n;
	ll m = (n+63)/64;
	int extra = 64*m - n;
	ll l = 0,r = m-1;
	while(r-l > 1){
		ll mid = (l+r)>>1;
		if(64*(mid+1) - head->cumulative[mid] - extra*(mid == m-1) <= num)l = mid;
		else r = mid;
	}
	if(64*(l+1) - head->cumulative[l] >= num)r = l;
	return r;
}

int wavelet_tree::bin_search1(int num,node *head){
	ll n = head->n;
	ll m = (n+63)/64;
	ll l = 0,r = m-1;
	while(r-l > 1){
		ll mid = (l+r)>>1;
		if(head->cumulative[mid] <= num)l = mid;
		else r = mid;
	}
	if(head->cumulative[l] >= num)r = l;
	return r;
}

int wavelet_tree::predecessor(int ind, int val){ // returns 0 if no predecessor
	return predecessor(ind, val, head);
}

int wavelet_tree::predecessor(int ind, int val, node* head){ // find the right most 'position' whose value < val and position <= ind & 1-based indexing
	if(head->left == NULL){
		if(head->div < val)return ind;
		else return 0;
	}

	if(val <= head->div){ // need to check only the left branch of tree
		int left;
		left = ind -1 - (ind-1)%64;
		left += (ind-1)%64 + 1 - rank_int(head->bits[(ind-1)/64],(ind-1)%64+1);
		if((ind-1)/64)left -= head->cumulative[(ind-1)/64 - 1];
		int idx = predecessor(left, val, head->left);
		ind = bin_search0(idx, head);
		ind = 64*ind + select_int(~(head->bits[ind]),idx - 64*(ind) + (ind > 0)*head->cumulative[ind-1]);
		return ind;
	}
	else{ // the rightmost character less than head->div and right side predecessor compete
		int left, right;
		right = rank_int(head->bits[(ind-1)/64],(ind-1)%64+1);
		if((ind-1)/64)right += head->cumulative[(ind-1)/64-1];
		int idx = predecessor(right, val, head->right);
		int ans = bin_search1(idx, head);
		ans = 64*ans + select_int(head->bits[ans],idx - (ans > 0)*head->cumulative[ans-1]);

		left = ind -1 - (ind-1)%64;
		left += (ind-1)%64 + 1 - rank_int(head->bits[(ind-1)/64],(ind-1)%64+1);
		ind = bin_search0(left, head);
		left = 64*ind + select_int(~(head->bits[ind]),left - 64*(ind) + (ind > 0)*head->cumulative[ind-1]);
		return max(left, ans);
	}
}

 int main(){
 	vector<int> v;
 	int n;
 	int x,y;
 	cin>>n;
 	for(int i = 0;i < n;i++){
 		cin>>x;
 		v.push_back(x);
 	}
 	wavelet_tree test(v);
 	cin>>x>>y;
 	int ind = test.predecessor(x,y);
 	cout<<ind<<" "<<v[ind-1]<<" "<<y<<endl;
}
