#include <bits/stdc++.h>
typedef unsigned long long ll;

int rank_int(ll v,int pos);
int select_int(ll v,int r);

class wavelet_tree{

	struct node{
		ll* bits;
		ll* cumulative;
		int n;
		int div;
		node *par , *left , *right;

		~node(){
			delete bits;
			delete cumulative;
		}
	};
	
	node *head;
	
	vector<int> arr;

public:

	// ll rank(unsigned char,ll);

	// ll rank(unsigned char,ll,node*);
	
	wavelet_tree(vector<int> &s);

	void construct(vector<int> &s, int l, int r, node* head);

	int twod(int l,int r,int c,int d);

	int twod(int x,int y,int l,int r,int c,int d,node*);

	~wavelet_tree();

	int bin_search0(int,node*);

	int bin_search1(int,node*);

	void destruct(node*);

	int predecessor(int, int);

	int predecessor(int, int, node*);
};