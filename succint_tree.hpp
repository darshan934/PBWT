#include <bits/stdc++.h>
typedef unsigned long long ll;

void build(int id , int x , int y);

void construct();

int left_deep(int id , int d);

int fwd_search(int ind , int d);

int right_deep(int id , int d);

int bwd_search(int ind , int d);

int sum_deep(int id , int j);

int sum(int i , int j);

int enclose(int i);

int parent(int i);

int rank_(int i);

int select1(int rnk);

int select0(int rnk);

int lmost_leaf(int i);

int rmost_leaf(int i);

void destroy();
