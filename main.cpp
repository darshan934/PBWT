//#include "C:\Users\Darshan\Desktop\Research\p-pattern matching\tree.hpp"
//#include "C:\Users\Darshan\Desktop\Research\p-pattern matching\succint_tree.hpp"
#include <C:\Users\Darshan\Desktop\Research\p-pattern matching\tree.hpp>
#include <C:\Users\Darshan\Desktop\Research\p-pattern matching\succint_tree.hpp>
bool param[256];

int main(){
	string s;
	cin>>s;
	// s = "sandeep";
	for(int i = 0;i < 26;i++)param['a'+i] = 1;
	P_suffixtree mytree(s,param,26);
	cout<<mytree.fsum()<<endl;
	cout<<bp_size<<endl;
	construct();
	cout<<enclose(320)<<endl;
	cout<<lmost_leaf(315)<<endl; //debug
	destroy();
}
