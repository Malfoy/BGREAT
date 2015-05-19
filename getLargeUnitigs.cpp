#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <algorithm>
#include <cmath>



using namespace std;

void getBigUnitigs(const string input, const string output, int length){
	ifstream in(input.c_str(),ios::in);
	ofstream out(output.c_str(),ios::out);
	string line;
	
	out<<">N"<<endl;
	
	while(!in.eof()){
		getline(in,line);	
		if(line.size()>length){
			line=line.substr(0,line.size()-1);
			transform(line.begin(),line.end(), line.begin(), ::toupper);
			out<<line;
		}
	}
}


int main(int argc, char ** argv)
{
	if(argc<4){
		cout<<"Need in, out and length..."<<endl;
	}
	if(argc==4){
		string arg1=argv[1];
		string arg2=argv[2];
		string arg3=argv[3];
		int l(atoi(arg3.c_str()));
		getBigUnitigs(arg1,arg2,l);
	}
}
