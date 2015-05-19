#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>



using namespace std;

void getBigUnitigs(const string input, const string output, int length){
	ifstream in(input);
	ofstream out(output);
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
		int l(stoi(arg3));
		getBigUnitigs(arg1,arg2,l);
	}
}
