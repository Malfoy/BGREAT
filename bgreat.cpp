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
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
//#include <sparsehash/sparse_hash_map>
//#include <sparsehash/dense_hash_map>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include "Aligner.h"


using namespace std;


int main(int argc, char ** argv)
{
	srand (time(NULL));
	if(argc<4){
		cout<<"Need k , reads file , unitig file, errors allowed and number of thread ..."<<endl;
		exit(0);
	}
	if(argc==6){
		string arg1=argv[1];
		string reads=argv[2];
		string unitigs=argv[3];
		string errors=argv[4];
		string threads=argv[5];
		int ka=stoi(arg1);
		//int nCores=stoi(arg3);
		auto start=chrono::system_clock::now();

		Aligner supervisor(reads,unitigs,"paths","noOverlap.fa","notAligned.fa",ka,stoi(threads),stoi(errors));
		supervisor.indexUnitigs();
		supervisor.alignAll(true);
		auto end=chrono::system_clock::now();auto waitedFor=end-start;cout<<"Last for "<<chrono::duration_cast<chrono::seconds>(waitedFor).count()<<" seconds"<<endl<<endl;
	}

	if(argc==7){
		string arg1=argv[1];
		string reads=argv[2];
		string unitigs=argv[3];
		string errors=argv[4];
		string threads=argv[5];
		string brute=argv[6];
		if(brute=="-b"){
			int ka=stoi(arg1);
			//int nCores=stoi(arg3);
			auto start=chrono::system_clock::now();

			Aligner supervisor(reads,unitigs,"paths","noOverlap.fa","notAligned.fa",ka,stoi(threads),stoi(errors));
			supervisor.indexUnitigs();
			supervisor.alignAll(false);
			auto end=chrono::system_clock::now();auto waitedFor=end-start;cout<<"Last for "<<chrono::duration_cast<chrono::seconds>(waitedFor).count()<<" seconds"<<endl<<endl;
		}else{
			cout<<"not recognized option "<<brute<<endl;
			cout<<"Did you mean -b for bruteforce ?"<<endl;
		}

	}
}
