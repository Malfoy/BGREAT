//  Aligner.cpp
//  BGREAT
//
//  Created by malfoy on 09/04/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "Aligner.h"
#include <thread>
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



using namespace std;

char revcomp (char s) {
	if (s == 'A') return 'T';
	else if (s == 'C') return 'G';
	else if (s == 'G') return 'C';
	else if (s == 'T') return 'A';
	// else if (s == 'a') return 't';
	// else if (s == 'c') return 'g';
	// else if (s == 'g') return 'c';
	// else if (s == 't') return 'a';
	cout<<"wtf rc"<<s<<endl;
	cin.get();
	return 'X';
}

string reverseComplement (const string& s){
	string rc;
	for (int i = (int)s.length() - 1; i >= 0; i--){
		rc += revcomp(s[i]);
	}
	return rc;

}

//Parse reads
void Aligner::getReads(vector<string>& reads, uint n){
	reads={};
	string read,header,inter;
	char c;

	for(int i(0);i<n;++i){
		getline(readFile,header);
		getline(readFile,read);
	point:
		c=readFile.peek();
		if(c=='>'){
			if(read.size()>2){
				bool fail(false);
				for(uint j(0);(j)<read.size();++j){
					if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G'){
						fail=true;
						break;
					}
				}
				if(!fail){
					reads.push_back(read);
				}
			}
			read="";
		}else{
			if(!readFile.eof()){
				getline(readFile,inter);
				read+=inter;
				goto point;
			}else{
				if(read.size()>2){
					bool fail(false);
					for(uint j(0);(j)<read.size();++j){
						if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G'){
							fail=true;
							break;
						}
					}
					if(!fail){
						reads.push_back(read);
					}
				}
				return;
			}
		}
	}
}


uint64_t Aligner::getRepresentNum(const string& str){
	string rc(reverseComplement(str));
	uint64_t a(str2num(str));
	uint64_t b(str2num(rc));
	if(b<=a){
		return b;
	}
	return a;
}

char randNuc(){
	switch (rand()%4){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	return 'A';
}

string mutate(string read,int n){
	for(int i(0); i<n; i++){
		int position(rand()%read.size());
		read[position]=randNuc();
	}
	return read;
}


string Aligner::num2str(uint64_t num){
	string str;
	int nuc;
	for(int i(0);i<k-1;i++){
		nuc=num%4;
		switch (nuc){
			case 0:str.push_back('A');break;
			case 1:str.push_back('C');break;
			case 2:str.push_back('G');break;
			case 3:str.push_back('T');break;
		}
		num>>=2;
	}
	reverse( str.begin(), str.end());
	return str;
}


uint64_t Aligner::str2num(const string& str){
	uint64_t res(0);
	for(int i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			case 'T':res+=3;break;
		}

	}
	return res;

}


char nuc2int(char c){
	switch(c){
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	cout<<"failnuc2int"<<endl;
	cout<<c<<endl;
	cin.get();
	return 0;
}

uint64_t nuc2intrc(char c){
	switch(c){
		case 'A': return 3;
		case 'C': return 2;
		case 'G': return 1;
		case 'T': return 0;
	}
	cout<<"failnuc2intRC"<<endl;
	cout<<c<<endl;
	cin.get();
	return 0;
}

//vector<uint32_t> Aligner::getBegin(uint64_t bin){
//	string str(num2str(bin));
//	string rc(reverseComplement(str));
//	if(str<=rc){
//		if(left.count(bin)!=0){
//			return left[bin];
//		}else{
//			return {};
//		}
//	}else{
//		uint64_t num(str2num(rc));
//		if(right.count(num)!=0){
//			return right[num];
//		}else{
//			return {};
//		}
//	}
//	return {};
//}
//
//
//
//
//vector<uint32_t> Aligner::getEnd(uint64_t bin){
//	string str(num2str(bin));
//	string rc(reverseComplement(str));
//	if(str<=rc){
//		if(right.count(bin)!=0){
//			return right[bin];
//		}else{
//			return {};
//		}
//	}else{
//		uint64_t num(str2num(rc));
//		if(left.count(num)!=0){
//			return left[num];
//		}else{
//			return {};
//		}
//	}
//	return {};
//}

vector<pair<string,uNumber>> Aligner::getEnd(uint64_t bin){
	vector<pair<string,uNumber>> result;
	string str(num2str(bin));
	string rc(reverseComplement(str));
	if(str<=rc){
		if(right.count(bin)!=0){
			for(uint i(0); i<right[bin].size(); ++i){
				string unitig(getUnitig(right.at(bin)[i]));
				if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
					result.push_back({unitig,right.at(bin)[i]});
				}else{
					result.push_back({reverseComplement(unitig),-right.at(bin)[i]});
				}
			}
		}else{
			return {};
		}
	}else{
		uint64_t num(str2num(rc));
		if(left.count(num)!=0){
			for(uint i(0); i<left[num].size(); ++i){
				string unitig(getUnitig(left.at(num)[i]));
				if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
					result.push_back({unitig,left.at(num)[i]});
				}else{
					result.push_back({reverseComplement(unitig),-left.at(num)[i]});
				}
			}
		}else{
			return {};
		}
	}
	return result;
}


vector<pair<string,uNumber>> Aligner::getBegin(uint64_t bin){
	vector<pair<string,uNumber>> result;
	string str(num2str(bin));
	string rc(reverseComplement(str));
	if(str<=rc){
		if(left.count(bin)!=0){
			for(uint i(0); i<left[bin].size(); ++i){
				string unitig(getUnitig(left.at(bin)[i]));
				if(str2num(unitig.substr(0,k-1))==bin){
					result.push_back({unitig,left.at(bin)[i]});

				}else{
					result.push_back({reverseComplement(unitig),-left.at(bin)[i]});
				}
			}
		}else{
			return {};
		}
	}else{
		uint64_t num(str2num(rc));
		if(right.count(num)!=0){
			for(uint i(0); i<right[num].size(); ++i){
				string unitig(getUnitig(right.at(num)[i]));
				if(str2num(unitig.substr(0,k-1))==bin){
					result.push_back({unitig,right.at(num)[i]});
				}else{
					result.push_back({reverseComplement(unitig),-right.at(num)[i]});
				}
			}
		}else{
			return {};
		}
	}
	return result;
}


vector<uNumber> Aligner::getEndNumber(uint64_t bin){
	vector<uNumber> result;
	string str(num2str(bin));
	string rc(reverseComplement(str));
	if(str<=rc){
		if(right.count(bin)!=0){
			for(uint i(0); i<right[bin].size(); ++i){
				result.push_back(right.at(bin)[i]);
			}
		}else{
			return {};
		}
	}else{
		uint64_t num(str2num(rc));
		if(left.count(num)!=0){
			for(uint i(0); i<left[num].size(); ++i){
				result.push_back(left.at(num)[i]);
			}
		}else{
			return {};
		}
	}
	return result;
}


vector<uNumber> Aligner::getBeginNumber(uint64_t bin){
	vector<uNumber> result;
	string str(num2str(bin));
	string rc(reverseComplement(str));
	if(str<=rc){
		if(left.count(bin)!=0){
			for(uint i(0); i<left[bin].size(); ++i){
				result.push_back(left.at(bin)[i]);
			}
		}else{
			return {};
		}
	}else{
		uint64_t num(str2num(rc));
		if(right.count(num)!=0){
			for(uint i(0); i<right[num].size(); ++i){
				result.push_back(right.at(num)[i]);
			}
		}else{
			return {};
		}
	}
	return result;
}


int missmatchNumber(const string& seq1, const string& seq2, unsigned int n){
	if(seq1.size()!=seq2.size()){
		cout<<seq1.size()<<" "<<seq2.size()<<endl;
		cout<<"missmatchNumber bug"<<endl;
		exit(0);
	}

	char miss(0);
	for(size_t i(0); i<seq2.size(); ++i){
		if(seq2[i]!=seq1[i]){
			if(++miss>n){
				return miss;
			}
		}
	}

	return miss;
}


int errorNumber(const string& seq1, const string& seq2, unsigned int n){
	if(seq1.size()>seq2.size()+n or seq2.size()>seq1.size()+n){
		return n+1;
	}
	size_t size(max(seq1.size(),seq2.size()));

	char miss(0);
	for(size_t i(0); i<size; ++i){
		if(seq2[i]!=seq1[i]){
			if(++miss>n){
				return miss;
			}
		}
	}

	return miss;
}


//vector<uint32_t> Aligner::alignRead(const string& read, bool& overlapFound){
//	vector<uint32_t> path;
//	vector<pair<uint64_t,uint>> listOverlap(getListOverlap(read));
//
//	if(listOverlap.empty()){
//		overlapFound=false;
//		return {};
//	}
//
//	overlapFound=true;
//
//	if(!checkBegin(read,listOverlap[0],path)){return {};}
//	//	cout<<"beginsuccess"<<endl;
//
//	for(size_t i(0); i<listOverlap.size()-1; ++i){
//		if(listOverlap[i+1].second-listOverlap[i].second>k){
//			if(!checkPair(listOverlap[i], listOverlap[i+1], read, path)){return {};}
//		}
//	}
//	//	cout<<"insidesuccess"<<endl;
//
//	if(!checkEnd(read,listOverlap[listOverlap.size()-1],path)){return {};}
//	//	cout<<"endsuccess"<<endl;
//
//	return path;
//}

//
//vector<uint32_t> Aligner::alignRead2(const string& read, bool& overlapFound){
//	vector<uint32_t> path;
//	vector<pair<uint64_t,uint>> listOverlap(getListOverlap(read));
//	if(listOverlap.empty()){
//		noOverlapRead++;
//		readNumber++;
//		return {};
//	}
//	overlaps+=listOverlap.size();
//	overlapFound=true;
//	for(int start(0); start<(int)listOverlap.size(); ++start){
//		if(checkBegin2(read,listOverlap[start],path)){
//			for(int end((int)listOverlap.size()-1); end>=start; --end){
//				if(checkEnd2(read,listOverlap[end],path)){
//					if(cover2(read,listOverlap,start,end,path)){
//						alignedRead++;
//						readNumber++;
//						return path;
//
//					}
//					//					if(mapOnLeftEndAll(read, path, listOverlap[end])){
//					//						alignedRead++;
//					//						readNumber++;
//					//						successMR2++;
//					//						return path;
//					//					}
//					break;
//				}
//			}
//			//			if(mapOnRightEndAll(read, path, listOverlap[start])){
//			//				alignedRead++;
//			//				readNumber++;
//			//				successMR2++;
//			//				return path;
//			//			}
//			break;
//		}
//	}
//	notAligned++;
//	readNumber++;
//	return {};
//}


vector<uNumber> Aligner::alignRead3(const string& read, bool& overlapFound, int errors, bool rc){
	vector<pair<uint64_t,uint>> listOverlap(getListOverlap(read));

	if(listOverlap.empty()){
		noOverlapRead++;
		readNumber++;
		return {};
	}

	overlaps+=listOverlap.size();
	overlapFound=true;
	for(int start(0); start<=min(tryNumber-1,(int)listOverlap.size()-1); ++start){
		vector<uNumber> pathBegin;
		int errorBegin(checkBeginExhaustive(read,listOverlap[start],pathBegin,errors));
		if(errorBegin<=errors){
			for(int end((int)listOverlap.size()-1); end>=max(start,(int)listOverlap.size()-tryNumber); --end){
				vector<uNumber> pathEnd;
				int errorsEnd(checkEndExhaustive(read,listOverlap[end],pathEnd,errors-errorBegin));
				if(errorsEnd+errorBegin<=errors){
					vector<uNumber> pathCover;
					bool ended(false);
					int errorCover(coverGreedy(read,listOverlap,start,end,pathCover,errors-errorsEnd-errorBegin,ended));

					if(errorCover<=errors-errorsEnd-errorBegin){
						alignedRead++;
						readNumber++;
						pathBegin.insert(pathBegin.end(), pathCover.begin(),pathCover.end());
						if (!ended) {
							pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
						}
						//						cout<<"read : "<<read<<endl;
						//												string recoveredpath=recoverPath(pathBegin,read.size());
						//
						//												if(missmatchNumber(read, recoveredpath, 1000)>errorsMax){
						//													cout<<missmatchNumber(read, recoveredpath, 1000)<<endl;
						//													cout<<read<<endl<<recoveredpath<<endl;
						//													cout<<"bug"<<endl;
						//													exit(0);
						//												}
						return pathBegin;
					}
				}
			}
		}
	}
	if(!rc){
		return alignRead3(reverseComplement(read), overlapFound,errors, true);
	}

	notAligned++;
	readNumber++;
	return {};
}


vector<uNumber> Aligner::alignRead4(const string& read, bool& overlapFound, int errors){
	vector<pair<uint64_t,uint>> listOverlap(getListOverlap(read));

	if(listOverlap.empty()){
		noOverlapRead++;
		readNumber++;
		return {};
	}

	overlaps+=listOverlap.size();
	overlapFound=true;
	for(int start(0); start<listOverlap.size(); ++start){
		vector<uNumber> pathBegin;
		int errorBegin(checkBeginExhaustive(read,listOverlap[start],pathBegin,errors));
		if(errorBegin<=errors){
			vector<uNumber> pathEnd;
			int errorsEnd(checkEndExhaustive(read,listOverlap[start],pathEnd,errors-errorBegin));
			if(errorsEnd+errorBegin<=errors){
				alignedRead++;
				readNumber++;
				pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
				//				string recoveredpath=recoverPath(pathBegin,read.size());
				//				if(missmatchNumber(read, recoveredpath, errorsMax)>errorsMax){
				//					cout<<recoveredpath<<endl;
				//					cout<<read<<endl;
				//						cout<<"bug"<<endl;
				//					exit(0);
				//					}
				return pathBegin;
			}
		}
	}
	notAligned++;
	readNumber++;
	return {};
}


string compactionBegin(const string& seq1,const string& seq2, size_t k){
	size_t s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}

	string rc2(reverseComplement(seq2));
	//	string rc1(reversecomplement(seq1));

	string beg1(seq1.substr(0,k));
	string end2(seq2.substr(s2-k,k));
	if(beg1==end2){
		return seq2+(seq1.substr(k));
	}

	string endrc2(rc2.substr(s2-k,k));
	if(beg1==endrc2){
		return rc2+(seq1.substr(k));
	}
	return "";
}


string compactionEnd(const string& seq1,const string& seq2, size_t k){
	size_t s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}

	string rc2(reverseComplement(seq2));
	//		string rc1(reversecomplement(seq1));

	string end1(seq1.substr(s1-k,k));
	string beg2(seq2.substr(0,k));
	//	cout<<end1<<endl<<beg2<<endl;
	if(end1==beg2){
		return seq1+(seq2.substr(k));
	}

	string begrc2(rc2.substr(0,k));
	if(end1==begrc2){
		return seq1+(rc2.substr(k));
	}
	return "";
}


string Aligner::recoverPath(vector<uNumber>& numbers,int size){
	//	for(size_t i(1); i<numbers.size(); ++i){
	//		string unitig=getUnitig(numbers[i]);
	////		cout<<unitig<<endl;
	//	}
	int offset(numbers[0]);
	string path(getUnitig(numbers[1]));
	//	cout<<"path"<<path<<endl;
	for(size_t i(2); i<numbers.size(); ++i){
		string unitig=getUnitig(numbers[i]);
		string inter=compactionEnd(path, unitig, k-1);
		if(inter.empty()){
			cout<<"bug compaction"<<endl;
			cout<<path<<" "<<unitig<<endl;
			exit(0);
		}else{
			path=inter;
			//			cout<<"path"<<path<<endl;
		}

	}
	path=path.substr(offset,size);
	return path;
}





//bool Aligner::cover2(const string& read, const vector<pair<uint64_t,uint>>& listOverlap, const size_t start, size_t end, vector<uNumber>& path){
//	if(start==end){return true;}
//
//	for(size_t i(1); i+start<=end; ++i){
//		if(checkPair(listOverlap[start], listOverlap[i+start], read, path)){
//			//						if(cover2(read, listOverlap, i+start, end, path)){
//			//							return true;
//			//						}
//			return (cover2(read, listOverlap, i+start, end, path));
//		}
//	}
//
//	bool endBool(false);
//	size_t next(mapOnRightAll(read,path,listOverlap[start],listOverlap,endBool,start));
//	if(endBool){
//		successMR++;
//		return true;
//	}
//	if(next>start){
//		if(cover2(read, listOverlap, next, end, path)){
//			successMR++;
//			return true;
//		};
//	}
//	return false;
//}

//seem OK, not exaustive
int Aligner::coverGreedy(const string& read, const vector<pair<uint64_t,uint>>& listOverlap, const size_t start, size_t end, vector<uNumber>& path, int  errors, bool& ended){
	if(start==end){return 0;}
	int indice(0);
	int minMissmatch(errors+1);
	uNumber number2keep;

	for(size_t i(1); i+start<=end and i<=tryNumber; ++i){
		uNumber number(0);
		//		cout<<i<<endl;
		//		cout<<listOverlap[i+start].first<<" "<<listOverlap[i+start].second<<endl;
		int missmatch(checkPair2(listOverlap[start], listOverlap[i+start], read, number,errors));
		if(missmatch<minMissmatch){
			indice=i;
			number2keep=number;
			minMissmatch=missmatch;
		}
	}

	if(minMissmatch<=errors){
		path.push_back(number2keep);
		return minMissmatch+coverGreedy(read, listOverlap, indice+start, end, path, errors-minMissmatch,ended);
	}

	auto pair(mapOnRight2(read,path,listOverlap[start],listOverlap,ended,start,errors));
	if(pair.second<=errors){
		if(ended){
			successMR2++;
			return pair.second;
		}

		if(pair.first>start){
			//			cout<<pair.first<<endl;
			//			cout<<num2str(listOverlap[pair.first].first)<<endl;
			int errorrecur(coverGreedy(read, listOverlap, pair.first, end, path,errors-pair.second,ended));
			if(errorrecur+pair.second<=errors){
				successMR2++;
				return errorrecur+pair.second;
			};
		}
	}
	return errors+1;
}


//bool Aligner::mapOnRightEnd(const string &read, vector<uint32_t>& path, pair<uint64_t, uint> overlap){
//	auto rangeUnitigs(overlap2unitigs.find(overlap.first)->second);
//	string unitig;
//	string readLeft(read.substr(overlap.second));
//	if(readLeft.size()<=k-1){return true;}
//
//	for(uint i(0); i<rangeUnitigs.size(); ++i){
//		unitig=getUnitig(rangeUnitigs[i]);
//		//case the rest of the read is too small
//		if(readLeft.size() <= unitig.size()){
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//		}else{
//			//case the read is big enough we want to recover a true overlap
//			if(missmatch(unitig, read.substr(overlap.second,unitig.size()), (unsigned char)(unitig.size()*errorRate))){
//				uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//				return (mapOnRightEnd(read, path, {overlapNum,overlap.second+unitig.size()-k+1}));
//			}else{
//				unitig=reverseComplement(unitig);
//				if(missmatch(unitig, read.substr(overlap.second,unitig.size()), (unsigned char)(unitig.size()*errorRate))){
//					uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//					return (mapOnRightEnd(read, path, {overlapNum,overlap.second+unitig.size()-k+1}));
//				}
//			}
//		}
//	}
//	return false;
//}


//bool Aligner::mapOnRightEndAll(const string &read, vector<uint32_t>& path, pair<uint64_t, uint> overlap){
//	auto rangeUnitigs(overlap2unitigs.find(overlap.first)->second);
//	string unitig;
//	string readLeft(read.substr(overlap.second));
//	if(readLeft.size()<=k-1){return true;}
//
//	for(uint i(0); i<rangeUnitigs.size(); ++i){
//		unitig=getUnitig(rangeUnitigs[i]);
//		//case the rest of the read is too small
//		if(readLeft.size() <= unitig.size()){
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//		}else{
//			//case the read is big enough we want to recover a true overlap
//			if(missmatch(unitig, read.substr(overlap.second,unitig.size()), 1+(unsigned char)(unitig.size()*errorRate))){
//				uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//				if(mapOnRightEndAll(read, path, {overlapNum,overlap.second+unitig.size()-k+1})){
//					return true;
//				}
//				//				return (mapOnRightEnd(read, path, {overlapNum,overlap.second+unitig.size()-k+1}));
//			}else{
//				unitig=reverseComplement(unitig);
//				if(missmatch(unitig, read.substr(overlap.second,unitig.size()), 1+(unsigned char)(unitig.size()*errorRate))){
//					uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//					if( mapOnRightEndAll(read, path, {overlapNum,overlap.second+unitig.size()-k+1})){
//						return true;
//					}
//					//					return (mapOnRightEnd(read, path, {overlapNum,overlap.second+unitig.size()-k+1}));
//				}
//			}
//		}
//	}
//	return false;
//}


//
//bool Aligner::mapOnLeftEndAll(const string &read, vector<uint32_t>& path,const pair<uint64_t, uint>& overlap){
//	auto rangeUnitigs(overlap2unitigs.find(overlap.first)->second);
//	string unitig;
//	string readLeft(read.substr(0,overlap.second));
//	if(readLeft.size()==0){return true;}
//
//	for(uint i(0); i<rangeUnitigs.size(); ++i){
//		unitig=getUnitig(rangeUnitigs[i]);
//		//case the rest of the read is too small
//		if(readLeft.size()+k-1 <= unitig.size()){
//			if(missmatch(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//		}else{
//			//case the read is big enough we want to recover a true overlap
//			if(missmatch(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), (unsigned char)((unitig.size()-k+1)*errorRate))){
//				uint64_t overlapNum(getRepresentNum(unitig.substr(0,k-1)));
//				if(mapOnLeftEndAll(read, path, {overlapNum,overlap.second-(unitig.size()-k+1)})){
//					return true;
//				}
//			}else{
//				unitig=reverseComplement(unitig);
//				if(missmatch(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), (unsigned char)((unitig.size()-k+1)*errorRate))){
//					uint64_t overlapNum(getRepresentNum(unitig.substr(0,k-1)));
//					if(mapOnLeftEndAll(read, path, {overlapNum,overlap.second-(unitig.size()-k+1)})){
//						return true;
//					}
//				}
//			}
//		}
//	}
//	return false;
//}


//bool Aligner::mapOnLeftEndFirst(const string &read, vector<uint32_t>& path, const  pair<uint64_t, uint>& overlap){
//	auto rangeUnitigs(overlap2unitigs.find(overlap.first)->second);
//	string unitig;
//	string readLeft(read.substr(0,overlap.second));
//	if(readLeft.size()==0){return true;}
//
//	for(uint i(0); i<rangeUnitigs.size(); ++i){
//		unitig=getUnitig(rangeUnitigs[i]);
//		//case the rest of the read is too small
//		if(readLeft.size()+k-1 <= unitig.size()){
//			if(missmatch(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//		}else{
//			//case the read is big enough we want to recover a true overlap
//			if(missmatch(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), (unsigned char)((unitig.size()-k+1)*errorRate))){
//				uint64_t overlapNum(getRepresentNum(unitig.substr(0,k-1)));
//				return(mapOnLeftEndFirst(read, path, {overlapNum,overlap.second-(unitig.size()-k+1)}));
//			}else{
//				unitig=reverseComplement(unitig);
//				if(missmatch(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), (unsigned char)((unitig.size()-k+1)*errorRate))){
//					uint64_t overlapNum(getRepresentNum(unitig.substr(0,k-1)));
//					return(mapOnLeftEndFirst(read, path, {overlapNum,overlap.second-(unitig.size()-k+1)}));
//				}
//			}
//		}
//	}
//	return false;
//}
//
//
//bool Aligner::mapOnLeftEndMax(const string &read, vector<uint32_t>& path, const  pair<uint64_t, uint>& overlap){
//	string unitig;
//	string readLeft(read.substr(0,overlap.second));
//	if(readLeft.size()==0){return true;}
//	auto rangeUnitigs(overlap2unitigs.find(overlap.first)->second);
//	double minimum(100);
//	int minimumindice(9);
//
//	for(uint i(0); i<rangeUnitigs.size(); ++i){
//		unitig=getUnitig(rangeUnitigs[i]);
//		//case the rest of the read is too small
//		if(readLeft.size()+k-1 <= unitig.size()){
//			if(missmatch(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//		}else{
//			//case the read is big enough we want to recover a true overlap
//			double miss(missMatchPercent(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1))));
//			if(miss<minimum){
//				minimum=miss;
//				minimumindice=2*i;
//			}
//			unitig=reverseComplement(unitig);
//			miss=(missMatchPercent(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1))));
//			if(miss<minimum){
//				minimum=miss;
//				minimumindice=2*i+1;
//			}
//		}
//	}
//	if(minimum<errorRate){
//		unitig=getUnitig(rangeUnitigs[minimumindice/2]);
//		if(minimumindice%2==1){
//			unitig=reverseComplement(unitig);
//		}
//		uint64_t overlapNum(getRepresentNum(unitig.substr(0,k-1)));
//		return (mapOnLeftEndMax(read, path, {overlapNum,overlap.second-(unitig.size()-k+1)}));
//	}
//	return false;
//}


//seem OK, exaustive
int Aligner::mapOnLeftEndExhaustive(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap , int errors){
	string unitig;
	string readLeft(read.substr(0,overlap.second));
	//	cout<<readLeft<<endl;
	vector<uNumber> path2keep;
	if(readLeft.size()==0){return 0;}
	auto rangeUnitigs(getEnd(overlap.first));
	int miniMiss(errors+1);
	int miniMissIndice(9);
	int offset(-2);
	bool ended(false);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()+k-1 <= unitig.size()){
			int miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				offset=unitig.size()-readLeft.size()-k+1;
				ended=true;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			int miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			if(miss<miniMiss){
				uint64_t overlapNum(str2num(unitig.substr(0,k-1)));
				vector<uNumber> possiblePath;
				miss+=mapOnLeftEndExhaustive(read , possiblePath, {overlapNum,overlap.second-(unitig.size()-k+1)}, errors-miss);
				if(miss<miniMiss){
					path2keep=possiblePath;
					miniMiss=miss;
					miniMissIndice=i;
					offset=-1;
					ended=false;
				}
			}
		}
	}

	if (miniMiss<=errors){
		if(ended){
			path.push_back(offset);
		}else{
			path.insert(path.end(), path2keep.begin(),path2keep.end());
		}
		path.push_back(rangeUnitigs[miniMissIndice].second);
	}

	return miniMiss;
}


int Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap , int errors){
	string unitig;

	string readLeft(read.substr(0,overlap.second));
	if(readLeft.empty()){path.push_back(0);return true;}
	auto rangeUnitigs(getEnd(overlap.first));
	int miniMiss(errors+1);
	int miniMissIndice(9);
	bool ended(false);
	int offset(0);
	uint64_t nextOverlap(0);
	string nextUnitig;


	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()+k-1 <= unitig.size()){
			int miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			int miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			if(miss<miniMiss){
				uint64_t overlapNum(str2num(unitig.substr(0,k-1)));
				//				miss+=mapOnLeftEnd2(read , path, {overlapNum,overlap.second-(unitig.size()-k+1)}, errors-miss);
				if(miss<miniMiss){
					ended=false;
					miniMiss=miss;
					miniMissIndice=i;
					nextUnitig=unitig;
					nextOverlap=overlapNum;
				}
			}
		}
	}

	if (miniMiss<=errors){
		if(ended){
			path.push_back(offset);
			path.push_back(rangeUnitigs[miniMissIndice].second);
			return miniMiss;
		}
		miniMiss+=mapOnLeftEndGreedy(read , path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, errors-miniMiss);
		path.push_back(rangeUnitigs[miniMissIndice].second);
	}

	return miniMiss;
}


int Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap , int errors){
	string unitig;
	string readLeft(read.substr(overlap.second));

	if(readLeft.size()<=k-1){return true;}

	auto rangeUnitigs(getBegin(overlap.first));
	int miniMiss(errors+1);
	int miniMissIndice(9);
	bool ended(false);
	uint64_t nextOverlap(0);
	string nextUnitig;

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()<= unitig.size()){
			int miss(missmatchNumber(unitig.substr(0,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
			}
			unitig=reverseComplement(unitig);
			miss=(missmatchNumber(unitig.substr(0,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			int miss(missmatchNumber(unitig, read.substr(overlap.second,unitig.size()), errors));
			if(miss<miniMiss){
				uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
				//				miss+=mapOnRightEnd2(read , path, {overlapNum,overlap.second+(unitig.size()-k+1)}, errors-miss);
				if(miss<miniMiss){
					miniMiss=miss;
					miniMissIndice=i;
					nextUnitig=unitig;
					nextOverlap=overlapNum;
				}
			}
			unitig=reverseComplement(unitig);
			miss=(missmatchNumber(unitig, read.substr(overlap.second,unitig.size()), errors));
			if(miss<miniMiss){
				uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
				//				miss+=mapOnRightEnd2(read , path, {overlapNum,overlap.second+(unitig.size()-k+1)}, errors-miss);
				if(miss<miniMiss){
					miniMiss=miss;
					miniMissIndice=i;
					nextUnitig=unitig;
					nextOverlap=overlapNum;
				}
			}
		}
	}

	if(miniMiss<errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if (ended){
			return miniMiss;
		}
		miniMiss+=mapOnRightEndGreedy(read , path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)}, errors-miniMiss);
	}

	return miniMiss;
}


//seem ok exaustive
int Aligner::mapOnRightEndExhaustive(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap , int errors){
	string unitig;
	string readLeft(read.substr(overlap.second+k-1));
	vector<uNumber> path2keep;


	if(readLeft.empty()){return 0;}
	auto rangeUnitigs(getBegin(overlap.first));
	int miniMiss(errors+1);
	int miniMissIndice(9);
	bool ended(false);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()<= unitig.size()-k+1){
			int miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
			}

		}else{
			//case the read is big enough we want to recover a true overlap
			int miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss<miniMiss){
				uint64_t overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				vector<uNumber> possiblePath;
				miss+=mapOnRightEndExhaustive(read , possiblePath, {overlapNum,overlap.second+(unitig.size()-k+1)}, errors-miss);
				if(miss<miniMiss){
					path2keep=possiblePath;
					miniMiss=miss;
					miniMissIndice=i;
					ended=false;
				}
			}
		}
	}

	if(miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(!ended){
			path.insert(path.end(), path2keep.begin(),path2keep.end());
		}
		successMR++;
	}

	return miniMiss;
}

//not exaustive
pair<size_t,int> Aligner::mapOnRight2(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap, const  vector<pair<uint64_t,uint>>& listOverlap, bool& ended,size_t start, int errors){
	//	cout<<"mOR2"<<endl;
	string unitig;
	string readLeft(read.substr(overlap.second+k-1));
	//	cout<<"rl : "<<readLeft<<endl;
	if(readLeft.empty()){cout<<"should not appears"<<endl;exit(0);return {start,0};}
	auto rangeUnitigs(getBegin(overlap.first));
	int miniMiss(errors+1);
	int miniMissIndice(9);
	size_t next(start);
	uint64_t nextOverlapNum(0);
	string nextUnitig;

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size() <= unitig.size()-k+1){
			int miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				ended=true;
				miniMiss=miss;
				miniMissIndice=i;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			int miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss<miniMiss){
				uint64_t overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				if(miss<miniMiss){
					//					cout<<"sub"<<unitig.substr(unitig.size()-k+1,k-1)<<endl;
					//					cout<<num2str(overlapNum)<<endl;
					ended=false;
					miniMiss=miss;
					miniMissIndice=i;
					nextOverlapNum=overlapNum;
					nextUnitig=unitig;
					next=start;
					for(uint j(start+1); j<listOverlap.size(); ++j){
						if(overlapNum==listOverlap[j].first and listOverlap[j].second==overlap.second+unitig.size()-k+1){
							//							cout<<listOverlap[j].second<<endl;
							//							cout<<num2str(listOverlap[j].first)<<endl;
							//							cout<<num2str(overlapNum)<<endl;
							next=j;
						}
					}
				}
			}
		}
	}
	if(ended){
		//		cout<<"end"<<endl;
		//		cout<<getUnitig(rangeUnitigs[miniMissIndice].second)<<endl;
		path.push_back(rangeUnitigs[miniMissIndice].second);
		return {start,miniMiss};
	}

	if(miniMiss<=errors){
		//		cout<<getUnitig(rangeUnitigs[miniMissIndice].second)<<endl;;
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(next>start){
			return {next,miniMiss};
		}
		auto res(mapOnRight2(read , path, {nextOverlapNum,overlap.second+(nextUnitig.size()-k+1)},listOverlap,ended,start, errors-miniMiss));
		return {res.first,res.second+miniMiss};
	}

	return {start,errors+1};
}



//bool Aligner::mapOnRightEndMax(const string &read, vector<uint32_t>& path, pair<uint64_t, uint> overlap){
//	string unitig;
//	string readLeft(read.substr(overlap.second));
//	if(readLeft.size()<=k-1){return true;}
//	vector<uint32_t> rangeUnitigs(getBegin(read.substr(overlap.second,k-1)));
//	double minimum(100);
//	int minimumindice(9);
//
//	for(uint i(0); i<rangeUnitigs.size(); ++i){
//		unitig=getUnitig(rangeUnitigs[i]);
//		//case the rest of the read is too small
//		if(readLeft.size() <= unitig.size()){
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				return true;
//			}
//		}else{
//			//case the read is big enough we want to recover a true overlap
//			double miss(missMatchPercent(unitig, read.substr(overlap.second,unitig.size())));
//			if(miss<minimum){
//				minimum=miss;
//				minimumindice=2*i;
//			}
//			unitig=reverseComplement(unitig);
//			miss=(missMatchPercent(unitig, read.substr(overlap.second,unitig.size())));
//			if(miss<minimum){
//				minimum=miss;
//				minimumindice=2*i+1;
//			}
//		}
//	}
//
//	if(minimum<errorRate){
//		unitig=getUnitig(rangeUnitigs[minimumindice/2]);
//		if(minimumindice%2==1){
//			unitig=reverseComplement(unitig);
//		}
//		uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//		return (mapOnRightEndMax(read, path, {overlapNum,overlap.second+unitig.size()-k+1}));
//	}
//
//	return false;
//}


//size_t Aligner::mapOnRight(const string &read, vector<uint32_t>& path, pair<uint64_t, uint> overlap, const  vector<pair<uint64_t,uint>>& listOverlap, bool& ended,size_t start){
//	//	cout<<"mor"<<endl;
//	vector<uint32_t> rangeUnitigs(getEnd(read.substr(overlap.second,k-1)));
//	string unitig;
//	string readLeft(read.substr(overlap.second));
//	if(readLeft.size()<=k-1){return true;}
//
//	for(uint i(0); i<rangeUnitigs.size(); ++i){
//		unitig=getUnitig(rangeUnitigs[i]);
//		//case the rest of the read is too small
//		if(readLeft.size() <= unitig.size()){
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				ended=true;
//				return 0;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				ended=true;
//				return 0;
//			}
//		}else{
//			//case the read is big enough we want to recover a true overlap
//			if(missmatch(unitig, read.substr(overlap.second,unitig.size()), (unsigned char)(unitig.size()*errorRate))){
//				uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//				for(uint j(start+1); j<listOverlap.size(); ++j){
//					if(overlapNum==listOverlap[j].first){
//						return j;
//					}
//				}
//				return mapOnRight(read, path, {overlapNum,overlap.second+unitig.size()-k+1}, listOverlap,ended,start);
//
//			}else{
//				unitig=reverseComplement(unitig);
//				if(missmatch(unitig, read.substr(overlap.second,unitig.size()), (unsigned char)(unitig.size()*errorRate))){
//					uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//					for(uint j(start+1); j<listOverlap.size(); ++j){
//						if(overlapNum==listOverlap[j].first){
//							return j;
//						}
//					}
//					return mapOnRight(read, path, {overlapNum,overlap.second+unitig.size()-k+1}, listOverlap,ended,start);
//				}
//			}
//		}
//	}
//	//mean fail
//	return 0;
//}


//size_t Aligner::mapOnRightAll(const string &read, vector<uint32_t>& path, pair<uint64_t, uint> overlap, const  vector<pair<uint64_t,uint>>& listOverlap, bool& ended,size_t start){
//	auto rangeUnitigs(overlap2unitigs.find(overlap.first)->second);
//	string unitig;
//	string readLeft(read.substr(overlap.second));
//	if(readLeft.size()<=k-1){return true;}
//
//	for(uint i(0); i<rangeUnitigs.size(); ++i){
//		unitig=getUnitig(rangeUnitigs[i]);
//		//case the rest of the read is too small
//		if(readLeft.size() <= unitig.size()){
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				ended=true;
//				return 0;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(0,readLeft.size()), readLeft, (unsigned char)((readLeft.size())*errorRate))){
//				ended=true;
//				return 0;
//			}
//		}else{
//			//case the read is big enough we want to recover a true overlap
//			if(missmatch(unitig, read.substr(overlap.second,unitig.size()), 1+(unsigned char)(unitig.size()*errorRate))){
//				uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//				for(uint j(start+1); j<listOverlap.size(); ++j){
//					if(overlapNum==listOverlap[j].first){
//						return j;
//					}
//				}
//				size_t mora(mapOnRightAll(read, path, {overlapNum,overlap.second+unitig.size()-k+1}, listOverlap,ended,start));
//				if(mora>start or ended){
//					return mora;
//				}
//
//			}else{
//				unitig=reverseComplement(unitig);
//				if(missmatch(unitig, read.substr(overlap.second,unitig.size()), 1+(unsigned char)(unitig.size()*errorRate))){
//					uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//					for(uint j(start+1); j<listOverlap.size(); ++j){
//						if(overlapNum==listOverlap[j].first){
//							return j;
//						}
//					}
//					size_t mora(mapOnRightAll(read, path, {overlapNum,overlap.second+unitig.size()-k+1}, listOverlap,ended,start));
//					if(mora>start or ended){
//						return mora;
//					}
//				}
//			}
//		}
//	}
//	//mean fail
//	return 0;
//}


string Aligner::getUnitig2(uint position){
	string unitig;
	unitigMutex.lock();
	{
		unitigFile.seekg(position,ios::beg);
		getline(unitigFile,unitig,';');
	}
	unitigMutex.unlock();
	transform(unitig.begin(), unitig.end(), unitig.begin(), ::toupper);
	if(unitig.size()<k){
		cout<<"bug"<<endl;
		exit(0);
	}
	return unitig;
}


string Aligner::readUnitig(uint position){
	string unitig;
	unitigMutex.lock();
	{
		unitigFile.seekg(position,ios::beg);
		getline(unitigFile,unitig,';');
	}
	unitigMutex.unlock();
	transform(unitig.begin(), unitig.end(), unitig.begin(), ::toupper);
	if(unitig.size()<k){
		cout<<"bug"<<endl;
		exit(0);
	}
	return unitig;
}


string Aligner::getUnitig(int position){
	if(fullMemory){
		if(position>0){
			return unitigs[position];
		}
		return reverseComplement(unitigs[-position]);
	}
	string unitig;
	unitigMutex.lock();
	{
		unitigFile.seekg(position,ios::beg);
		getline(unitigFile,unitig,';');
	}
	unitigMutex.unlock();
	transform(unitig.begin(), unitig.end(), unitig.begin(), ::toupper);
	if(unitig.size()<k){
		cout<<"bug"<<endl;
		exit(0);
	}
	return unitig;
}


//bool Aligner::checkPair(const pair<uint64_t, uint>& overlap1, const pair<uint64_t, uint>& overlap2, const string& read, vector<uint32_t>& path){
//	auto rangeunitigs1 = overlap2unitigs.find(overlap1.first)->second;
//	auto rangeunitigs2 = overlap2unitigs.find(overlap2.first)->second;
//	string unitig,readLeft(read.substr(overlap1.second,overlap2.second-overlap1.second+k-1));
//	uint32_t positionget1, positionget2;
//	for(uint i(0); i<rangeunitigs1.size(); ++i){
//		positionget1=rangeunitigs1[i];
//		for(uint j(0); j<rangeunitigs2.size(); ++j){
//			positionget2=rangeunitigs2[j];
//			if(positionget2==positionget1){
//				unitig=getUnitig(positionget1);
//				if(unitig.size()==readLeft.size()){
//					if(missmatch(unitig, readLeft, (unsigned char)(unitig.size()*errorRate))){
//						path.push_back(positionget1);
//						return true;
//					}
//					unitig=reverseComplement(unitig);
//					if(missmatch(unitig, readLeft, (unsigned char)(unitig.size()*errorRate))){
//						path.push_back(positionget1);
//						return true;
//					}
//				}
//			}
//		}
//	}
//	return false;
//}


//seem ok, exaustive
int Aligner::checkPair2(const pair<uint64_t, uint>& overlap1, const pair<uint64_t, uint>& overlap2, const string& read, uNumber& number, int errorsAllowed){
	//	cout<<overlap1.second<<" "<<overlap2.second<<endl;
	if(overlap2.second-overlap1.second<k){
		int32_t positionget1, positionget2;
		auto rangeUnitigs1(getBegin(overlap1.first));
		auto rangeUnitigs2(getEnd(overlap2.first));
		for(uint i(0); i<rangeUnitigs1.size(); ++i){
			positionget1=rangeUnitigs1[i].second;
			for(uint j(0); j<rangeUnitigs2.size(); ++j){
				positionget2=rangeUnitigs2[j].second;
				if(positionget2==positionget1 and rangeUnitigs1[i].first.size()==(overlap2.second-overlap1.second+k-1)){
					number=positionget1;
					//					cout<<getUnitig(number)<<endl;
					return 0;
				}
			}
		}
		return errorsAllowed+1;
	}
	string unitig,subRead(read.substr(overlap1.second+k-1,overlap2.second-overlap1.second-(k-1)));
	auto rangeUnitigs1(getBegin(overlap1.first));
	auto rangeUnitigs2(getEnd(overlap2.first));
	int minMissMatch(errorsAllowed+1);
	int indice(0);
	int32_t positionget1, positionget2;
	for(uint i(0); i<rangeUnitigs1.size(); ++i){
		positionget1=rangeUnitigs1[i].second;
		for(uint j(0); j<rangeUnitigs2.size(); ++j){
			positionget2=rangeUnitigs2[j].second;
			if(positionget2==positionget1){
				unitig=getUnitig(rangeUnitigs1[i].second);
				if(unitig.size()-2*(k-1)==subRead.size()){
					int missmatch(missmatchNumber(unitig.substr(k-1,subRead.size()), subRead, errorsAllowed));
					if(missmatch<minMissMatch){
						minMissMatch=missmatch;
						indice=i;
					}
				}
			}
		}
	}

	if(minMissMatch<=errorsAllowed){
		//		cout<<getUnitig(number)<<endl;;
		number=(rangeUnitigs1[indice].second);
	}
	return minMissMatch;
}


//bool Aligner::checkBegin(const string& read, pair<uint64_t, uint>& overlap, vector<uint32_t>& path){
//	if(overlap.second==0){return true;}
//	auto rangeunitigs = overlap2unitigs.find(overlap.first)->second;
//	string readLeft(read.substr(0,overlap.second));
//	string unitig;
//
//	for(uint i(0);i<rangeunitigs.size();++i){
//		unitig=getUnitig(rangeunitigs[i]);
//		if(unitig.size()>=overlap.second+k-1){
//			if(missmatch(unitig.substr(unitig.size()-overlap.second-k+1,overlap.second),readLeft,(unsigned int)(errorRate*overlap.second))){
//				path.push_back(rangeunitigs[i]);
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(unitig.size()-overlap.second-k+1,overlap.second),readLeft,(unsigned int)(overlap.second*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				return true;
//			}
//		}
//	}
//	return false;
//}
//
//
//bool Aligner::checkBegin2(const string& read, pair<uint64_t, uint>& overlap, vector<uint32_t>& path){
//	if(overlap.second==0){return true;}
//	string readLeft(read.substr(0,overlap.second));
//	auto rangeunitigs = overlap2unitigs.find(overlap.first)->second;
//	string unitig;
//
//	for(uint i(0); i<rangeunitigs.size(); ++i){
//		unitig=getUnitig(rangeunitigs[i]);
//		if(unitig.size()>=readLeft.size()+k-1){
//			if(missmatch(unitig.substr(unitig.size()-overlap.second-k+1,overlap.second),readLeft, (unsigned int)(errorRate*readLeft.size()))){
//				path.push_back(rangeunitigs[i]);
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(unitig.size()-overlap.second-k+1,overlap.second),readLeft, (unsigned int)(readLeft.size()*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				return true;
//			}
//		}else{
//			if(missmatch(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), 1+(unsigned int)((unitig.size()-k+1)*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				uint64_t overlapNum(getRepresentNum(unitig.substr(0,k-1)));
//				if(mapOnLeftEndMax(read, path, {overlapNum,overlap.second-(unitig.size()-k+1)})){
//					sucessML++;
//					return true;
//				}
//				//				return mapOnLeftEndFirst(read, path, {overlapNum,overlap.second-(unitig.size()-k+1)});
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(0,unitig.size()-k+1),readLeft.substr(readLeft.size()+k-1-unitig.size()), 1+(unsigned int)((unitig.size()-k+1)*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				uint64_t overlapNum(getRepresentNum(unitig.substr(0,k-1)));
//				if(mapOnLeftEndMax(read, path, {overlapNum,overlap.second-(unitig.size()-k+1)})){
//					sucessML++;
//					return true;
//				}
//				//				return mapOnLeftEndFirst(read, path, {overlapNum,overlap.second-(unitig.size()-k+1)});
//			}
//		}
//	}
//	return false;
//}

//Seem ok;exhaustive;
int Aligner::checkBeginExhaustive(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, int errors){

	if(overlap.second==0){
		path.push_back(0);
		return 0;
	}

	string readLeft(read.substr(0,overlap.second));
	//	auto rangeunitigs(overlap2unitigs.find(overlap.first)->second);
	auto rangeUnitigs(getEnd(overlap.first));
	string unitig;
	int minMiss(errors+1);
	int indiceMinMiss(0);
	int offset(-2);
	bool ended(false);
	vector<uNumber> path2keep;
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);

		if(unitig.size()-k+1>=readLeft.size()){
			int miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=(unitig.size()-readLeft.size()-k+1);
			}
		}else{
			int miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			if(miss<minMiss){
				uint64_t overlapNum(str2num(unitig.substr(0,k-1)));
				vector<uNumber> possiblePath;
				miss+=mapOnLeftEndExhaustive(read, possiblePath, {overlapNum,overlap.second-(unitig.size()-k+1)},errors-miss);
				if(miss<minMiss){
					sucessML++;
					minMiss=miss;
					indiceMinMiss=i;
					path2keep=possiblePath;
					ended=false;
				}
			}
		}
	}

	if(minMiss<=errors){
		if(ended){
			path.push_back(offset);
			path.push_back(rangeUnitigs[indiceMinMiss].second);
		}else{
			path.insert(path.end(), path2keep.begin(), path2keep.end());
			path.push_back(rangeUnitigs[indiceMinMiss].second);
		}
	}

	return minMiss;
}


int Aligner::checkBeginGreedy(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, int errors){
	if(overlap.second==0){path.push_back(0);return 0;}

	string readLeft(read.substr(0,overlap.second));
	//	auto rangeunitigs(overlap2unitigs.find(overlap.first)->second);
	auto rangeUnitigs(getEnd(overlap.first));
	string unitig;
	int minMiss(errors+1);
	int indiceMinMiss(0);
	bool ended(false);
	int offset(0);
	uint64_t nextOverlap(0);
	string nextUnitig;

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		if(unitig.size()-k+1>=readLeft.size()){
			int miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			int miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			if(miss<minMiss){
				uint64_t overlapNum(str2num(unitig.substr(0,k-1)));
				if(miss<minMiss){
					ended=false;
					minMiss=miss;
					indiceMinMiss=i;
					nextOverlap=overlapNum;
					nextUnitig=unitig;
				}
			}

		}
	}


	if(minMiss<=errors){
		if(ended){
			path.push_back(offset);
			path.push_back(rangeUnitigs[indiceMinMiss].second);
			return minMiss;
		}
		minMiss+=mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)},errors-minMiss);
		if(minMiss<=errors){
			path.push_back(rangeUnitigs[indiceMinMiss].second);
			sucessML++;
		}
	}

	return minMiss;
}



//seem ok; exaustive
int Aligner::checkEndExhaustive(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, int errors){
	string readLeft(read.substr(overlap.second+k-1));
	vector<uNumber> path2keep;
	if(readLeft.empty()){
		return 0;
	}
	//	auto rangeunitigs = overlap2unitigs.find(overlap.first)->second;
	auto rangeUnitigs(getBegin(overlap.first));
	string unitig;
	int minMiss(errors+1);
	int indiceMinMiss(9);
	bool ended(false);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		if(unitig.size()-k+1>=readLeft.size()){
			int miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
			}
		}else{
			int miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss<minMiss){
				uint64_t overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				vector<uNumber> possiblePath;
				miss+=mapOnRightEndExhaustive(read, possiblePath, {overlapNum,overlap.second+(unitig.size()-k+1)},errors-miss);
				if(miss<minMiss){
					path2keep=possiblePath;
					minMiss=miss;
					indiceMinMiss=i;
					ended=false;
				}
			}
		}
	}
	if(minMiss<=errors){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(!ended){
			path.insert(path.end(), path2keep.begin(),path2keep.end());
		}
	}

	return minMiss;
}


//seem ok; exaustive
int Aligner::checkEndGreedy(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, int errors){
	string readLeft(read.substr(overlap.second+k-1));

	if(readLeft.size()<=k-1){return 0;}

	//	auto rangeunitigs = overlap2unitigs.find(overlap.first)->second;
	auto rangeUnitigs(getBegin(overlap.first));
	string unitig;
	int minMiss(errors+1);
	int indiceMinMiss(9);
	bool ended(false);
	uint64_t nextOverlap(0);
	string nextUnitig;

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);

		if(unitig.size()-k+1>=readLeft.size()){
			int miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
			}
		}else{
			int miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss<minMiss){
				uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
				//				miss+=mapOnRightEndExhaustive(read, path, {overlapNum,overlap.second+(unitig.size()-k+1)},errors-miss);
				if(miss<minMiss){
					minMiss=miss;
					indiceMinMiss=i;
					nextOverlap=overlapNum;
					nextUnitig=unitig;
				}
			}
		}
	}
	if(minMiss<=errors){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){
			return minMiss;
		}
		minMiss+=mapOnRightEndGreedy(read, path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)},errors-minMiss);
		if(minMiss<=errors){
			successMR++;
		}
	}

	return minMiss;
}

//
//bool Aligner::checkEnd(const string& read, pair<uint64_t, uint>& overlap, vector<uint32_t>& path){
//	if(read.size()<=overlap.second+k-1){return true;}
////	auto rangeunitigs = overlap2unitigs.find(overlap.first)->second;
//	vector<uint32_t> rangeUnitigs(getEnd(read.substr(overlap.second,k-1)));
//	string unitig;
//
//	for(uint i(0);i<rangeunitigs.size();++i){
//		unitig=getUnitig(rangeunitigs[i]);
//		if(unitig.size()>=read.size()-overlap.second){
//			if(missmatch(unitig.substr(0,read.size()-overlap.second),read.substr(overlap.second),(unsigned int)(overlap.second*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(0,read.size()-overlap.second),read.substr(overlap.second),(unsigned int)(overlap.second*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				return true;
//			}
//		}
//	}
//	return false;
//}

//
//bool Aligner::checkEnd2(const string& read, pair<uint64_t, uint>& overlap, vector<uint32_t>& path){
//	if(read.size()<=overlap.second+k-1){return true;}
//	auto rangeunitigs = overlap2unitigs.find(overlap.first)->second;
//	string unitig;
//
//	for(uint i(0);i<rangeunitigs.size();++i){
//		unitig=getUnitig(rangeunitigs[i]);
//		if(unitig.size()>=read.size()-overlap.second){
//			if(missmatch(unitig.substr(0,read.size()-overlap.second),read.substr(overlap.second),(unsigned int)(overlap.second*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				return true;
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig.substr(0,read.size()-overlap.second),read.substr(overlap.second),(unsigned int)(overlap.second*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				return true;
//			}
//		}else{
//			if(missmatch(unitig,read.substr(overlap.second,unitig.size()),1+(unsigned int)(unitig.size()*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//				if(mapOnRightEndMax(read, path, {overlapNum,overlap.second+unitig.size()-k+1})){
//					successMR++;
//					return true;
//				}
//			}
//			unitig=reverseComplement(unitig);
//			if(missmatch(unitig,read.substr(overlap.second,unitig.size()),1+(unsigned int)(unitig.size()*errorRate))){
//				path.push_back(rangeunitigs[i]);
//				uint64_t overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
//				if(mapOnRightEndMax(read, path, {overlapNum,overlap.second+unitig.size()-k+1})){
//					successMR++;
//					return true;
//				}
//			}
//		}
//	}
//
//	return false;
//}

//Seem correct indices checked
//vector<pair<uint64_t,uint>> Aligner::getListOverlap(const string& read){
//	vector<pair<uint64_t,uint>> listOverlap;
//	string overlap;
//	uint64_t num,bin;
//
//
//	for(uint i(0);(i+k-1)<=read.size();++i){
//		overlap=read.substr(i,k-1);
//		bin=str2num(overlap);
//
//		num=getRepresentNum(overlap);
//		if(left.count(num)!=0){
//			listOverlap.push_back({bin,i});
//		}else{
//			if(right.count(num)!=0){
//				listOverlap.push_back({bin,i});
//			}
//		}
//	}
//	return listOverlap;
//}


void Aligner::update(uint64_t&	min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdate;
}


void Aligner::updateRC(uint64_t&	min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-4));
	//	cout<<k<<endl;
	//	cout<<(2*k-2)<<endl;
	//	cout<<((nuc2intrc(nuc))<<(2))<<endl;
	//	cout<<"        "<<num2str(((nuc2intrc(nuc))<<(2*k-4)));
	//	cin.get();
}

vector<pair<uint64_t,uint>> Aligner::getListOverlap(const string& read){
	//	cout<<read<<endl;
	vector<pair<uint64_t,uint>> listOverlap;
	string overlap(read.substr(0,k-1));
	uint64_t num(str2num(overlap)),rcnum(str2num(reverseComplement(overlap)));
	uint64_t rep(min(num, rcnum));


	for(uint i(0);;++i){
		if(left.count(rep)!=0){
			listOverlap.push_back({num,i});
		}else{
			if(right.count(rep)!=0){
				listOverlap.push_back({num,i});
			}
		}
		if(i+k-1<read.size()){
			//			cout<<read[i+k-1]<<endl;;
			update(num,read[i+k-1]);
			updateRC(rcnum,read[i+k-1]);
			rep=(min(num, rcnum));
			//			if(num!=str2num(read.substr(i+1,k-1))){
			//				cout<<read.substr(i+1,k-1)<<endl;
			//				cout<<num2str(num)<<endl;;
			//				cout<<"bug2"<<endl;
			//				exit(0);
			//			}
			//			if(rcnum!=str2num((reverseComplement(read.substr(i+1,k-1))))){
			//				cout<<"bug1"<<endl;
			//				cout<<(reverseComplement(read.substr(i+1,k-1)))<<endl;
			//				cout<<num2str(rcnum)<<endl;
			//				exit(0);
			//			}
		}else{
			//			cin.get();
			return listOverlap;
		}
		//		cout<<"overlap "<<num2str(num)<<endl;
		//		cout<<"rc      "<<reverseComplement(num2str(num))<<endl;;
		//		cout<<"rctry   "<<num2str(rcnum)<<endl;
	}
	return listOverlap;
}


void Aligner::alignPartExhaustive(){
	vector<string> multiread;
	vector<uNumber> path;
	string read;
	bool overlapFound(false);
	while(!readFile.eof()){
		readMutex.lock();
		{
			getReads(multiread,10000);
		}
		readMutex.unlock();
		for(size_t i(0);i<multiread.size();++i){
			read=multiread[i],0;
			overlapFound=false;
			path=alignRead4(read,overlapFound,errorsMax);
			if(path.size()!=0){
				//				pathMutex.lock();
				//				{
				//					printPath(path,&pathFile);
				//				}
				//				pathMutex.unlock();
			}else{
				//				noOverlapMutex.lock();
				notMappedFile<<">A"<<endl<<read<<endl;
				//				noOverlapMutex.unlock();
				if(!overlapFound){
					//					noOverlapMutex.lock();
					//					{
					//						noOverlapFile<<">A"<<endl<<read<<endl;
					//					}
					//					noOverlapMutex.unlock();
				}else{
					//					notMappedMutex.lock();
					//					{
					//						notMappedFile<<">A"<<endl<<read<<endl;
					//					}
					//					notMappedMutex.unlock();
				}
			}
		}

		cout<<"Read : "<<readNumber<<endl;
		cout<<"No Overlap : "<<noOverlapRead<<" Percent : "<<(100*float(noOverlapRead))/readNumber<<endl;
		cout<<"Got Overlap : "<<alignedRead+notAligned<<" Percent : "<<(100*float(alignedRead+notAligned))/readNumber<<endl;
		cout<<"Overlap and Aligned : "<<alignedRead<<" Percent : "<<(100*float(alignedRead))/(alignedRead+notAligned)<<endl;
		cout<<"Overlap but no aligne: "<<notAligned<<" Percent : "<<(100*float(notAligned))/(alignedRead+notAligned)<<endl;
		auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
		cout<<"Reads/seconds : "<<readNumber/(chrono::duration_cast<chrono::seconds>(waitedFor).count()+1)<<endl;
		cout<<"Overlap per reads : "<<(overlaps)/(alignedRead+notAligned)<<endl;

		cout<<endl;
	}
}


void Aligner::alignPartGreedy(){
	vector<string> multiread;
	vector<uNumber> path;
	string read;
	bool overlapFound(false);
	while(!readFile.eof()){
		readMutex.lock();
		{
			getReads(multiread,1000000);
		}
		readMutex.unlock();
		for(size_t i(0);i<multiread.size();++i){
			overlapFound=false;
			read=multiread[i];
			path=alignRead3(read,overlapFound,errorsMax,false);
			if(path.size()!=0){
				//				pathMutex.lock();
				//				{
				//					printPath(path,&pathFile);
				//				}
				//				pathMutex.unlock();
			}else{
				//				noOverlapMutex.lock();
				notMappedFile<<">A"<<endl<<read<<endl;
				//				noOverlapMutex.unlock();
				if(!overlapFound){
					//					noOverlapMutex.lock();
					//					{
					//						noOverlapFile<<">A"<<endl<<read<<endl;
					//					}
					//					noOverlapMutex.unlock();
				}else{
					//					notMappedMutex.lock();
					//					{
					//						notMappedFile<<">A"<<endl<<read<<endl;
					//					}
					//					notMappedMutex.unlock();
				}
			}
		}
		cout<<"Read : "<<readNumber<<endl;
		cout<<"No Overlap : "<<noOverlapRead<<" Percent : "<<(100*float(noOverlapRead))/readNumber<<endl;
		cout<<"Got Overlap : "<<alignedRead+notAligned<<" Percent : "<<(100*float(alignedRead+notAligned))/readNumber<<endl;
		cout<<"Overlap and Aligned : "<<alignedRead<<" Percent : "<<(100*float(alignedRead))/(alignedRead+notAligned)<<endl;
		cout<<"Overlap but no aligne: "<<notAligned<<" Percent : "<<(100*float(notAligned))/(alignedRead+notAligned)<<endl;
		auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
		cout<<"Reads/seconds : "<<readNumber/(chrono::duration_cast<chrono::seconds>(waitedFor).count()+1)<<endl;
		cout<<"Overlap per reads : "<<(overlaps)/(alignedRead+notAligned)<<endl;
		cout<<endl;
	}
}

//
//void Aligner::indexUnitigsAux2(){
//	uint32_t position;
//	string line, overlap1, overlap2, waste;
//	uint64_t size, num1, num2;
//	while(!unitigFile.eof()){
////		if(unitigNumber%10000==0){
////			cout<<"Number of unitig k : "<<unitigNumber/1000<<endl;
////		}
////		unitigMutex.lock();
//		position=unitigFile.tellg();
//		getline(unitigFile,line,';');
//		unitigFile.seekg(1, ios::cur);
////		unitigMutex.unlock();
//		transform(line.begin(), line.end(), line.begin(), ::toupper);
//		size=line.size();
//		if(size<k){
//			break;
//		}else{
//			unitigNumber++;
//			num1=getRepresentNum(line.substr(0,k-1));
//			num2=getRepresentNum(line.substr(size-k+1,k-1));
//			if(fullMemory){
////				indexMutex.lock();
//				{
//					unitigCache[position]=line;
//				}
////				indexMutex.unlock();
//			}
////			indexMutex.lock();
//			{
//				overlap2unitigs[num1].push_back(position);
//				if(num1!=num2){
//					overlap2unitigs[num2].push_back(position);
//				}
//			}
////			indexMutex.unlock();
//		}
//	}
//}

void Aligner::indexUnitigsAux(){
	uint32_t position(0);
	string line, overlap1, overlap2, waste;
	uint64_t size;
	unitigs.push_back("");
	while(!unitigFile.eof()){
		if(!fullMemory){
			position=unitigFile.tellg();
		}
		getline(unitigFile,line,';');
		unitigFile.seekg(1, ios::cur);
		transform(line.begin(), line.end(), line.begin(), ::toupper);
		size=line.size();
		if(size<k){
			break;
		}else{
			if(fullMemory){
				unitigs.push_back(line);
				position=unitigs.size()-1;
			}
			unitigNumber++;
			string beg(line.substr(0,k-1));
			string rcBeg(reverseComplement(beg));
			if(beg<=rcBeg){
				left[str2num(beg)].push_back(position);
			}else{
				right[str2num(rcBeg)].push_back(position);
			}
			string end(line.substr(size-k+1,k-1));
			string rcEnd(reverseComplement(end));
			if(end<=rcEnd){
				right[str2num(end)].push_back(position);
			}else{
				left[str2num(rcEnd)].push_back(position);
			}
		}
	}
}


void Aligner::indexUnitigs(){
	unsigned char nbThreads(1);
	vector<thread> threads;
	for (size_t i(0); i<nbThreads; ++i){
		threads.push_back(thread(&Aligner::indexUnitigsAux,this));
	}
	for(auto &t : threads){t.join();}
	cout<<"Number of unitig: "<<unitigNumber<<endl;
	unitigFile.close();
	//	ofstream out(unitigFileName, ios::app);
	//	out<<':';
	//	out.close();
	unitigFile.open(unitigFileName);
}


void Aligner::alignAll(bool greedy){
	startChrono=chrono::system_clock::now();
	unsigned char nbThreads(coreNumber);
	vector<thread> threads;
	for (size_t i(0); i<nbThreads; ++i){
		if(greedy){
			threads.push_back(thread(&Aligner::alignPartGreedy,this));
		}else{
			threads.push_back(thread(&Aligner::alignPartExhaustive,this));
		}

	}

	for(auto &t : threads){t.join();}
	cout<<"Read : "<<readNumber<<endl;
	cout<<"No Overlap : "<<noOverlapRead<<" Percent : "<<(100*float(noOverlapRead))/readNumber<<endl;
	cout<<"Got Overlap : "<<alignedRead+notAligned<<" Percent : "<<(100*float(alignedRead+notAligned))/readNumber<<endl;
	cout<<"Overlap and Aligned : "<<alignedRead<<" Percent : "<<(100*float(alignedRead))/(alignedRead+notAligned)<<endl;
	cout<<"Overlap but no aligne: "<<notAligned<<" Percent : "<<(100*float(notAligned))/(alignedRead+notAligned)<<endl;
	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
	cout<<"Reads/seconds : "<<readNumber/(chrono::duration_cast<chrono::seconds>(waitedFor).count()+1)<<endl;
	cout<<"Overlap per reads : "<<(overlaps)/readNumber<<endl;
	cout<<endl;
}
