/*****************************************************************************
 *   Bgreat : De Bruijn Graph Read Mapping Tool
 *   Copyright (C) 2014  INRIA
 *   Authors: Antoine Limasset
 *   Contact: antoine.limasset@inria.fr, INRIA/IRISA/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
 *   Source: https://github.com/Malfoy/BGREAT
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/



#include "aligner.h"
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


uint64_t transform_to_size_t(__uint128_t n){
	return n>>64;
}


void printPath(const vector<int32_t>& path, ofstream* file){
	for(size_t i(0); i<path.size(); ++i){
		*file<<path[i]<<'.';
	}
	*file<<endl;
}


char revCompChar(char s) {
	if (s == 'A') return 'T';
	else if (s == 'C') return 'G';
	else if (s == 'G') return 'C';
	else if (s == 'T') return 'A';
	// else if (s == 'a') return 't';
	// else if (s == 'c') return 'g';
	// else if (s == 'g') return 'c';
	// else if (s == 't') return 'a';
	//~ cout<<"wtf rc"<<s<<endl;
	//~ cin.get();
	return 'X';
}


string reverseComplement(const string& s){
	string rc(s.size(),0);
	for (int i = (int)s.length() - 1; i >= 0; i--){
		rc[s.size()-1-i]= revCompChar(s[i]);
	}
	return rc;
}


//Parse reads
void Aligner::getReads(vector<pair<string,string>>& reads, uint n){
	reads={};
	string read,header,inter;
	char c;
	if(fastq){
		for(uint i(0);i<n;++i){
			getline(readFile,header);
			getline(readFile,read);
			if(read.size()>2){
				bool fail(false);
				for(uint j(0);(j)<read.size();++j){
					if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G'){
						fail=true;
						break;
					}
				}
				if(!fail){
					reads.push_back({header,read});
				}
			}
			getline(readFile,header);
			getline(readFile,header);
			if(readFile.eof()){return;}
		}
	}else{

		for(uint i(0);i<n;++i){
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
						reads.push_back({header,read});
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
							reads.push_back({header,read});
						}
					}
					return;
				}
			}
		}
	}
}


kmer Aligner::getRepresentNum(const string& str){
	string rc(reverseComplement(str));
	kmer a(str2num(str)),b(str2num(rc));
	return ((b<=a) ? b : a);
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


string mutate(string& read, int n){
	for(int i(0); i<n; ++i){
		int position(rand()%read.size());
		read[position]=randNuc();
	}
	return read;
}


string Aligner::num2str(kmer num){
	string str;
	int nuc;
	for(uint i(0);i<k-1;++i){
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


kmer Aligner::str2num(const string& str){
	kmer res(0);
	for(uint i(0);i<str.size();i++){
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
		//~ case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	//~ cout<<"failnuc2int"<<endl;
	//~ cout<<c<<endl;
	//~ cin.get();
	return 0;
}


char nuc2intrc(char c){
	switch(c){
		case 'A': return 3;
		case 'C': return 2;
		case 'G': return 1;
		//~ case 'T': return 0;
	}
	//~ cout<<"failnuc2intRC"<<endl;
	//~ cout<<c<<endl;
	//~ cin.get();
	return 0;
}


vector<pair<string,uNumber>> Aligner::getEnd(kmer bin){
	vector<pair<string,uNumber>> result;
	string str(num2str(bin)),rc(reverseComplement(str)),unitig;
	if(str<=rc){
		if(right.count(bin)!=0){
			for(uint i(0); i<right[bin].size(); ++i){
				unitig=(getUnitig(right.at(bin)[i]));
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
		kmer num(str2num(rc));
		if(left.count(num)!=0){
			for(uint i(0); i<left[num].size(); ++i){
				unitig=(getUnitig(left.at(num)[i]));
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


vector<pair<string,uNumber>> Aligner::getBegin(kmer bin){
	vector<pair<string,uNumber>> result;
	string str(num2str(bin)),rc(reverseComplement(str)),unitig;
	if(str<=rc){
		if(left.count(bin)!=0){
			for(uint i(0); i<left[bin].size(); ++i){
				unitig=(getUnitig(left.at(bin)[i]));
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
		kmer num(str2num(rc));
		if(right.count(num)!=0){
			for(uint i(0); i<right[num].size(); ++i){
				unitig=(getUnitig(right.at(num)[i]));
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


vector<uNumber> Aligner::getEndNumber(kmer bin){
	vector<uNumber> result;
	string str(num2str(bin)),rc(reverseComplement(str));
	if(str<=rc){
		if(right.count(bin)!=0){
			for(uint i(0); i<right[bin].size(); ++i){
				result.push_back(right.at(bin)[i]);
			}
		}else{
			return {};
		}
	}else{
		kmer num(str2num(rc));
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


vector<uNumber> Aligner::getBeginNumber(kmer bin){
	vector<uNumber> result;
	string str(num2str(bin)), rc(reverseComplement(str));
	if(str<=rc){
		if(left.count(bin)!=0){
			for(uint i(0); i<left[bin].size(); ++i){
				result.push_back(left.at(bin)[i]);
			}
		}else{
			return {};
		}
	}else{
		kmer num(str2num(rc));
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


uint8_t missmatchNumber(const string& seq1, const string& seq2, unsigned int n){
	uint8_t miss(0);
	for(size_t i(0); i<seq2.size(); ++i){
		if(seq2[i]!=seq1[i]){
			if(++miss>n){
				return miss;
			}
		}
	}
	return miss;
}


vector<uNumber> Aligner::alignReadGreedy(const string& read, bool& overlapFound, uint8_t errors, bool rc){
	vector<pair<kmer,uint>> listOverlap(getListOverlap(read));
	if(listOverlap.empty()){
		noOverlapRead++;
		readNumber++;
		return {};
	}

	overlaps+=listOverlap.size();
	overlapFound=true;
	for(uint start(0); start<=min(tryNumber-1,(uint)listOverlap.size()-1); ++start){
		vector<uNumber> pathBegin;
		uint8_t errorBegin(checkBeginExhaustive(read,listOverlap[start],pathBegin,errors));
		if(errorBegin<=errors){
			for(int end((int)listOverlap.size()-1); end>=max((int)start,(int)listOverlap.size()-(int)tryNumber); --end){
				vector<uNumber> pathEnd;
				uint8_t errorsEnd(checkEndExhaustive(read,listOverlap[end],pathEnd,errors-errorBegin));
				if(errorsEnd+errorBegin<=errors){
					vector<uNumber> pathCover;
					bool ended(false);
					uint8_t errorCover(coverGreedy(read,listOverlap,start,end,pathCover,errors-errorsEnd-errorBegin,ended));
					if(errorCover<=errors-errorsEnd-errorBegin){
						++alignedRead;
						++readNumber;
						pathBegin.insert(pathBegin.end(), pathCover.begin(),pathCover.end());
						if (!ended){pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());}
						return pathBegin;
					}
				}
			}
		}
	}
	if(!rc){return alignReadGreedy(reverseComplement(read), overlapFound,errors, true);}
	++notAligned;
	++readNumber;
	return {};
}


vector<uNumber> Aligner::alignReadGreedyPath(const string& read, bool& overlapFound, uint8_t errors, bool rc){
	vector<pair<kmer,uint>> listOverlap(getListOverlap(read));
	if(listOverlap.empty()){
		++noOverlapRead;
		++readNumber;
		return {};
	}

	overlaps+=listOverlap.size();
	overlapFound=true;
	for(uint start(0); start<=min(tryNumber-1,(uint)listOverlap.size()-1); ++start){
		vector<uNumber> path;
		uint8_t errorBegin(checkBeginExhaustivePath(read,listOverlap[start],path,errors));
		if(errorBegin<=errors){
			for(int end((int)listOverlap.size()-1); end>=max((int)start,(int)listOverlap.size()-(int)tryNumber); --end){
				uint8_t errorsEnd(checkEndExhaustivePath(read,listOverlap[end],path,errors-errorBegin));
				if(errorsEnd+errorBegin<=errors){
					bool ended(false);
					uint8_t errorCover(coverGreedyPath(read,listOverlap,start,end,path,errors-errorsEnd-errorBegin,ended));
					if(errorCover<=errors-errorsEnd-errorBegin){
						++alignedRead;
						++readNumber;
						//~ if (!ended){pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());}
						return path;
					}
				}
			}
		}
	}
	if(!rc){return alignReadGreedyPath(reverseComplement(read), overlapFound,errors, true);}
	++notAligned;
	++readNumber;
	return {};
}


vector<uNumber> Aligner::alignReadExhaustive(const string& read, bool& overlapFound, uint8_t errors){
	vector<pair<kmer,uint>> listOverlap(getListOverlap(read));
	if(listOverlap.empty()){
		noOverlapRead++;
		readNumber++;
		return {};
	}

	overlaps+=listOverlap.size();
	overlapFound=true;
	for(uint start(0); start<listOverlap.size(); ++start){
		vector<uNumber> pathBegin;
		uint8_t errorBegin(checkBeginExhaustive(read,listOverlap[start],pathBegin,errors));
		if(errorBegin<=errors){
			vector<uNumber> pathEnd;
			uint8_t errorsEnd(checkEndExhaustive(read,listOverlap[start],pathEnd,errors-errorBegin));
			if(errorsEnd+errorBegin<=errors){
				alignedRead++;
				readNumber++;
				pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
				return pathBegin;
			}
		}
	}
	notAligned++;
	readNumber++;
	return {};
}


string compactionEnd(const string& seq1,const string& seq2, size_t k){
	size_t s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}
	string rc2(reverseComplement(seq2)),end1(seq1.substr(s1-k,k)), beg2(seq2.substr(0,k));
	if(end1==beg2){return seq1+(seq2.substr(k));}
	string begrc2(rc2.substr(0,k));
	if(end1==begrc2){return seq1+(rc2.substr(k));}
	return "";
}


string Aligner::recoverPath(vector<uNumber>& numbers,int size){
	int offset(numbers[0]);
	string path(getUnitig(numbers[1]));
	for(size_t i(2); i<numbers.size(); ++i){
		string unitig(getUnitig(numbers[i])),inter(compactionEnd(path, unitig, k-1));
		if(inter.empty()){
			cout<<"bug compaction"<<endl;
			cout<<path<<" "<<unitig<<endl;
			exit(0);
		}else{
			path=inter;
		}
	}
	path=path.substr(offset,size);
	return path;
}


uint8_t Aligner::coverGreedy(const string& read, const vector<pair<kmer,uint>>& listOverlap, const size_t start, size_t end, vector<uNumber>& path, uint8_t  errors, bool& ended){
	if(start==end){return 0;}
	int indice(0);
	uint8_t minMissmatch(errors+1);
	uNumber number2keep;

	for(size_t i(1); i+start<=end and i<=tryNumber; ++i){
		uNumber number(0);
		uint8_t missmatch(checkPair(listOverlap[start], listOverlap[i+start], read, number,errors));
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

	auto pair(mapOnRight(read,path,listOverlap[start],listOverlap,ended,start,errors));
	if(pair.second<=errors){
		if(ended){
			successMR2++;
			return pair.second;
		}
		if(pair.first>start){
			uint8_t errorrecur(coverGreedy(read, listOverlap, pair.first, end, path,errors-pair.second,ended));
			if(errorrecur+pair.second<=errors){
				successMR2++;
				return errorrecur+pair.second;
			};
		}
	}
	return errors+1;
}


uint8_t Aligner::coverGreedyPath(const string& read, const vector<pair<kmer,uint>>& listOverlap, const size_t start, size_t end, vector<uNumber>& path, uint8_t  errors, bool& ended){
	if(start==end){return 0;}
	int indice(0);
	uint8_t minMissmatch(errors+1);
	uNumber number2keep;

	for(size_t i(1); i+start<=end and i<=tryNumber; ++i){
		uNumber number(0);
		uint8_t missmatch(checkPairPaths(listOverlap[start], listOverlap[i+start], read, number,errors,path));
		if(missmatch<minMissmatch){
			indice=i;
			number2keep=number;
			minMissmatch=missmatch;
		}
	}
	if(minMissmatch<=errors){
		path.push_back(number2keep);
		return minMissmatch+coverGreedyPath(read, listOverlap, indice+start, end, path, errors-minMissmatch,ended);
	}

	auto pair(mapOnRightPath(read,path,listOverlap[start],listOverlap,ended,start,errors));
	if(pair.second<=errors){
		if(ended){
			successMR2++;
			return pair.second;
		}
		if(pair.first>start){
			uint8_t errorrecur(coverGreedyPath(read, listOverlap, pair.first, end, path,errors-pair.second,ended));
			if(errorrecur+pair.second<=errors){
				successMR2++;
				return errorrecur+pair.second;
			};
		}
	}
	return errors+1;
}


uint8_t Aligner::mapOnLeftEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors){
	string unitig, readLeft(read.substr(0,overlap.second));
	vector<uNumber> path2keep;
	if(readLeft.size()==0){return 0;}
	auto rangeUnitigs(getEnd(overlap.first));
	uint8_t miniMiss(errors+1),miniMissIndice(9);
	int offset(-2);
	bool ended(false);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()+k-1 <= unitig.size()){
			uint8_t miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				offset=unitig.size()-readLeft.size()-k+1;
				ended=true;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint8_t miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			if(miss<miniMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
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


uint8_t Aligner::mapOnLeftEndExhaustivePaths(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors){
	string unitig, readLeft(read.substr(0,overlap.second));
	vector<uNumber> path2keep;
	if(readLeft.size()==0){return 0;}
	auto rangeUnitigs(getEnd(overlap.first));
	uint8_t miniMiss(errors+1),miniMissIndice(9);
	int offset(-2);
	bool ended(false);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		auto it = find (path.begin(), path.end(), rangeUnitigs[i].second);
		if(it==path.end()){
			unitig=(rangeUnitigs[i].first);
			//case the rest of the read is too small
			if(readLeft.size()+k-1 <= unitig.size()){
				uint8_t miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
				if(miss<miniMiss){
					miniMiss=miss;
					miniMissIndice=i;
					offset=unitig.size()-readLeft.size()-k+1;
					ended=true;
				}
			}else{
				//case the read is big enough we want to recover a true overlap
				uint8_t miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
				if(miss<miniMiss){
					kmer overlapNum(str2num(unitig.substr(0,k-1)));
					vector<uNumber> possiblePath;
					miss+=mapOnLeftEndExhaustivePaths(read , possiblePath, {overlapNum,overlap.second-(unitig.size()-k+1)}, errors-miss);
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


uint8_t Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors){
	string unitig,readLeft(read.substr(0,overlap.second)),nextUnitig;
	if(readLeft.empty()){path.push_back(0);return true;}
	auto rangeUnitigs(getEnd(overlap.first));
	uint8_t miniMiss(errors+1),miniMissIndice(9);
	bool ended(false);
	int offset(0);
	kmer nextOverlap(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()+k-1 <= unitig.size()){
			uint8_t miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint8_t miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			if(miss<miniMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
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


uint8_t Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors){
	string unitig,readLeft(read.substr(overlap.second)),nextUnitig;
	if(readLeft.size()<=k-1){return true;}
	auto rangeUnitigs(getBegin(overlap.first));
	uint8_t miniMiss(errors+1), miniMissIndice(9);
	bool ended(false);
	kmer nextOverlap(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()<= unitig.size()){
			uint8_t miss(missmatchNumber(unitig.substr(0,readLeft.size()), readLeft, errors));
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
			uint8_t miss(missmatchNumber(unitig, read.substr(overlap.second,unitig.size()), errors));
			if(miss<miniMiss){
				kmer overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
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
				kmer overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
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


uint8_t Aligner::mapOnRightEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors){
	string unitig,readLeft(read.substr(overlap.second+k-1));
	vector<uNumber> path2keep;
	if(readLeft.empty()){path.push_back(0);return 0;}
	auto rangeUnitigs(getBegin(overlap.first));
	uint8_t miniMiss(errors+1), miniMissIndice(9);
	bool ended(false);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()<= unitig.size()-k+1){
			uint8_t miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint8_t miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss<miniMiss){
				kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
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
		if(ended){
			path.push_back(rangeUnitigs[miniMissIndice].second);
			path.push_back(readLeft.size()+k-1);
		}else{
			path.push_back(rangeUnitigs[miniMissIndice].second);
			path.insert(path.end(), path2keep.begin(),path2keep.end());
		}
		successMR++;
	}
	return miniMiss;
}


uint8_t Aligner::mapOnRightEndExhaustivePath(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors){
	string unitig,readLeft(read.substr(overlap.second+k-1));
	vector<uNumber> path2keep;
	if(readLeft.empty()){path.push_back(0);return 0;}
	auto rangeUnitigs(getBegin(overlap.first));
	uint8_t miniMiss(errors+1), miniMissIndice(9);
	bool ended(false);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		auto it = find (path.begin(), path.end(), rangeUnitigs[i].second);
		if(it==path.end()){
			unitig=(rangeUnitigs[i].first);
			//case the rest of the read is too small
			if(readLeft.size()<= unitig.size()-k+1){
				uint8_t miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
				if(miss<miniMiss){
					miniMiss=miss;
					miniMissIndice=i;
					ended=true;
				}
			}else{
				//case the read is big enough we want to recover a true overlap
				uint8_t miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
				if(miss<miniMiss){
					kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
					vector<uNumber> possiblePath;
					miss+=mapOnRightEndExhaustivePath(read , possiblePath, {overlapNum,overlap.second+(unitig.size()-k+1)}, errors-miss);
					if(miss<miniMiss){
						path2keep=possiblePath;
						miniMiss=miss;
						miniMissIndice=i;
						ended=false;
					}
				}
			}
		}
	}

	if(miniMiss<=errors){
		if(ended){
			path.push_back(rangeUnitigs[miniMissIndice].second);
			path.push_back(readLeft.size()+k-1);
		}else{
			path.push_back(rangeUnitigs[miniMissIndice].second);
			path.insert(path.end(), path2keep.begin(),path2keep.end());
		}
		successMR++;
	}
	return miniMiss;
}


pair<size_t,uint8_t> Aligner::mapOnRight(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,size_t start, uint8_t errors){
	string unitig, readLeft(read.substr(overlap.second+k-1)),nextUnitig;
	if(readLeft.empty()){cout<<"should not appears"<<endl;exit(0);return {start,0};}
	auto rangeUnitigs(getBegin(overlap.first));
	uint8_t miniMiss(errors+1),miniMissIndice(9);
	size_t next(start);
	kmer nextOverlapNum(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size() <= unitig.size()-k+1){
			uint8_t miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				ended=true;
				miniMiss=miss;
				miniMissIndice=i;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint8_t miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss<miniMiss){
				kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				if(miss<miniMiss){
					ended=false;
					miniMiss=miss;
					miniMissIndice=i;
					nextOverlapNum=overlapNum;
					nextUnitig=unitig;
					next=start;
					for(uint j(start+1); j<listOverlap.size(); ++j){
						if(overlapNum==listOverlap[j].first and listOverlap[j].second==overlap.second+unitig.size()-k+1){
							next=j;
						}
					}
				}
			}
		}
	}

	if(ended){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		return {start,miniMiss};
	}
	if(miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(next>start){
			return {next,miniMiss};
		}
		auto res(mapOnRight(read , path, {nextOverlapNum,overlap.second+(nextUnitig.size()-k+1)},listOverlap,ended,start, errors-miniMiss));
		return {res.first,res.second+miniMiss};
	}
	return {start,errors+1};
}



pair<size_t,uint8_t> Aligner::mapOnRightPath(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,size_t start, uint8_t errors){
	string unitig, readLeft(read.substr(overlap.second+k-1)),nextUnitig;
	if(readLeft.empty()){cout<<"should not appears"<<endl;exit(0);return {start,0};}
	auto rangeUnitigs(getBegin(overlap.first));
	uint8_t miniMiss(errors+1),miniMissIndice(9);
	size_t next(start);
	kmer nextOverlapNum(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		auto it = find (path.begin(), path.end(), rangeUnitigs[i].second);
		if(it==path.end()){
			unitig=(rangeUnitigs[i].first);
			//case the rest of the read is too small
			if(readLeft.size() <= unitig.size()-k+1){
				uint8_t miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
				if(miss<miniMiss){
					ended=true;
					miniMiss=miss;
					miniMissIndice=i;
				}
			}else{
				//case the read is big enough we want to recover a true overlap
				uint8_t miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
				if(miss<miniMiss){
					kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
					if(miss<miniMiss){
						ended=false;
						miniMiss=miss;
						miniMissIndice=i;
						nextOverlapNum=overlapNum;
						nextUnitig=unitig;
						next=start;
						for(uint j(start+1); j<listOverlap.size(); ++j){
							if(overlapNum==listOverlap[j].first and listOverlap[j].second==overlap.second+unitig.size()-k+1){
								next=j;
							}
						}
					}
				}
			}
		}
	}

	if(ended){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		return {start,miniMiss};
	}
	if(miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(next>start){
			return {next,miniMiss};
		}
		auto res(mapOnRight(read , path, {nextOverlapNum,overlap.second+(nextUnitig.size()-k+1)},listOverlap,ended,start, errors-miniMiss));
		return {res.first,res.second+miniMiss};
	}
	return {start,errors+1};
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


uint8_t Aligner::checkPair(const pair<kmer, uint>& overlap1, const pair<kmer, uint>& overlap2, const string& read, uNumber& number, uint8_t errorsAllowed){
//case where we just have to check if such a unitig exist
	if(overlap2.second-overlap1.second<k){
		int32_t positionget1, positionget2;
		auto rangeUnitigs1(getBegin(overlap1.first));
		auto rangeUnitigs2(getEnd(overlap2.first));
		for(uint i(0); i<rangeUnitigs1.size(); ++i){
			positionget1=rangeUnitigs1[i].second;
			for(uint j(0); j<rangeUnitigs2.size(); ++j){
				positionget2=rangeUnitigs2[j].second;
				if(positionget2==positionget1){
					number=positionget1;
					return 0;
				}
			}
		}
		return errorsAllowed+1;
	}

	string unitig,subRead(read.substr(overlap1.second+k-1,overlap2.second-overlap1.second-(k-1)));
	auto rangeUnitigs1(getBegin(overlap1.first));
	auto rangeUnitigs2(getEnd(overlap2.first));
	uint8_t minMissMatch(errorsAllowed+1),indice(0);
	int32_t positionget1, positionget2;
	for(uint i(0); i<rangeUnitigs1.size(); ++i){
		positionget1=rangeUnitigs1[i].second;
		for(uint j(0); j<rangeUnitigs2.size(); ++j){
			positionget2=rangeUnitigs2[j].second;
			if(positionget2==positionget1){
				unitig=getUnitig(rangeUnitigs1[i].second);
				if(unitig.size()-2*(k-1)==subRead.size()){
					uint8_t missmatch(missmatchNumber(unitig.substr(k-1,subRead.size()), subRead, errorsAllowed));
					if(missmatch<minMissMatch){
						minMissMatch=missmatch;
						indice=i;
					}
				}
			}
		}
	}

	if(minMissMatch<=errorsAllowed){
		number=(rangeUnitigs1[indice].second);
	}
	return minMissMatch;
}


uint8_t Aligner::checkPairPaths(const pair<kmer, uint>& overlap1, const pair<kmer, uint>& overlap2, const string& read, uNumber& number, uint8_t errorsAllowed, const vector<uNumber> path){
	if(overlap2.second-overlap1.second<k){
		int32_t positionget1, positionget2;
		auto rangeUnitigs1(getBegin(overlap1.first));
		auto rangeUnitigs2(getEnd(overlap2.first));
		for(uint i(0); i<rangeUnitigs1.size(); ++i){
			positionget1=rangeUnitigs1[i].second;
			for(uint j(0); j<rangeUnitigs2.size(); ++j){
				positionget2=rangeUnitigs2[j].second;
				if(positionget2==positionget1){
					auto it = find (path.begin(), path.end(), positionget1);
					if(it==path.end()){
						number=positionget1;
					}
					return 0;
				}
			}
		}
		return errorsAllowed+1;
	}
	string unitig,subRead(read.substr(overlap1.second+k-1,overlap2.second-overlap1.second-(k-1)));
	auto rangeUnitigs1(getBegin(overlap1.first));
	auto rangeUnitigs2(getEnd(overlap2.first));
	uint8_t minMissMatch(errorsAllowed+1),indice(0);
	int32_t positionget1, positionget2;
	for(uint i(0); i<rangeUnitigs1.size(); ++i){
		positionget1=rangeUnitigs1[i].second;
		for(uint j(0); j<rangeUnitigs2.size(); ++j){
			positionget2=rangeUnitigs2[j].second;
			if(positionget2==positionget1){
				auto it = find (path.begin(), path.end(), positionget1);
				if(it==path.end()){
					unitig=getUnitig(rangeUnitigs1[i].second);
					if(unitig.size()-2*(k-1)==subRead.size()){
						uint8_t missmatch(missmatchNumber(unitig.substr(k-1,subRead.size()), subRead, errorsAllowed));
						if(missmatch<minMissMatch){
							minMissMatch=missmatch;
							indice=i;
						}
					}
				}
			}
		}
	}

	if(minMissMatch<=errorsAllowed){
		number=(rangeUnitigs1[indice].second);
	}
	return minMissMatch;
}


uint8_t Aligner::checkBeginExhaustive(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors){
	if(overlap.second==0){
		path.push_back(0);
		return 0;
	}
	string readLeft(read.substr(0,overlap.second)),unitig;
	auto rangeUnitigs(getEnd(overlap.first));
	uint8_t minMiss(errors+1),indiceMinMiss(0);
	int offset(-2);
	bool ended(false);
	vector<uNumber> path2keep;

	if(partial & rangeUnitigs.empty()){
		//if(!path.empty()){
			return 0;
		//}
	}

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		if(unitig.size()-k+1>=readLeft.size()){
			uint8_t miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=(unitig.size()-readLeft.size()-k+1);
			}
		}else{
			uint8_t miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			if(miss<minMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
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


uint8_t Aligner::checkBeginExhaustivePath(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors){
	if(overlap.second==0){
		path.push_back(0);
		return 0;
	}
	string readLeft(read.substr(0,overlap.second)),unitig;
	auto rangeUnitigs(getEnd(overlap.first));
	uint8_t minMiss(errors+1),indiceMinMiss(0);
	int offset(-2);
	bool ended(false);
	vector<uNumber> path2keep;

	if(partial & rangeUnitigs.empty()){
		//if(!path.empty()){
			return 0;
		//}
	}

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		auto it = find (path.begin(), path.end(), rangeUnitigs[i].second);
		if(it==path.end()){
			unitig=(rangeUnitigs[i].first);
			if(unitig.size()-k+1>=readLeft.size()){
				uint8_t miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
				if(miss<minMiss){
					minMiss=miss;
					indiceMinMiss=i;
					ended=true;
					offset=(unitig.size()-readLeft.size()-k+1);
				}
			}else{
				uint8_t miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
				if(miss<minMiss){
					kmer overlapNum(str2num(unitig.substr(0,k-1)));
					vector<uNumber> possiblePath;
					miss+=mapOnLeftEndExhaustivePaths(read, possiblePath, {overlapNum,overlap.second-(unitig.size()-k+1)},errors-miss);
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
	}

	if(minMiss<=errors){
		if(ended){
			auto it = find (path.begin(), path.end(), rangeUnitigs[indiceMinMiss].second);
			if(it==path.end()){
				path.push_back(offset);
				path.push_back(rangeUnitigs[indiceMinMiss].second);
			}
		}else{
			path.insert(path.end(), path2keep.begin(), path2keep.end());
			path.push_back(rangeUnitigs[indiceMinMiss].second);
		}
	}
	return minMiss;
}


uint8_t Aligner::checkBeginGreedy(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors){
	if(overlap.second==0){path.push_back(0);return 0;}
	string readLeft(read.substr(0,overlap.second)),unitig,nextUnitig;
	auto rangeUnitigs(getEnd(overlap.first));
	uint8_t minMiss(errors+1),indiceMinMiss(0);
	bool ended(false);
	int offset(0);
	kmer nextOverlap(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		if(unitig.size()-k+1>=readLeft.size()){
			uint8_t miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			uint8_t miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			if(miss<minMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
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


uint8_t Aligner::checkEndExhaustive(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors){
	string readLeft(read.substr(overlap.second+k-1)),unitig;
	vector<uNumber> path2keep;
	if(readLeft.empty()){
		path.push_back(0);
		return 0;
	}
	auto rangeUnitigs(getBegin(overlap.first));
	uint8_t minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
	int offset(-2);
	if(partial & rangeUnitigs.empty()){
		//if(!path.empty()){
			return 0;
		//}
	}

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		if(unitig.size()-k+1>=readLeft.size()){
			uint8_t miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=readLeft.size()+k-1;
			}
		}else{
			uint8_t miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss<minMiss){
				kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
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
		if(ended){
			path.push_back(rangeUnitigs[indiceMinMiss].second);
			path.push_back(offset);
		}else{
			path.push_back(rangeUnitigs[indiceMinMiss].second);
			path.insert(path.end(), path2keep.begin(),path2keep.end());
		}
	}
	return minMiss;
}


uint8_t Aligner::checkEndExhaustivePath(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors){
	string readLeft(read.substr(overlap.second+k-1)),unitig;
	vector<uNumber> path2keep;
	if(readLeft.empty()){
		path.push_back(0);
		return 0;
	}
	auto rangeUnitigs(getBegin(overlap.first));
	uint8_t minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
	int offset(-2);
	if(partial & rangeUnitigs.empty()){
		//if(!path.empty()){
			return 0;
		//}
	}

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		auto it = find (path.begin(), path.end(), rangeUnitigs[i].second);
		if(it==path.end()){
			unitig=(rangeUnitigs[i].first);
			if(unitig.size()-k+1>=readLeft.size()){
				uint8_t miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
				if(miss<minMiss){
					minMiss=miss;
					indiceMinMiss=i;
					ended=true;
					offset=readLeft.size()+k-1;
				}
			}else{
				uint8_t miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
				if(miss<minMiss){
					kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
					vector<uNumber> possiblePath;
					miss+=mapOnRightEndExhaustivePath(read, possiblePath, {overlapNum,overlap.second+(unitig.size()-k+1)},errors-miss);
					if(miss<minMiss){
						path2keep=possiblePath;
						minMiss=miss;
						indiceMinMiss=i;
						ended=false;
					}
				}
			}
		}
	}

	if(minMiss<=errors){
		if(ended){
			auto it = find (path.begin(), path.end(), rangeUnitigs[indiceMinMiss].second);
			if(it==path.end()){
				path.push_back(rangeUnitigs[indiceMinMiss].second);
				path.push_back(offset);
			}
		}else{
			path.push_back(rangeUnitigs[indiceMinMiss].second);
			path.insert(path.end(), path2keep.begin(),path2keep.end());
		}
	}
	return minMiss;
}



uint8_t Aligner::checkEndGreedy(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors){
	string readLeft(read.substr(overlap.second+k-1)),unitig,nextUnitig;
	if(readLeft.size()<=k-1){return 0;}
	auto rangeUnitigs(getBegin(overlap.first));
	uint8_t minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
	kmer nextOverlap(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		if(unitig.size()-k+1>=readLeft.size()){
			uint8_t miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
			}
		}else{
			uint8_t miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss<minMiss){
				kmer overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
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


void Aligner::update(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdate;
}


void Aligner::updateRC(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-4));
}


vector<pair<kmer,uint>> Aligner::getListOverlap(const string& read){
	vector<pair<kmer,uint>> listOverlap;
	string overlap(read.substr(0,k-1));
	kmer num(str2num(overlap)),rcnum(str2num(reverseComplement(overlap))), rep(min(num, rcnum));

	for(uint i(0);;++i){
		if(left.count(rep)!=0){
			listOverlap.push_back({num,i});
		}else{
			if(right.count(rep)!=0){
				listOverlap.push_back({num,i});
			}
		}
		if(i+k-1<read.size()){
			update(num,read[i+k-1]);
			updateRC(rcnum,read[i+k-1]);
			rep=(min(num, rcnum));
		}else{
			return listOverlap;
		}
	}
	return listOverlap;
}


void Aligner::alignPartExhaustive(){
	vector<pair<string,string>> multiread;
	vector<uNumber> path;
	string read,header;
	bool overlapFound(false);
	while(!readFile.eof()){
		readMutex.lock();
		{
			getReads(multiread,10000);
		}
		readMutex.unlock();
		for(size_t i(0);i<multiread.size();++i){
			header=multiread[i].first;
			read=multiread[i].second;
			overlapFound=false;
			path=alignReadExhaustive(read,overlapFound,errorsMax);
			if(path.size()!=0){
				pathMutex.lock();
				{
					printPath(path,&pathFile);
				}
				pathMutex.unlock();
			}else{
				if(!overlapFound){
						noOverlapMutex.lock();
						{
							noOverlapFile<<header<<endl<<read<<endl;
						}
						noOverlapMutex.unlock();
				}else{
					notMappedMutex.lock();
					{
						notMappedFile<<header<<endl<<read<<endl;
					}
					notMappedMutex.unlock();
				}
			}
		}
		if(iter++%100==0){
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
}


void Aligner::alignPartGreedy(){
	vector<pair<string,string>> multiread;
	vector<uNumber> path;
	string read,header;
	bool overlapFound(false);
	while(!readFile.eof()){
		readMutex.lock();
		{
			getReads(multiread,10000);
		}
		readMutex.unlock();
		for(size_t i(0);i<multiread.size();++i){
			overlapFound=false;
			header=multiread[i].first;
			read=multiread[i].second;
			if(pathOption){
				path=alignReadGreedyPath(read,overlapFound,errorsMax,false);
			}else{
				path=alignReadGreedy(read,overlapFound,errorsMax,false);
			}
			if(path.size()!=0){
				pathMutex.lock();
				{
					printPath(path,&pathFile);
				}
				pathMutex.unlock();
			}else{
				if(!overlapFound){
					noOverlapMutex.lock();
					{
						notMappedFile<<header<<endl<<read<<endl;
					}
					noOverlapMutex.unlock();
				}else{
					notMappedMutex.lock();
					{
						notMappedFile<<header<<endl<<read<<endl;
					}
					notMappedMutex.unlock();
				}
			}
		}
		if(iter++%100==0){
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
}


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
			string beg(line.substr(0,k-1)),rcBeg(reverseComplement(beg));
			if(beg<=rcBeg){
				left[str2num(beg)].push_back(position);
			}else{
				right[str2num(rcBeg)].push_back(position);
			}
			string end(line.substr(size-k+1,k-1)), rcEnd(reverseComplement(end));
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
}


void Aligner::alignAll(bool greedy, const string& reads){
	startChrono=chrono::system_clock::now();
	uint last(0);
	string file;
	for(uint i(0);i<reads.size();++i){
		if(reads[i]==','){
			file=reads.substr(last,i-last);
			readFile.close();
			readFile.open(file);
			cout<<file<<endl;
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
			last=i+1;
		}
	}
	file=reads.substr(last);
	readFile.close();
	readFile.open(file);
	cout<<file<<endl;
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
	cout<<"Overlap per reads : "<<(overlaps)/(alignedRead+notAligned)<<endl;
	cout<<endl;
}
