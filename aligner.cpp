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
#include "utils.h"
#include "alignerPaths.cpp"
#include "alignerExhaustive.cpp"
#include "alignerGreedy.cpp"
#include <thread>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>


using namespace std;


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
					if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G' and read[j]!='N'){
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
						if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G' and read[j]!='N'){
							fail=true;
							break;
						}
					}
					if(!fail){
						if(read.size()>k){
							reads.push_back({header,read});
						}
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
							if(read[j]!='A' and read[j]!='C' and read[j]!='T' and read[j]!='G' and read[j]!='N'){
								fail=true;
								break;
							}
						}
						if(!fail){
							if(read.size()>k){
								reads.push_back({header,read});
							}
						}
					}
					return;
				}
			}
		}
	}
}


kmer Aligner::getRepresentNum(const string& str){
	// string rc(reverseComplement(str));
	kmer a(str2num(str));
	kmer b(rcb(a,k-1));
	return ((b<=a) ? b : a);
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


vector<pair<string,uNumber>> Aligner::getEnd(kmer bin){
	vector<pair<string,uNumber>> result;
	kmer rc(rcb(bin,k-1));
	string unitig;
	uNumber num;
	if(bin<=rc){
		if(right.count(bin)!=0){
			for(uint i(0); i<right[bin].size(); ++i){
				num=right[bin][i];
				unitig=(getUnitig(num));
				if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
					result.push_back({unitig,num});
				}else{
					result.push_back({reverseComplements(unitig),-num});
				}
			}
		}else{
			return {};
		}
	}else{
		if(left.count(rc)!=0){
			for(uint i(0); i<left[rc].size(); ++i){
				num=left[rc][i];
				unitig=(getUnitig(num));
				if(str2num(unitig.substr(unitig.size()-k+1,k-1))==bin){
					result.push_back({unitig,num});
				}else{
					result.push_back({reverseComplements(unitig),-num});
				}
			}
		}else{
			return {};
		}
	}
	return result;
}


//TODO can count be removed with the use of at ?
vector<pair<string,uNumber>> Aligner::getBegin(kmer bin){
	vector<pair<string,uNumber>> result;
	kmer rc(rcb(bin,k-1));
	string unitig;
	if(bin<=rc){
		if(left.count(bin)!=0){
			for(uint i(0); i<left[bin].size(); ++i){
				unitig=(getUnitig(left.at(bin)[i]));
				if(str2num(unitig.substr(0,k-1))==bin){
					result.push_back({unitig,left.at(bin)[i]});
				}else{
					result.push_back({reverseComplements(unitig),-left.at(bin)[i]});
				}
			}
		}else{
			return {};
		}
	}else{
		if(right.count(rc)!=0){
			for(uint i(0); i<right[rc].size(); ++i){
				unitig=(getUnitig(right.at(rc)[i]));
				if(str2num(unitig.substr(0,k-1))==bin){
					result.push_back({unitig,right.at(rc)[i]});
				}else{
					result.push_back({reverseComplements(unitig),-right.at(rc)[i]});
				}
			}
		}else{
			return {};
		}
	}
	return result;
}

//
// vector<pair<string,uNumber>> Aligner::getBeginOpti(kmer bin, uNumber last){
// 	// cout<<"getbegin"<<endl;
// 	vector<pair<string,uNumber>> result;
// 	// kmer rc(rcb(bin,k-1));
// 	uNumber* begin(unitigs[(last>=0) ? last : -last].begin);
// 	uNumber* end(unitigs[(last>=0) ? last : -last].end);
// 	// cout<<lastUnitig.str<<endl;
// 	// kmer beg(str2num(lastUnitig.str.substr(0,k-1)));
// 	// kmer end(str2num(lastUnitig.str.substr(lastUnitig.str.size()-k+1)));
// 	// cout<<lastUnitig.str.substr(lastUnitig.str.size()-k+1)<<endl;
// 	// if(bin<=rc){
// 	if(last>0){
// 		// cout<<1;
// 		for(uint i(0);i<4;++i){
// 			// cout<<lastUnitig.begin[i]<<endl;
// 			if(begin[i]!=0){
// 				if(begin[i]>0){
// 					result.push_back({unitigs[begin[i]].str,begin[i]});
// 				// cout<<"add"<<unitigs[lastUnitig.begin[i]].str<<endl;
// 				}else{
// 					result.push_back({reverseComplements(unitigs[-begin[i]].str),begin[i]});
// 				}
// 			}
// 		}
// 	}else{
// 		// cout<<2<<endl;;
// 		for(uint i(0);i<4;++i){
// 			if(end[i]!=0){
// 				if(end[i]>0){
// 					// cout<<
// 					result.push_back({reverseComplements(unitigs[end[i]].str),-end[i]});
// 				}else{
// 					result.push_back({(unitigs[-end[i]].str),-end[i]});
// 				}
// 			}
// 		}
// 	}
// 	return result;
// }
//
//
// vector<pair<string,uNumber>> Aligner::getEndOpti(kmer bin, uNumber last){
// 		// cout<<"getend"<<endl;
// 	vector<pair<string,uNumber>> result;
// 	uNumber* begin(unitigs[(last>=0) ? last : -last].begin);
// 	uNumber* end(unitigs[(last>=0) ? last : -last].end);
// 		// auto lastUnitig(unitigs[last]);
// 	// kmer rc(rcb(bin,k-1));
// 	// cout<<lastUnitig.str<<endl;
// 	// kmer end(str2num(lastUnitig.str.substr(lastUnitig.str.size()-k+1)));
// 	// kmer beg(str2num(lastUnitig.str.substr(0,k-1)));
// 	// cout<<lastUnitig.str.substr(lastUnitig.str.size()-k)<<endl;
// 	// cout<<bin<<endl;
// 	// if(bin<=rc){
// 	if(last>0){
// 		for(uint i(0);i<4;++i){
// 			if(end[i]!=0){
// 				if(end[i]>0){
// 					result.push_back({unitigs[end[i]].str,end[i]});
// 				}else{
// 					result.push_back({reverseComplements(unitigs[-end[i]].str),end[i]});
// 				}
// 			}
// 		}
// 	}else{
// 		// cout<<"lol"<<endl;
// 		for(uint i(0);i<4;++i){
// 			if(begin[i]!=0){
// 				if(begin[i]>0){
// 					result.push_back({reverseComplements(unitigs[begin[i]].str),-begin[i]});
// 				// cout<<"add"<<unitigs[lastUnitig.begin[i]].str<<endl;
// 			}else{
// 				result.push_back({(unitigs[-begin[i]].str),-begin[i]});
// 				}
// 			}
// 		}
// 	}
// 		// cout<<4<<endl;
// 	return result;
// }
//

string Aligner::recoverPath(vector<uNumber>& numbers,uint size){
	int offset(numbers[0]);
	string path(getUnitig(numbers[1]));
	for(uint i(2); i<numbers.size(); ++i){
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


string Aligner::getUnitig(int position){
	// if(fullMemory){
	// 	if(position>0){
	// 		return unitigs[position].str;
	// 	}
	// 	return reverseComplements(unitigs[-position].str);
	// }
	if(fullMemory){
		if(position>0){
			return unitigs[position];
		}
		return reverseComplements(unitigs[-position]);
	}
	return "";

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
	kmer num(str2num(overlap)),rcnum(rcb(num,k-1)), rep(min(num, rcnum));

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


vector<pair<kmer,uint>> Aligner::getNOverlap(const string& read, uint n){
	vector<pair<kmer,uint>> listOverlap;
	kmer num(str2num(read.substr(0,k-1))),rcnum(rcb(num,k-1)), rep(min(num, rcnum));
	for(uint i(0);;++i){
		if(left.unordered_map::count(rep)!=0){
			listOverlap.push_back({num,i});
		}else{
			if(right.unordered_map::count(rep)!=0){
				listOverlap.push_back({num,i});
			}
		}
		if(listOverlap.size()>=n){
			return listOverlap;
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


vector<overlapStruct> Aligner::getListOverlapCache(const string& read){
	vector<overlapStruct> listOverlap;
	string overlap(read.substr(0,k-1));
	kmer num(str2num(overlap)),rcnum(rcb(num,k-1)), rep(min(num, rcnum));
	overlapStruct over;
	bool push;
	for(uint i(0);;++i){
		push=false;
		over.unitig.clear();
		over.unitigNumbers.clear();
		if(right.count(rep)!=0){
			for(uint ii(0); ii<right[rep].size(); ++ii){
				over.seq=num;
				over.pos=i;
				over.unitig.push_back(getUnitig(right[rep][ii]));
				over.unitigNumbers.push_back(right[rep][ii]);
			}
			push=true;
		}
		if(left.count(rep)!=0){
			for(uint ii(0); ii<left[rep].size(); ++ii){
				over.seq=num;
				over.pos=i;
				over.unitig.push_back(getUnitig(left[rep][ii]));
				over.unitigNumbers.push_back(left[rep][ii]);
			}
			push=true;
		}

		if(push){listOverlap.push_back(over);}
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


//TODO can be multithreaded with N hash table and a xorshift
void Aligner::indexUnitigsAux(){
	uint position(0),count(0);
	string line, overlap1, overlap2, waste;
	// unitigs.push_back({"",{0,0,0,0},{0,0,0,0}});
	unitigs.push_back("");
	while(!unitigFile.eof()){
		if(!fullMemory){
			position=unitigFile.tellg();
		}
		getline(unitigFile,line);
		getline(unitigFile,line);
		// unitigFile.seekg(1, ios::cur);
		// transform(line.begin(), line.end(), line.begin(), ::toupper);
		if(line.size()<k){
			break;
		}else{
			if(++count%1000000==0){cout<<count/1000000<<"M unitigs treated"<<endl;}
			if(fullMemory){
				// unitigs.push_back({line,{0,0,0,0},{0,0,0,0}});
				// cout<<line<<endl;
				unitigs.push_back(line);
				position=unitigs.size()-1;
			}
			++unitigNumber;

			kmer beg(str2num(line.substr(0,k-1))),rcBeg(rcb(beg,k-1));
			if(beg<=rcBeg){
				left[(beg)].push_back(position);
			}else{
				right[(rcBeg)].push_back(position);
			}
			kmer end(str2num(line.substr(line.size()-k+1,k-1))),rcEnd(rcb(end,k-1));
			if(end<=rcEnd){
				right[(end)].push_back(position);
			}else{
				left[(rcEnd)].push_back(position);
			}
		}
	}
}


void Aligner::indexUnitigs(){
	uint nbThreads(1);
	vector<thread> threads;
	for (uint i(0); i<nbThreads; ++i){
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

	cout<<"Ended"<<endl;
	cout<<"Read : "<<readNumber<<endl;
	cout<<"No Overlap : "<<noOverlapRead<<" Percent : "<<(100*float(noOverlapRead))/readNumber<<endl;
	cout<<"Got Overlap : "<<alignedRead+notAligned<<" Percent : "<<(100*float(alignedRead+notAligned))/readNumber<<endl;
	cout<<"Overlap and Aligned : "<<alignedRead<<" Percent : "<<(100*float(alignedRead))/(alignedRead+notAligned)<<endl;
	cout<<"Overlap but no aligne: "<<notAligned<<" Percent : "<<(100*float(notAligned))/(alignedRead+notAligned)<<endl;
	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
	cout<<"Reads/seconds : "<<readNumber/(chrono::duration_cast<chrono::seconds>(waitedFor).count()+1)<<endl;
	cout<<endl;
}


// void Aligner::knowNeighbour(){
// 	string unitig;
// 	vector<pair<string,uNumber>> rangeUnitigs;
// 	for(uint i(1); i<unitigs.size();++i){
// 		unitig=unitigs[i].str;
// 		// cout<<unitig<<endl;
// 		// cout<<unitig.substr(0,k-1)<<" "<<unitig.substr(unitig.size()-k+1)<<endl;
// 		kmer beg(str2num((unitig.substr(0,k-1))));
// 		kmer end(str2num((unitig.substr(unitig.size()-k+1))));
// 		rangeUnitigs=(getBegin(end));
// 		// cout<<"beg"<<endl;
// 		for(uint ii(0); ii<rangeUnitigs.size(); ++ii){
// 			unitigs[i].begin[ii]=rangeUnitigs[ii].second;
// 			// cout<<rangeUnitigs1[ii].second<<endl;;
// 		}
// 		// cout<<"end"<<endl;
// 		rangeUnitigs=(getEnd(beg));
// 		for(uint ii(0); ii<rangeUnitigs.size(); ++ii){
// 			unitigs[i].end[ii]=rangeUnitigs[ii].second;
// 			// cout<<rangeUnitigs2[ii].second<<endl;;
// 		}
// 		// cout<<"end"<<endl;
// 	}
// 	cout<<endl<<endl;
// }
