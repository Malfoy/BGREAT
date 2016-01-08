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





vector<uNumber> Aligner::alignReadExhaustive(const string& read, bool& overlapFound, uint errors){
	vector<pair<kmer,uint>> listOverlap(getListOverlap(read));
	if(listOverlap.empty()){++noOverlapRead;++readNumber;return {};}
	overlaps+=listOverlap.size();
	overlapFound=true;
	vector<uNumber> pathBegin,pathEnd;
	for(uint start(0); start<(uint)listOverlap.size(); ++start){
		pathBegin={};
		uint errorBegin(checkBeginExhaustive(read,listOverlap[start],pathBegin,errors));
		if(errorBegin<=errors){
			pathEnd={};
			uint errorsEnd(checkEndExhaustive(read,listOverlap[start],pathEnd,errors-errorBegin));
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


uint Aligner::mapOnRightEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	string unitig,readLeft(read.substr(overlap.second+k-1));
	vector<uNumber> path2keep;
	if(readLeft.empty()){path.push_back(0);return 0;}
	auto rangeUnitigs(getBegin(overlap.first));
	uint miniMiss(errors+1), miniMissIndice(9);
	bool ended(false);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()<= unitig.size()-k+1){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
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


uint Aligner::mapOnLeftEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	string unitig, readLeft(read.substr(0,overlap.second));
	vector<uNumber> path2keep;
	if(readLeft.size()==0){return 0;}
	auto rangeUnitigs(getEnd(overlap.first));
	uint miniMiss(errors+1),miniMissIndice(9);
	int offset(-2);
	bool ended(false);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size()+k-1 <= unitig.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				offset=unitig.size()-readLeft.size()-k+1;
				ended=true;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
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


uint Aligner::checkBeginExhaustive(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	if(overlap.second==0){path.push_back(0);return 0;}
	string readLeft(read.substr(0,overlap.second)),unitig;
	auto rangeUnitigs(getEnd(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(0);
	int offset(-2);
	bool ended(false);
	vector<uNumber> path2keep;

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=(unitig.size()-readLeft.size()-k+1);
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
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


uint Aligner::checkEndExhaustive(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	string readLeft(read.substr(overlap.second+k-1)),unitig;
	vector<uNumber> path2keep;
	if(readLeft.empty()){
		path.push_back(0);
		return 0;
	}
	auto rangeUnitigs(getBegin(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(9);
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
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=readLeft.size()+k-1;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
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
		for(uint i(0);i<multiread.size();++i){
			header=multiread[i].first;
			read=multiread[i].second;
			overlapFound=false;
			if(pathOption){
				path=alignReadExhaustivePath(read,overlapFound,errorsMax);
			}else{
				path=alignReadExhaustive(read,overlapFound,errorsMax);
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
						noOverlapFile<<header<<endl<<read<<endl;
						//~ notMappedFile<<read<<endl;
					}
					noOverlapMutex.unlock();
				}else{
					notMappedMutex.lock();
					{
						notMappedFile<<header<<endl<<read<<endl;
						//~ notMappedFile<<read<<endl;
					}
					notMappedMutex.unlock();
				}
			}
		}
		if(iter++%10==0){
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
