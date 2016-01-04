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

//PATH ALGORIGHMS
vector<uNumber> Aligner::alignReadGreedyPath(const string& read, bool& overlapFound, uint8_t errors, bool rc){
	vector<pair<kmer,uint>> listOverlap(getListOverlap(read));
	if(listOverlap.empty()){++noOverlapRead;++readNumber;return {};}
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
	// if(!rc){return alignReadGreedyPath(reverseComplement(read), overlapFound,errors, true);}
	++notAligned;
	++readNumber;
	return {};
}


vector<uNumber> Aligner::alignReadExhaustivePath(const string& read, bool& overlapFound, uint8_t errors){
	vector<pair<kmer,uint>> listOverlap(getListOverlap(read));
	if(listOverlap.empty()){++noOverlapRead;++readNumber;return {};}
	overlaps+=listOverlap.size();
	overlapFound=true;
	for(uint start(0); start<listOverlap.size(); ++start){
		vector<uNumber> pathBegin;
		uint8_t errorBegin(checkBeginExhaustivePath(read,listOverlap[start],pathBegin,errors));
		if(errorBegin<=errors){
			vector<uNumber> pathEnd;
			uint8_t errorsEnd(checkEndExhaustivePath(read,listOverlap[start],pathEnd,errors-errorBegin));
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


uint8_t Aligner::mapOnLeftEndExhaustivePaths(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors){
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
