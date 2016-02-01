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


vector<uNumber> Aligner::alignReadGreedy(const string& read, bool& overlapFound, uint errors, bool r){
	vector<pair<kmer,uint>> listOverlap(getNOverlap(read,tryNumber));
	if(listOverlap.empty()){++noOverlapRead;++readNumber;return {};}
	overlapFound=true;
	vector<uNumber> pathBegin,pathEnd;
	for(uint start(0); start<(uint)listOverlap.size(); ++start){
		pathBegin={};
		uint errorBegin(checkBeginGreedy(read,listOverlap[start],pathBegin,errors));
		if(errorBegin<=errors){
			pathEnd={};
			uint errorsEnd(checkEndGreedy(read,listOverlap[start],pathEnd,errors-errorBegin));
			if(errorsEnd+errorBegin<=errors){
				alignedRead++;
				readNumber++;
				reverse(pathBegin.begin(),pathBegin.end());
				pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
				return pathBegin;
			}
		}
	}
	//TODO TEST this
	//if(!rc){return alignReadGreedy(reverseComplements(read), overlapFound,errors, true);}
	notAligned++;
	readNumber++;
	return {};
}


// vector<uNumber> Aligner::alignReadGreedyCover(const string& read, bool& overlapFound, uint errors, bool rc){
// 	vector<pair<kmer,uint>> listOverlap(getListOverlap(read));
// 	if(listOverlap.empty()){++noOverlapRead;++readNumber;return {};}
// 	overlaps+=listOverlap.size();
// 	overlapFound=true;
// 	for(uint start(0); start<min(tryNumber,(uint)listOverlap.size()); ++start){
// 		vector<uNumber> pathBegin;
// 		uint errorBegin(checkBeginGreedy(read,listOverlap[start],pathBegin,errors));
// 		if(errorBegin<=errors){
// 			for(int end((int)listOverlap.size()-1); end>=max((int)start,(int)listOverlap.size()-(int)tryNumber); --end){
// 				vector<uNumber> pathEnd;
// 				uint errorsEnd(checkEndGreedy(read,listOverlap[end],pathEnd,errors-errorBegin));
// 				if(errorsEnd+errorBegin<=errors){
// 					vector<uNumber> pathCover;
// 					bool ended(false);
// 					uint errorCover(coverGreedy(read,listOverlap,start,end,pathCover,errors-errorsEnd-errorBegin,ended));
// 					if(errorCover<=errors-errorsEnd-errorBegin){
// 						++alignedRead;
// 						++readNumber;
// 						pathBegin.insert(pathBegin.end(), pathCover.begin(),pathCover.end());
// 						if (!ended){pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());}
// 						return pathBegin;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	//TODO TEST this
// 	//if(!rc){return alignReadGreedy(reverseComplements(read), overlapFound,errors, true);}
// 	++notAligned;
// 	++readNumber;
// 	return {};
// }

//
// vector<uNumber> Aligner::alignReadGreedyCache(const string& read, bool& overlapFound, uint errors, bool rc){
// 	auto listOverlap(getListOverlapCache(read));
// 	if(listOverlap.empty()){++noOverlapRead;++readNumber;return {};}
// 	overlaps+=listOverlap.size();
// 	overlapFound=true;
// 	for(uint start(0); start<=min(tryNumber-1,(uint)listOverlap.size()-1); ++start){
// 		vector<uNumber> pathBegin;
// 		uint errorBegin(checkBeginGreedyCache(read,listOverlap[start],pathBegin,errors));
// 		if(errorBegin<=errors){
// 			for(int end((int)listOverlap.size()-1); end>=max((int)start,(int)listOverlap.size()-(int)tryNumber); --end){
// 				vector<uNumber> pathEnd;
// 				uint errorsEnd(checkEndGreedyCache(read,listOverlap[end],pathEnd,errors-errorBegin));
// 				if(errorsEnd+errorBegin<=errors){
// 					vector<uNumber> pathCover;
// 					bool ended(false);
// 					uint errorCover(coverGreedyCache(read,listOverlap,start,end,pathCover,errors-errorsEnd-errorBegin,ended));
// 					if(errorCover<=errors-errorsEnd-errorBegin){
// 						++alignedRead;
// 						++readNumber;
// 						pathBegin.insert(pathBegin.end(), pathCover.begin(),pathCover.end());
// 						if (!ended){pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());}
// 						return pathBegin;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	//TODO TEST this
// 	//if(!rc){return alignReadGreedy(reverseComplements(read), overlapFound,errors, true);}
// 	++notAligned;
// 	++readNumber;
// 	return {};
// }
//
//
// uint Aligner::coverGreedy(const string& read, const vector<pair<kmer,uint>>& listOverlap, const uint start, uint end, vector<uNumber>& path, uint  errors, bool& ended){
// 	if(start==end){return 0;}
// 	int indice(0);
// 	uint minMissmatch(errors+1);
// 	uNumber number2keep;
//
// 	for(uint i(1); i+start<=end and i<=tryNumber; ++i){
// 		uNumber number(0);
// 		uint missmatch(checkPair(listOverlap[start], listOverlap[i+start], read, number,errors));
// 		if(missmatch<minMissmatch){
// 			indice=i;
// 			number2keep=number;
// 			minMissmatch=missmatch;
// 		}
// 	}
// 	if(minMissmatch<=errors){
// 		path.push_back(number2keep);
// 		return minMissmatch+coverGreedy(read, listOverlap, indice+start, end, path, errors-minMissmatch,ended);
// 	}
// 	auto pair(mapOnRight(read,path,listOverlap[start],listOverlap,ended,start,errors));
// 	if(pair.second<=errors){
// 		if(ended){
// 			successMR2++;
// 			return pair.second;
// 		}
// 		if(pair.first>start){
// 			uint errorrecur(coverGreedy(read, listOverlap, pair.first, end, path,errors-pair.second,ended));
// 			if(errorrecur+pair.second<=errors){
// 				successMR2++;
// 				return errorrecur+pair.second;
// 			};
// 		}
// 	}
// 	return errors+1;
// }
//
//
// uint Aligner::coverGreedyCache(const string& read, const vector<overlapStruct>& listOverlap, const uint start, uint end, vector<uNumber>& path, uint  errors, bool& ended){
// 	if(start==end){return 0;}
// 	int indice(0);
// 	uint minMissmatch(errors+1);
// 	uNumber number2keep;
// 	for(uint i(1); i+start<=end and i<=tryNumber; ++i){
// 		uNumber number(0);
// 		uint missmatch(checkPairCache(listOverlap[start], listOverlap[i+start], read, number,errors));
// 		if(missmatch<minMissmatch){
// 			indice=i;
// 			number2keep=number;
// 			minMissmatch=missmatch;
// 		}
// 	}
// 	if(minMissmatch<=errors){
// 		path.push_back(number2keep);
// 		return minMissmatch+coverGreedyCache(read, listOverlap, indice+start, end, path, errors-minMissmatch,ended);
// 	}
// 	auto pair(mapOnRightCache(read,path,listOverlap[start],listOverlap,ended,start,errors));
// 	if(pair.second<=errors){
// 		if(ended){
// 			successMR2++;
// 			return pair.second;
// 		}
// 		if(pair.first>start){
// 			uint errorrecur(coverGreedyCache(read, listOverlap, pair.first, end, path,errors-pair.second,ended));
// 			if(errorrecur+pair.second<=errors){
// 				successMR2++;
// 				return errorrecur+pair.second;
// 			};
// 		}
// 	}
// 	return errors+1;
// }
//

uint Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	// if(overlap.second==0){path.push_back(0);return 0;}
	string unitig,readLeft(read.substr(0,overlap.second)),nextUnitig;
	// cout<<"moleg"<<endl;
	// cout<<read.substr(0,overlap.second+k)<<endl;
	auto rangeUnitigs(getEndOpti(overlap.first,path.back()));
	// auto rangeUnitigs(getEnd(overlap.first));
	// cout<<unitigs[path.back()].str<<endl;
	// auto rangeUnitigs(getEnd(overlap.first));
	uint miniMiss(errors+1),miniMissIndice(9);
	bool ended(false);
	int offset(0);
	kmer nextOverlap(0);
	// cout<<"go"<<endl;
	// for(uint i(0); i<rangeUnitigs2.size(); ++i){
	// 	cout<<(rangeUnitigs2[i].first)<<endl;
	// }
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		// if(rangeUnitigs[i].first!=rangeUnitigs2[i].first){
		// 	cout<<"lol1"<<endl;
		// 	// cin.get();
		// }
		// cout<<(rangeUnitigs[i].first)<<endl;

		//case the rest of the read is too small
		if(readLeft.size()+k-1 <= unitig.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
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
	// cout<<"end"<<endl;

	if (miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){
			path.push_back(offset);
			return miniMiss;
		}
		miniMiss+=mapOnLeftEndGreedy(read , path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, errors-miniMiss);
	}
	return miniMiss;
}


uint Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	// cout<<"moreg"<<endl;
	string unitig,readLeft(read.substr(overlap.second)),nextUnitig;
	auto rangeUnitigs(getBeginOpti(overlap.first,path.back()));
	// auto rangeUnitigs(getBegin(overlap.first));
	uint miniMiss(errors+1), miniMissIndice(9);
	bool ended(false);
	// int offset(0);
	kmer nextOverlap(0);
	// cout<<"go"<<endl;
	// for(uint i(0); i<rangeUnitigs2.size(); ++i){
	// 	cout<<(rangeUnitigs2[i].first)<<endl;
	// }
	// cout<<"true"<<endl;
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		// bool stop(false);
		// if(rangeUnitigs[i].first!=rangeUnitigs2[i].first){
		// 	cout<<read<<endl;
		// 	cout<<"lol2"<<endl;
		// 	cout<<rangeUnitigs[i].first<<endl<<rangeUnitigs2[i].first<<endl;
		// 	cout<<rangeUnitigs[i].second<<" "<<rangeUnitigs2[i].second<<endl;
		// 	stop=true;
		// }
		// if(stop){cin.get();}
		// cout<<unitig<<endl;
		// cout<<rangeUnitigs[i].first<<endl;
		//case the rest of the read is too small
		if(readLeft.size()<=unitig.size()){
			uint miss(missmatchNumber(unitig.substr(0,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				// offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig, read.substr(overlap.second,unitig.size()), errors));
			if(miss<miniMiss){
				if(miss<miniMiss){
					kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
					miniMiss=miss;
					miniMissIndice=i;
					nextUnitig=unitig;
					nextOverlap=overlapNum;
				}
			}
		}
	}
	// cout<<"end"<<endl;
	if(miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if (ended){return miniMiss;}
		miniMiss+=mapOnRightEndGreedy(read , path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)}, errors-miniMiss);
	}
	return miniMiss;
}


pair<uint,uint> Aligner::mapOnRight(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,uint start, uint errors){
	string unitig, readLeft(read.substr(overlap.second+k-1)),nextUnitig;
	if(readLeft.empty()){cout<<"should not appears"<<endl;exit(0);return {start,0};}
	auto rangeUnitigs(getBegin(overlap.first));
	uint miniMiss(errors+1),miniMissIndice(9);
	uint next(start);
	kmer nextOverlapNum(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(readLeft.size() <= unitig.size()-k+1){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			if(miss<miniMiss){
				ended=true;
				miniMiss=miss;
				miniMissIndice=i;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
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

//
// pair<uint,uint> Aligner::mapOnRightCache(const string &read, vector<uNumber>& path, const overlapStruct& overlap, const  vector<overlapStruct>& listOverlap, bool& ended,uint start, uint errors){
// 	string unitig, readLeft(read.substr(overlap.pos+k-1)),nextUnitig;
// 	if(readLeft.empty()){cout<<"should not appears"<<endl;exit(0);return {start,0};}
// 	auto rangeUnitigs(overlap.unitig);
// 	uint miniMiss(errors+1),miniMissIndice(9);
// 	uint next(start);
// 	kmer nextOverlapNum(0);
//
// 	for(uint i(0); i<rangeUnitigs.size(); ++i){
// 		unitig=(rangeUnitigs[i]);
// 		//case the rest of the read is too small
// 		if(readLeft.size() <= unitig.size()-k+1){
// 			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
// 			if(miss<miniMiss){
// 				ended=true;
// 				miniMiss=miss;
// 				miniMissIndice=i;
// 			}
// 		}else{
// 			//case the read is big enough we want to recover a true overlap
// 			uint miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
// 			if(miss<miniMiss){
// 				kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
// 				if(miss<miniMiss){
// 					ended=false;
// 					miniMiss=miss;
// 					miniMissIndice=i;
// 					nextOverlapNum=overlapNum;
// 					nextUnitig=unitig;
// 					next=start;
// 					for(uint j(start+1); j<listOverlap.size(); ++j){
// 						if(overlapNum==listOverlap[j].seq and listOverlap[j].pos==overlap.pos+unitig.size()-k+1){
// 							next=j;
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// 	if(ended){
// 		path.push_back(overlap.unitigNumbers[miniMissIndice]);
// 		return {start,miniMiss};
// 	}
// 	if(miniMiss<=errors){
// 		path.push_back(overlap.unitigNumbers[miniMissIndice]);
// 		if(next>start){
// 			return {next,miniMiss};
// 		}
// 		auto res(mapOnRightCache(read , path, {nextOverlapNum,overlap.pos+((uint)nextUnitig.size()-k+1)},listOverlap,ended,start, errors-miniMiss));
// 		return {res.first,res.second+miniMiss};
// 	}
// 	return {start,errors+1};
// }
//
// //
// uint Aligner::checkPair(const pair<kmer, uint>& overlap1, const pair<kmer, uint>& overlap2, const string& read, uNumber& number, uint errorsAllowed){
// 	if(overlap2.second-overlap1.second<k){
// 		int32_t positionget1, positionget2;
// 		auto rangeUnitigs1(getBegin(overlap1.first));
// 		auto rangeUnitigs2(getEnd(overlap2.first));
// 		for(uint i(0); i<rangeUnitigs1.size(); ++i){
// 			positionget1=rangeUnitigs1[i].second;
// 			for(uint j(0); j<rangeUnitigs2.size(); ++j){
// 				positionget2=rangeUnitigs2[j].second;
// 				if(positionget2==positionget1){
// 					number=positionget1;
// 					return 0;
// 				}
// 			}
// 		}
// 		return errorsAllowed+1;
// 	}
// 	string unitig,subRead(read.substr(overlap1.second+k-1,overlap2.second-overlap1.second-(k-1)));
// 	auto rangeUnitigs1(getBegin(overlap1.first));
// 	auto rangeUnitigs2(getEnd(overlap2.first));
// 	uint minMissMatch(errorsAllowed+1),indice(0);
// 	int32_t positionget1, positionget2;
// 	for(uint i(0); i<rangeUnitigs1.size(); ++i){
// 		positionget1=rangeUnitigs1[i].second;
// 		for(uint j(0); j<rangeUnitigs2.size(); ++j){
// 			positionget2=rangeUnitigs2[j].second;
// 			if(positionget2==positionget1){
// 				unitig=getUnitig(rangeUnitigs1[i].second);
// 				if(unitig.size()-2*(k-1)==subRead.size()){
// 					uint missmatch(missmatchNumber(unitig.substr(k-1,subRead.size()), subRead, errorsAllowed));
// 					if(missmatch<minMissMatch){
// 						minMissMatch=missmatch;
// 						indice=i;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	if(minMissMatch<=errorsAllowed){number=(rangeUnitigs1[indice].second);}
// 	return minMissMatch;
// }
//
//
// uint Aligner::checkPairCache(const overlapStruct& overlap1, const overlapStruct& overlap2, const string& read, uNumber& number, uint errorsAllowed){
// 	if(overlap2.pos-overlap1.pos<k){
// 		//TODO maybe it is a bad idea ...
// 		return 0;
// 	}
// 	string unitig,subRead(read.substr(overlap1.pos+k-1,overlap2.pos-overlap1.pos-(k-1)));
// 	auto rangeUnitigs1(overlap1.unitig);
// 	auto rangeUnitigs2(overlap2.unitig);
// 	uint minMissMatch(errorsAllowed+1),indice(0);
// 	string unitig1, unitig2;
// 	for(uint i(0); i<rangeUnitigs1.size(); ++i){
// 		unitig1=rangeUnitigs1[i];
// 		for(uint j(0); j<rangeUnitigs2.size(); ++j){
// 			unitig2=rangeUnitigs2[j];
// 			if(unitig2==unitig1){
// 				if(unitig1.size()-2*(k-1)==subRead.size()){
// 					uint missmatch(missmatchNumber(unitig1.substr(k-1,subRead.size()), subRead, errorsAllowed));
// 					if(missmatch<minMissMatch){
// 						minMissMatch=missmatch;
// 						indice=i;
// 					}
// 				}
// 			}
// 		}
// 	}
// 	if(minMissMatch<=errorsAllowed){number=(overlap1.unitigNumbers[indice]);}
// 	return minMissMatch;
// }
//

uint Aligner::checkBeginGreedy(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	if(overlap.second==0){path.push_back(0);return 0;}
	string readLeft(read.substr(0,overlap.second)),unitig;
	auto rangeUnitigs(getEnd(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(0);
	bool ended(false);
	int offset(0);
	kmer nextOverlap(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			// if(miss==0){
			// 	path.push_back(unitig.size()-readLeft.size()-k+1);
			// 	path.push_back(rangeUnitigs[i].second);
			// 	return 0;
			// }
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			// if(miss==0){
			// 	minMiss+=mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)},errors);
			// 	if(minMiss<=errors){
			// 		path.push_back(rangeUnitigs[indiceMinMiss].second);
			// 		sucessML++;
			// 	}
			// }
			if(miss<minMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
				if(miss<minMiss){
					ended=false;
					minMiss=miss;
					indiceMinMiss=i;
					nextOverlap=overlapNum;
				}
			}

		}
	}

	if(minMiss<=errors){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){
			path.push_back(offset);
			return minMiss;
		}
		minMiss+=mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.second-(rangeUnitigs[indiceMinMiss].first.size()-k+1)},errors-minMiss);
		if(minMiss<=errors){
			++sucessML;
		}
	}
	return minMiss;
}

//
// uint Aligner::checkBeginGreedyCache(const string& read, overlapStruct& overlap, vector<uNumber>& path, uint errors){
// 	if(overlap.pos==0){path.push_back(0);return 0;}
// 	string readLeft(read.substr(0,overlap.pos)),unitig,nextUnitig;
// 	auto rangeUnitigs(overlap.unitig);
// 	uint minMiss(errors+1),indiceMinMiss(0);
// 	bool ended(false);
// 	int offset(0);
// 	kmer nextOverlap(0);
// 	for(uint i(0); i<rangeUnitigs.size(); ++i){
// 		unitig=(rangeUnitigs[i]);
// 		if(unitig.size()-k+1>=readLeft.size()){
// 			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
// 			if(miss<minMiss){
// 				minMiss=miss;
// 				indiceMinMiss=i;
// 				ended=true;
// 				offset=unitig.size()-readLeft.size()-k+1;
// 			}
// 		}else{
// 			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
// 			if(miss<minMiss){
// 				kmer overlapNum(str2num(unitig.substr(0,k-1)));
// 				if(miss<minMiss){
// 					ended=false;
// 					minMiss=miss;
// 					indiceMinMiss=i;
// 					nextOverlap=overlapNum;
// 					nextUnitig=unitig;
// 				}
// 			}
//
// 		}
// 	}
// 	if(minMiss<=errors){
// 		if(ended){
// 			path.push_back(offset);
// 			path.push_back(overlap.unitigNumbers[indiceMinMiss]);
// 			return minMiss;
// 		}
// 		minMiss+=mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.pos-(nextUnitig.size()-k+1)},errors-minMiss);
// 		if(minMiss<=errors){
// 			path.push_back(overlap.unitigNumbers[indiceMinMiss]);
// 			sucessML++;
// 		}
// 	}
// 	return minMiss;
// }


uint Aligner::checkEndGreedy(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	string readLeft(read.substr(overlap.second+k-1)),unitig,nextUnitig;
	// if(readLeft.size()<=k-1){return 0;}
	auto rangeUnitigs(getBegin(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
	kmer nextOverlap(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss<minMiss){
				if(miss<minMiss){
					kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
					minMiss=miss;
					indiceMinMiss=i;
					nextOverlap=overlapNum;
					nextUnitig=unitig;
					// ended=false;
				}
			}
		}
	}

	if(minMiss<=errors){
		// cout<<"hey"<<endl;
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){return minMiss;}
		minMiss+=mapOnRightEndGreedy(read, path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)},errors-minMiss);
		if(minMiss<=errors){++successMR;}
	}
	return minMiss;
}


// uint Aligner::checkEndGreedyCache(const string& read, overlapStruct& overlap, vector<uNumber>& path, uint errors){
// 	string readLeft(read.substr(overlap.pos+k-1)),unitig,nextUnitig;
// 	if(readLeft.size()<=k-1){return 0;}
// 	auto rangeUnitigs(overlap.unitig);
// 	uint minMiss(errors+1),indiceMinMiss(9);
// 	bool ended(false);
// 	kmer nextOverlap(0);
// 	for(uint i(0); i<rangeUnitigs.size(); ++i){
// 		unitig=(rangeUnitigs[i]);
// 		if(unitig.size()-k+1>=readLeft.size()){
// 			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
// 			if(miss<minMiss){
// 				minMiss=miss;
// 				indiceMinMiss=i;
// 				ended=true;
// 			}
// 		}else{
// 			uint miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
// 			if(miss<minMiss){
// 				kmer overlapNum(getRepresentNum(unitig.substr(unitig.size()-k+1,k-1)));
// 				if(miss<minMiss){
// 					minMiss=miss;
// 					indiceMinMiss=i;
// 					nextOverlap=overlapNum;
// 					nextUnitig=unitig;
// 				}
// 			}
// 		}
// 	}
//
// 	if(minMiss<=errors){
// 		path.push_back(overlap.unitigNumbers[indiceMinMiss]);
// 		if(ended){
// 			return minMiss;
// 		}
// 		minMiss+=mapOnRightEndGreedy(read, path, {nextOverlap,overlap.pos+(nextUnitig.size()-k+1)},errors-minMiss);
// 		if(minMiss<=errors){
// 			successMR++;
// 		}
// 	}
// 	return minMiss;
// }


void Aligner::alignPartGreedy(){
	vector<pair<string,string>> multiread;
	vector<uNumber> path;
	string read,header,corrected;
	bool overlapFound(false);
	while(!readFile.eof()){
		readMutex.lock();
		{
			getReads(multiread,10000);
		}
		readMutex.unlock();
		for(uint i(0);i<multiread.size();++i){
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
					if(correctionMode){
						corrected=(recoverPath(path,read.size()));
						// pathFile<<header<<'\n';
						pathFile<<corrected<<'\n';
					}else{
						// pathFile<<header<<'\n';
						// cout<<header<<'\n';
						printPath(path,&pathFile);
					}
				}
				pathMutex.unlock();
			}else{
				if(false){
					noOverlapMutex.lock();
					{
						noOverlapFile<<header<<'\n'<<read<<'\n';
						// cout<<header<<'\n'<<read<<'\n';
					}
					noOverlapMutex.unlock();
				}else{
					notMappedMutex.lock();
					{
						notMappedFile<<header<<'\n'<<read<<'\n';
						// notMappedFile<<readNumber<<'\n';
						// cout<<header<<'\n'<<read<<'\n';
					}
					notMappedMutex.unlock();
				}
			}
		}
		// if(iter++%10==0){
		// 	cout<<"Read : "<<readNumber<<endl;
		// 	cout<<"No Overlap : "<<noOverlapRead<<" Percent : "<<(100*float(noOverlapRead))/readNumber<<endl;
		// 	cout<<"Got Overlap : "<<alignedRead+notAligned<<" Percent : "<<(100*float(alignedRead+notAligned))/readNumber<<endl;
		// 	cout<<"Overlap and Aligned : "<<alignedRead<<" Percent : "<<(100*float(alignedRead))/(alignedRead+notAligned)<<endl;
		// 	cout<<"Overlap but no aligne: "<<notAligned<<" Percent : "<<(100*float(notAligned))/(alignedRead+notAligned)<<endl;
		// 	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
		// 	cout<<"Reads/seconds : "<<readNumber/(chrono::duration_cast<chrono::seconds>(waitedFor).count()+1)<<endl;
		// 	// cout<<"Overlap per reads : "<<(overlaps)/(alignedRead+notAligned)<<endl;
		// 	cout<<endl;
		// }
	}
}
