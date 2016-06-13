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


vector<uNumber> Aligner::alignReadGreedy(const string& read, bool& overlapFound, uint errors, bool& rc){
	vector<pair<kmer,uint>> listOverlap(getNOverlap(read,tryNumber));
	if(listOverlap.empty()){++noOverlapRead;return {};}
	overlapFound=true;
	vector<uNumber> pathBegin,pathEnd;
	for(uint start(0); start<(uint)listOverlap.size(); ++start){
		pathBegin={};
		uint errorBegin(checkBeginGreedy(read,listOverlap[start],pathBegin,errors));
		if(errorBegin<=errors){
			pathEnd={};
			uint errorsEnd(checkEndGreedy(read,listOverlap[start],pathEnd,errors-errorBegin));
			if(errorsEnd+errorBegin<=errors){
				++alignedRead;
				reverse(pathBegin.begin(),pathBegin.end());
				pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
				return pathBegin;
			}
		}
	}
	if(!rc){rc=true;return alignReadGreedy(reverseComplements(read), overlapFound,errors, rc);}
	++notAligned;
	return {};
}


vector<uNumber> Aligner::alignReadGreedyAnchors(const string& read, bool& overlapFound, uint errorMax, bool& rc){
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	if(listAnchors.empty()){++noOverlapRead;return {};}
	//~ cout<<"go?"<<endl;
	overlapFound=true;
	vector<uNumber> pathBegin,pathEnd;
	string unitig;
	bool returned(false);
	for(uint start(0); start<(uint)listAnchors.size(); ++start){
		uint unitigNumber(listAnchors[start].first.first),positionUnitig(listAnchors[start].first.second),positionRead(listAnchors[start].second);
		unitig=(unitigs[unitigNumber]);
		if(unitig.size()<k){
			//TODO wtf happen here ?
			continue;
		}
		if(str2num(unitig.substr(positionUnitig,k))!=str2num(read.substr(positionRead,k))){
			//TODO IS THIS USEFULL AND HOW TO IMPROVE
			unitig=reverseComplements(unitig);
			positionUnitig=unitig.size()-k-positionUnitig;
			returned=true;
		}
		if(positionRead>=positionUnitig){
			if(read.size()-positionRead>=unitig.size()-positionUnitig){
				
				//CASE 1 : unitig included in read
				//~ cout<<"1:"<<endl;
				uint errors(missmatchNumber(read.substr(positionRead-positionUnitig,unitig.size()),unitig,errorMax));
				if(errors<=errorMax){
					pathBegin={};
					uint errorBegin(checkBeginGreedy(read,{str2num(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
					if(errorBegin+errors<=errorMax){
						if(returned){
							pathEnd={-(int)unitigNumber};
						}else{
							pathEnd={(int)unitigNumber};	
						}
						uint errorsEnd(checkEndGreedy(read,{str2num(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()-k+1},pathEnd,errorMax-errors-errorBegin));
						if(errorBegin+errors+errorsEnd<=errorMax){	
							++alignedRead;
							reverse(pathBegin.begin(),pathBegin.end());
							pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
							return pathBegin;
						}
					}
				}
			}else{
				
				//CASE 2 : unitig overap read
				//~ cout<<"2:"<<endl;
				uint errors(missmatchNumber(read.substr(positionRead-positionUnitig),unitig.substr(0,read.size()-positionRead+positionUnitig),errorMax));
				if(errors<=errorMax){
					pathBegin={};
					uint errorBegin(checkBeginGreedy(read,{str2num(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
					if(errorBegin+errors<=errorMax){
						++alignedRead;
						reverse(pathBegin.begin(),pathBegin.end());
						if(returned){
							pathBegin.push_back(-(int)unitigNumber);
						}else{
							pathBegin.push_back((int)unitigNumber);
						}
						return pathBegin;
					}
				}
			}
		}else{
			if(read.size()-positionRead>=unitig.size()-positionUnitig){
				
				//CASE 3 : read overlap unitig
				//~ cout<<"3:"<<endl;
				uint errors(missmatchNumber(unitig.substr(positionUnitig-positionRead),read.substr(0,unitig.size()+positionRead-positionUnitig),errorMax));
				if(errors<=errorMax){
					if(returned){
						pathEnd={(int)positionUnitig-(int)positionRead,-(int)unitigNumber};
					}else{
						pathEnd={(int)positionUnitig-(int)positionRead,(int)unitigNumber};
					}//TODO one line
					uint errorsEnd(checkEndGreedy(read,{str2num(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()-k+1},pathEnd,errorMax-errors));
					if(errors+errorsEnd<=errorMax){
						++alignedRead;
						return pathEnd;
					}
				}
			}else{
				
				//CASE 4 : read included in unitig
				//~ cout<<"4:"<<endl;				
				uint errors(missmatchNumber(unitig.substr(positionUnitig-positionRead,read.size()),read,errorMax));
				if(errors<=errorMax){
					++alignedRead;
					if(returned){
						return {(int)positionUnitig-(int)positionRead,-(int)unitigNumber};
					}//TODO one line
					return {(int)positionUnitig-(int)positionRead,(int)unitigNumber};
				}
			}
		}
	}
	if(!rc){rc=true;return alignReadGreedyAnchors(reverseComplements(read), overlapFound,errorMax, rc);}
	++notAligned;
	return {};
}


uint Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	if(overlap.second==0){path.push_back(0);return 0;}
	string unitig,readLeft(read.substr(0,overlap.second)),nextUnitig;
	auto rangeUnitigs(getEnd(overlap.first));
	uint miniMiss(errors+1),miniMissIndice(9);
	bool ended(false);
	int offset(0);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				path.push_back(unitig.size()-readLeft.size()-k+1);
				return 0;
			}else if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read , path, {str2num(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)}, errors);
			}else if(miss<miniMiss){
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
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){
			path.push_back(offset);
			return miniMiss;
		}
		return miniMiss+mapOnLeftEndGreedy(read , path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, errors-miniMiss);
	}
	return miniMiss;
}


uint Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	string unitig,readLeft(read.substr(overlap.second)),nextUnitig;
	if(readLeft.size()<k){return 0;}
	auto rangeUnitigs(getBegin(overlap.first));
	uint miniMiss(errors+1), miniMissIndice(9);
	bool ended(false);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=(rangeUnitigs[i].first);
		//case the rest of the read is too small
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(0,readLeft.size()), readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}else if(miss<miniMiss){
				miniMiss=miss;
				miniMissIndice=i;
				ended=true;
			}
		}else{
			//case the read is big enough we want to recover a true overlap
			uint miss(missmatchNumber(unitig, read.substr(overlap.second,unitig.size()), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return (mapOnRightEndGreedy(read , path, {str2num(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)}, errors));
			}else if(miss<miniMiss){
				if(miss<miniMiss){
					ended=false;
					kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
					miniMiss=miss;
					miniMissIndice=i;
					nextUnitig=unitig;
					nextOverlap=overlapNum;
				}
			}
		}
	}
	if(miniMiss<=errors){
		path.push_back(rangeUnitigs[miniMissIndice].second);
		if(ended){return miniMiss;}
		return (miniMiss+mapOnRightEndGreedy(read , path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)}, errors-miniMiss));
	}
	return miniMiss;
}


uint Aligner::checkBeginGreedy(const string& read,const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	if(overlap.second==0){path.push_back(0);return 0;}
	string readLeft(read.substr(0,overlap.second)),unitig,nextUnitig;
	auto rangeUnitigs(getEnd(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
	int offset(0);
	kmer nextOverlap(0);

	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=rangeUnitigs[i].first;
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				path.push_back(unitig.size()-readLeft.size()-k+1);
				return 0;
			}else if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnLeftEndGreedy(read, path, {str2num(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)},errors);
			}else if(miss<minMiss){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
				if(miss<minMiss){
					ended=false;
					minMiss=miss;
					indiceMinMiss=i;
					nextUnitig=unitig;
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
		return minMiss+mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)},errors-minMiss);
	}
	return minMiss;
}


uint Aligner::checkEndGreedy(const string& read,const pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors){
	string readLeft(read.substr(overlap.second+k-1)),unitig,nextUnitig;
	if(readLeft.empty()){return 0;}
	auto rangeUnitigs(getBegin(overlap.first));
	uint minMiss(errors+1),indiceMinMiss(9);
	bool ended(false);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		unitig=rangeUnitigs[i].first;
		if(unitig.size()-k+1>=readLeft.size()){
			uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return 0;
			}else if(miss<minMiss){
				minMiss=miss;
				indiceMinMiss=i;
				ended=true;
			}
		}else{
			uint miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			if(miss==0){
				path.push_back(rangeUnitigs[i].second);
				return mapOnRightEndGreedy(read, path, {str2num(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)},errors);
			}else if(miss<minMiss){
				if(miss<minMiss){
					kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
					minMiss=miss;
					indiceMinMiss=i;
					nextOverlap=overlapNum;
					nextUnitig=unitig;
					ended=false;
				}
			}
		}
	}
	if(minMiss<=errors){
		path.push_back(rangeUnitigs[indiceMinMiss].second);
		if(ended){return minMiss;}
		return minMiss+mapOnRightEndGreedy(read, path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)},errors-minMiss);
	}
	return minMiss;
}


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
			++readNumber;
			bool rc=false;
			if(pathOption){
				path=alignReadGreedyPath(read,overlapFound,errorsMax,rc);
			}else{
				if(dogMode){
					path=alignReadGreedyAnchors(read,overlapFound,errorsMax,rc);
				}else{
					path=alignReadGreedy(read,overlapFound,errorsMax,rc);
				}
			}
			if(path.size()!=0){
				if(correctionMode){
					corrected=(recoverPath(path,read.size()));
					if(rc){
						corrected=reverseComplements(corrected);
					}
					header+='\n'+corrected+'\n';
					pathMutex.lock();
					{
						fwrite((header).c_str(), sizeof(char), header.size(), pathFilef);
					}
					pathMutex.unlock();
			}else{
					header+='\n'+printPath(path);;
					pathMutex.lock();
					{
						fwrite((header).c_str(), sizeof(char), header.size(), pathFilef);
					}
					pathMutex.unlock();
				}
			}else{
				if(false){
					noOverlapMutex.lock();
					{
						noOverlapFile<<header<<'\n'<<read<<'\n';
					}
					noOverlapMutex.unlock();
				}else{
					header+='\n'+read+'\n';
					notMappedMutex.lock();
					{
						fwrite(header.c_str(), sizeof(char), header.size(), notMappedFilef);
					}
					notMappedMutex.unlock();
				}
			}
		}
	}
}
