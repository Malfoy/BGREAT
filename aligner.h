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


#ifndef __BGREAT__Aligner__
#define __BGREAT__Aligner__

#include <stdio.h>
#include <fstream>
#include <vector>
#include <atomic>
#include <mutex>
#include <unordered_map>

//#include <sparsehash/sparse_hash_map>

using namespace std;
//using namespace google;

#define uNumber int32_t

class Aligner{
public:
	ifstream unitigFile, readFile;
	ofstream pathFile, noOverlapFile, notMappedFile;
	atomic<size_t> alignedRead, readNumber, noOverlapRead, notAligned, unitigNumber, overlaps, successMR, successMR2, sucessML,iter;
	uint k;
//	unordered_map<uint64_t,vector<uint32_t>> overlap2unitigs;
//	unordered_map<uint,string> unitigCache;
	unordered_map <uint64_t,vector<uint32_t>> right;
	unordered_map <uint64_t,vector<uint32_t>> left;
	vector<string> unitigs;
	uint64_t offsetUpdate;
	unsigned char coreNumber;
	uint errorsMax,tryNumber;
	mutex unitigMutex, readMutex, indexMutex, pathMutex, noOverlapMutex, notMappedMutex;
	string unitigFileName;
	chrono::system_clock::time_point startChrono;
	bool fullMemory;

	Aligner(const string& reads, const string& Unitigs, const string& paths, const string& noOverlaps, const string& notMapped, uint kValue, unsigned char cores,unsigned int errorsAllowed){
		unitigFileName=Unitigs;
		unitigFile.open(unitigFileName);
		readFile.open(reads);
		pathFile.open(paths);
		noOverlapFile.open(noOverlaps);
		notMappedFile.open(notMapped);
		k=kValue;
		coreNumber=cores;
		errorsMax=errorsAllowed;
		tryNumber=2;
		fullMemory=true;
		alignedRead=readNumber=noOverlapRead=notAligned=unitigNumber=overlaps=successMR=successMR2=sucessML=0;
		offsetUpdate=1;
		offsetUpdate<<=(2*(k-1));
		iter=1;
	}

	void indexUnitigs();
	void alignAll(bool);
	void alignPartGreedy();
	void alignPartExhaustive();
	void indexUnitigsAux();
	vector<pair<uint64_t,uint>> getListOverlap(const string& read);
	uint8_t checkEndExhaustive(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	bool cover(const string& read, vector<pair<uint64_t,uint>>& listOverlap, size_t start, size_t end,vector<uNumber>& path);
	bool cover2(const string& read, const vector<pair<uint64_t,uint>>& listOverlap, const size_t start, size_t end,vector<uNumber>& path);
	size_t mapOnRightAll(const string &read, vector<uNumber>& path, pair<uint64_t, uint> overlap, const vector<pair<uint64_t,uint>>& listOverlap, bool& end, size_t start);
	bool mapOnRightEndAll(const string &read, vector<uNumber>& path, pair<uint64_t, uint> overlap);
	bool mapOnRightEndMax(const string &read, vector<uNumber>& path, pair<uint64_t, uint> overlap);
	bool mapOnLeftEndAll(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap);
	bool mapOnLeftEndFirst(const string &read, vector<uNumber>& path,const pair<uint64_t, uint>& overlap);
	bool mapOnLeftEndMax(const string &read, vector<uNumber>& path, const  pair<uint64_t, uint>& overlap);
	uint8_t mapOnLeftEndExhaustive(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& , uint8_t errors);
	uint8_t mapOnRightEndExhaustive(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& , uint8_t errors);
	string getUnitig(int position);
	string getUnitigFile(uint position);
	void getReads(vector<pair<string,string>>& reads, uint n);
	vector<uNumber> alignReadGreedy(const string& read, bool& overlapFound, uint8_t errors,bool rc);
	uint8_t checkBeginExhaustive(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	uint8_t checkPair(const pair<uint64_t, uint>& overlap1, const pair<uint64_t, uint>& overlap2, const string& read,uNumber& path,uint8_t errorsAllowed);
	uint8_t coverGreedy(const string& read, const vector<pair<uint64_t,uint>>& listOverlap, const size_t start, size_t end, vector<uNumber>& path, uint8_t  errors, bool&);
	pair<size_t,int> mapOnRight2(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap, const  vector<pair<uint64_t,uint>>& listOverlap, bool& ended,size_t start, uint8_t errors);
	uint8_t mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap , uint8_t errors);
	uint8_t checkBeginGreedy(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	uint8_t mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap , uint8_t errors);
	vector<uNumber> alignReadExhaustive(const string& read, bool& overlapFound, uint8_t errors);
	uint8_t checkEndGreedy(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	string readUnitig(uint position);
	vector<pair<string,uNumber>> getBegin(uint64_t);
	vector<pair<string,uNumber>> getEnd(uint64_t);
	string num2str(uint64_t num);
	uint64_t str2num(const string& str);
	uint64_t getRepresentNum(const string& str);
	string recoverPath(vector<uNumber>& numbers,int size);
	vector<uNumber> getBeginNumber(uint64_t bin);
	vector<uNumber> getEndNumber(uint64_t bin);
	void updateRC(uint64_t&	min, char nuc);
	void update(uint64_t&	min, char nuc);
	pair<size_t,uint8_t> mapOnRight(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap, const  vector<pair<uint64_t,uint>>& listOverlap, bool& ended,size_t start, uint8_t errors);

};

#endif /* defined(__BGREAT__Aligner__) */
