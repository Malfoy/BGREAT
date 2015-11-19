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
#define kmer uint64_t
//#define kmer __uint128_t

uint64_t transform_to_size_t(__uint128_t n);

//~ namespace std { template <> struct hash<__uint128_t> {
	//~ typedef __uint128_t argument_type;
	//~ typedef uint64_t result_type; uint64_t operator()(__uint128_t key) const { return transform_to_size_t(key); } }; }

class Aligner{
public:
	bool partial,fastq,pathOption;
	ifstream unitigFile, readFile;
	ofstream pathFile, noOverlapFile, notMappedFile;
	atomic<size_t> alignedRead, readNumber, noOverlapRead, notAligned, unitigNumber, overlaps, successMR, successMR2, sucessML,iter;
	uint k;
//	unordered_map<uint64_t,vector<uint32_t>> overlap2unitigs;
//	unordered_map<uint,string> unitigCache;
	unordered_map <kmer,vector<uint32_t>> right;
	unordered_map <kmer,vector<uint32_t>> left;
	vector<string> unitigs;
	kmer offsetUpdate;
	unsigned char coreNumber;
	uint errorsMax,tryNumber;
	mutex unitigMutex, readMutex, indexMutex, pathMutex, noOverlapMutex, notMappedMutex;
	string unitigFileName;
	chrono::system_clock::time_point startChrono;
	bool fullMemory;

	Aligner(const string& Unitigs, const string& paths, const string& noOverlaps, const string& notMapped, uint kValue, unsigned char cores,unsigned int errorsAllowed, bool bpartial,bool bfastq,bool bpath){
		unitigFileName=Unitigs;
		unitigFile.open(unitigFileName);
		pathFile.open(paths);
		noOverlapFile.open(noOverlaps);
		notMappedFile.open(notMapped);
		k=kValue;
		coreNumber=cores;
		errorsMax=errorsAllowed;
		tryNumber=2;
		fullMemory=true;
		partial=bpartial;
		fastq=bfastq;
		alignedRead=readNumber=noOverlapRead=notAligned=unitigNumber=overlaps=successMR=successMR2=sucessML=0;
		offsetUpdate=1;
		offsetUpdate<<=(2*(k-1));
		iter=1;
		pathOption=bpath;
	}

	void indexUnitigs();
	void alignAll(bool, const string&);
	void alignPartGreedy();
	void alignPartExhaustive();
	void indexUnitigsAux();
	vector<pair<kmer,uint>> getListOverlap(const string& read);
	uint8_t checkEndExhaustive(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	bool cover(const string& read, vector<pair<kmer,uint>>& listOverlap, size_t start, size_t end,vector<uNumber>& path);
	bool cover2(const string& read, const vector<pair<kmer,uint>>& listOverlap, const size_t start, size_t end,vector<uNumber>& path);
	size_t mapOnRightAll(const string &read, vector<uNumber>& path, pair<kmer, uint> overlap, const vector<pair<kmer,uint>>& listOverlap, bool& end, size_t start);
	bool mapOnRightEndAll(const string &read, vector<uNumber>& path, pair<kmer, uint> overlap);
	bool mapOnRightEndMax(const string &read, vector<uNumber>& path, pair<kmer, uint> overlap);
	bool mapOnLeftEndAll(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap);
	bool mapOnLeftEndFirst(const string &read, vector<uNumber>& path,const pair<kmer, uint>& overlap);
	bool mapOnLeftEndMax(const string &read, vector<uNumber>& path, const  pair<kmer, uint>& overlap);
	uint8_t mapOnLeftEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& , uint8_t errors);
	uint8_t mapOnRightEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& , uint8_t errors);
	string getUnitig(int position);
	string getUnitigFile(uint position);
	void getReads(vector<pair<string,string>>& reads, uint n);
	vector<uNumber> alignReadGreedy(const string& read, bool& overlapFound, uint8_t errors,bool rc);
	uint8_t checkBeginExhaustive(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	uint8_t checkPair(const pair<kmer, uint>& overlap1, const pair<kmer, uint>& overlap2, const string& read,uNumber& path,uint8_t errorsAllowed);
	uint8_t coverGreedy(const string& read, const vector<pair<kmer,uint>>& listOverlap, const size_t start, size_t end, vector<uNumber>& path, uint8_t  errors, bool&);
	pair<size_t,int> mapOnRight2(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,size_t start, uint8_t errors);
	uint8_t mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors);
	uint8_t checkBeginGreedy(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	uint8_t mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors);
	vector<uNumber> alignReadExhaustive(const string& read, bool& overlapFound, uint8_t errors);
	uint8_t checkEndGreedy(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	string readUnitig(uint position);
	vector<pair<string,uNumber>> getBegin(kmer);
	vector<pair<string,uNumber>> getEnd(kmer);
	string num2str(kmer num);
	kmer str2num(const string& str);
	kmer getRepresentNum(const string& str);
	string recoverPath(vector<uNumber>& numbers,int size);
	vector<uNumber> getBeginNumber(kmer bin);
	vector<uNumber> getEndNumber(kmer bin);
	void updateRC(kmer&	min, char nuc);
	void update(kmer&	min, char nuc);
	pair<size_t,uint8_t> mapOnRight(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,size_t start, uint8_t errors);
	uint8_t mapOnRightEndExhaustivePartial(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors);
	uint8_t mapOnLeftEndExhaustivePartial(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors);
	uint8_t coverGreedyPath(const string& read, const vector<pair<kmer,uint>>& listOverlap, const size_t start, size_t end, vector<uNumber>& path, uint8_t  errors, bool& ended);
	uint8_t checkPairPaths(const pair<kmer, uint>& overlap1, const pair<kmer, uint>& overlap2, const string& read, uNumber& number, uint8_t errorsAllowed, const vector<uNumber> path);
	vector<uNumber> alignReadGreedyPath(const string& read, bool& overlapFound, uint8_t errors, bool rc);
	uint8_t checkBeginExhaustivePath(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	uint8_t checkEndExhaustivePath(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint8_t errors);
	pair<size_t,uint8_t> mapOnRightPath(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,size_t start, uint8_t errors);
	uint8_t mapOnLeftEndExhaustivePaths(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors);
	uint8_t mapOnRightEndExhaustivePath(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint8_t errors);

};

#endif /* defined(__BGREAT__Aligner__) */
