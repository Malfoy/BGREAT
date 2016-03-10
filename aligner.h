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
#include "utils.h"
#include "BooPHF.h"


using namespace std;


typedef boomphf::SingleHashFunctor<u_int64_t>  hasher;
typedef boomphf::mphf<  u_int64_t, hasher  > MPHF;


namespace std { template <> struct hash<__uint128_t> {
	typedef __uint128_t argument_type;
	typedef uint64_t result_type; uint64_t operator()(__uint128_t key) const { return transform_to_size_t(key); } };
}

struct unitigIndices{
	kmer overlap;
	uint32_t indice1;
	uint32_t indice2;
	uint32_t indice3;
	uint32_t indice4;
};


class Aligner{
public:
	bool partial,fastq,pathOption;
	ifstream unitigFile, readFile;
	ofstream pathFile, noOverlapFile, notMappedFile;
	FILE * pathFilef;
	FILE * notMappedFilef;
	MPHF leftMPHF,rightMPHF;
	vector<unitigIndices> leftIndices,rightIndices;
	atomic<uint> alignedRead, readNumber, noOverlapRead, notAligned, unitigNumber, overlaps,iter;
	uint k;
	// unordered_map <kmer,vector<uint32_t>> right;
	// unordered_map <kmer,vector<uint32_t>> left;
	vector<string> unitigs;
	kmer offsetUpdate;
	uint coreNumber;
	uint gammaFactor;
	uint errorsMax,tryNumber;
	mutex unitigMutex,unitigMutex2, readMutex, indexMutex, pathMutex, noOverlapMutex, notMappedMutex;
	string unitigFileName,pathToWrite;
	chrono::system_clock::time_point startChrono;
	bool fullMemory,correctionMode;

	Aligner(const string& Unitigs, const string& paths, const string& noOverlaps, const string& notMapped, uint kValue, unsigned char cores,unsigned int errorsAllowed, bool bpartial,bool bfastq,bool bpath,bool bcorrectionMode){
		unitigFileName=Unitigs;
		unitigFile.open(unitigFileName);
		// pathFile.open(paths);
		pathFilef=fopen(paths.c_str(),"wb");
		notMappedFilef=fopen(notMapped.c_str(),"wb");
		// noOverlapFile.open(noOverlaps);
		// notMappedFile.open(notMapped);
		// notMappedFile.open(notMapped,ios::binary);
		k=kValue;
		coreNumber=cores;
		errorsMax=errorsAllowed;
		tryNumber=2;
		gammaFactor=2;
		fullMemory=true;
		correctionMode=bcorrectionMode;
		partial=bpartial;
		fastq=bfastq;
		alignedRead=readNumber=noOverlapRead=notAligned=unitigNumber=overlaps=0;
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
	uint checkEndExhaustive(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	bool cover(const string& read, vector<pair<kmer,uint>>& listOverlap, uint start, uint end,vector<uNumber>& path);
	bool cover2(const string& read, const vector<pair<kmer,uint>>& listOverlap, const uint start, uint end,vector<uNumber>& path);
	uint mapOnRightAll(const string &read, vector<uNumber>& path, pair<kmer, uint> overlap, const vector<pair<kmer,uint>>& listOverlap, bool& end, uint start);
	bool mapOnRightEndAll(const string &read, vector<uNumber>& path, pair<kmer, uint> overlap);
	bool mapOnRightEndMax(const string &read, vector<uNumber>& path, pair<kmer, uint> overlap);
	bool mapOnLeftEndAll(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap);
	bool mapOnLeftEndFirst(const string &read, vector<uNumber>& path,const pair<kmer, uint>& overlap);
	bool mapOnLeftEndMax(const string &read, vector<uNumber>& path, const  pair<kmer, uint>& overlap);
	uint mapOnLeftEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& , uint errors);
	uint mapOnRightEndExhaustive(const string &read, vector<uNumber>& path, const pair<kmer, uint>& , uint errors);
	string getUnitig(int position);
	string getUnitigFile(uint position);
	void getReads(vector<pair<string,string>>& reads, uint n);
	vector<uNumber> alignReadGreedy(const string& read, bool& overlapFound, uint errors,bool rc);
	uint checkBeginExhaustive(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	uint checkPair(const pair<kmer, uint>& overlap1, const pair<kmer, uint>& overlap2, const string& read,uNumber& path,uint errorsAllowed);
	uint coverGreedy(const string& read, const vector<pair<kmer,uint>>& listOverlap, const uint start, uint end, vector<uNumber>& path, uint  errors, bool&);
	pair<uint,int> mapOnRight2(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,uint start, uint errors);
	uint mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	uint checkBeginGreedy(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	uint mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	vector<uNumber> alignReadExhaustive(const string& read, bool& overlapFound, uint errors);
	uint checkEndGreedy(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	string readUnitig(uint position);
	vector<pair<string,uNumber>> getBegin(kmer);
	vector<pair<string,uNumber>> getEnd(kmer);
	string num2str(kmer num);
	kmer getRepresentNum(const string& str);
	string recoverPath(vector<uNumber>& numbers,uint size);
	vector<uNumber> getBeginNumber(kmer bin);
	vector<uNumber> getEndNumber(kmer bin);
	void updateRC(kmer&	min, char nuc);
	void update(kmer&	min, char nuc);
	pair<uint,uint> mapOnRight(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,uint start, uint errors);
	uint mapOnRightEndExhaustivePartial(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	uint mapOnLeftEndExhaustivePartial(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	uint coverGreedyPath(const string& read, const vector<pair<kmer,uint>>& listOverlap, const uint start, uint end, vector<uNumber>& path, uint  errors, bool& ended);
	uint checkPairPaths(const pair<kmer, uint>& overlap1, const pair<kmer, uint>& overlap2, const string& read, uNumber& number, uint errorsAllowed, const vector<uNumber> path);
	vector<uNumber> alignReadGreedyPath(const string& read, bool& overlapFound, uint errors, bool rc);
	uint checkBeginExhaustivePath(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	uint checkEndExhaustivePath(const string& read, pair<kmer, uint>& overlap, vector<uNumber>& path, uint errors);
	pair<uint,uint> mapOnRightPath(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap, const  vector<pair<kmer,uint>>& listOverlap, bool& ended,uint start, uint errors);
	uint mapOnLeftEndExhaustivePaths(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	uint mapOnRightEndExhaustivePath(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors);
	vector<uNumber> alignReadExhaustivePath(const string& read, bool& overlapFound, uint errors);
	vector<overlapStruct> getListOverlapStuct(const string& read);
	vector<uNumber> alignReadGreedyCache(const string& read, bool& overlapFound, uint errors, bool rc);
	uint coverGreedyCache(const string& read, const vector<overlapStruct>& listOverlap, const uint start, uint end, vector<uNumber>& path, uint  errors, bool& ended);
	vector<overlapStruct> getListOverlapCache(const string& read);
	uint checkEndGreedyCache(const string& read, overlapStruct& overlap, vector<uNumber>& path, uint errors);
	uint checkBeginGreedyCache(const string& read, overlapStruct& overlap, vector<uNumber>& path, uint errors);
	uint checkPairCache(const overlapStruct& overlap1, const overlapStruct& overlap2, const string& read, uNumber& number, uint errorsAllowed);
	pair<uint,uint> mapOnRightCache(const string &read, vector<uNumber>& path, const overlapStruct& overlap, const  vector<overlapStruct>& listOverlap, bool& ended,uint start, uint errors);
	vector<pair<kmer,uint>> getNOverlap(const string& read, uint n);
	vector<uNumber> alignReadExhaustiveR(const string& read, bool& overlapFound, uint errors);
	vector<uNumber> alignReadGreedyCover(const string& read, bool& overlapFound, uint errors, bool rc);
	vector<string> getEndStr(kmer bin);
	vector<string> getBeginStr(kmer bin);
	void knowNeighbour();
	vector<pair<string,uNumber>> getBeginOpti(kmer bin, uNumber last);
	vector<pair<string,uNumber>> getEndOpti(kmer bin, uNumber last);
	string printPath(const vector<int32_t>& path);

};

#endif /* defined(__BGREAT__Aligner__) */
