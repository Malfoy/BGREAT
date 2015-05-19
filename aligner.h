
//  Aligner.h
//  BGREAT
//
//  Created by malfoy on 09/04/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

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
	int errorsMax,tryNumber;
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
	int checkEndExhaustive(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, int errors);
	bool cover(const string& read, vector<pair<uint64_t,uint>>& listOverlap, size_t start, size_t end,vector<uNumber>& path);
	bool cover2(const string& read, const vector<pair<uint64_t,uint>>& listOverlap, const size_t start, size_t end,vector<uNumber>& path);
	size_t mapOnRightAll(const string &read, vector<uNumber>& path, pair<uint64_t, uint> overlap, const vector<pair<uint64_t,uint>>& listOverlap, bool& end, size_t start);
	bool mapOnRightEndAll(const string &read, vector<uNumber>& path, pair<uint64_t, uint> overlap);
	bool mapOnRightEndMax(const string &read, vector<uNumber>& path, pair<uint64_t, uint> overlap);
	bool mapOnLeftEndAll(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap);
	bool mapOnLeftEndFirst(const string &read, vector<uNumber>& path,const pair<uint64_t, uint>& overlap);
	bool mapOnLeftEndMax(const string &read, vector<uNumber>& path, const  pair<uint64_t, uint>& overlap);
	int mapOnLeftEndExhaustive(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& , int errors);
	int mapOnRightEndExhaustive(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& , int errors);
	string getUnitig(int position);
	string getUnitig2(uint position);
	void getReads(vector<string>& reads,uint);
	vector<uNumber> alignRead3(const string& read, bool& overlapFound, int errors,bool rc);
	int checkBeginExhaustive(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, int errors);
	int checkPair2(const pair<uint64_t, uint>& overlap1, const pair<uint64_t, uint>& overlap2, const string& read,uNumber& path,int errorsAllowed);
	int coverGreedy(const string& read, const vector<pair<uint64_t,uint>>& listOverlap, const size_t start, size_t end, vector<uNumber>& path, int  errors, bool&);
	pair<size_t,int> mapOnRight2(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap, const  vector<pair<uint64_t,uint>>& listOverlap, bool& ended,size_t start, int errors);
	int mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap , int errors);
	int checkBeginGreedy(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, int errors);
	int mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<uint64_t, uint>& overlap , int errors);
	vector<uNumber> alignRead4(const string& read, bool& overlapFound, int errors);
	int checkEndGreedy(const string& read, pair<uint64_t, uint>& overlap, vector<uNumber>& path, int errors);
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

};

#endif /* defined(__BGREAT__Aligner__) */
