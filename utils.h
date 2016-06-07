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
#ifndef BGREAT_UTILS
#define BGREAT_UTILS


#define uNumber int32_t
#define kmer uint64_t
// #define kmer __uint128_t


using namespace std;


struct overlapStruct{
    kmer seq;
    uint pos;
    vector<string> unitig;
    vector<uint32_t> unitigNumbers;
};


static char char2int[85];


uint64_t transform_to_size_t(__uint128_t& n);
void printPath(const vector<int32_t>& path, ofstream* file);
char revCompChar(char s);
char randNuc();
string mutate(string& read, int n);
kmer str2num(const string& str);
kmer nuc2int(char c);
kmer nuc2intrc(char c);
uint missmatchNumber(const string& seq1, const string& seq2, unsigned int n);
string compactionEnd(const string& seq1,const string& seq2, uint k);
kmer rcb(kmer min,uint n);
string reverseComplements(const string& s);
string getRepresent(const string& s);
void initRc();




#endif /* defined(__UTILS__) */
