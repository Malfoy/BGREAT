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


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "aligner.h"



using namespace std;

int main(int argc, char ** argv){
	// initRc();
	string reads;
	string unitigs("unitig.fa");
	string pathFile("paths");
	string noOverlapFile("noOverlap.fa");
	string notAlignedFile("notAligned.fa");
	int errors(2);
	int threads(1);
	int ka(31);
	int c;
	bool brute(false),incomplete(false),fastq(false),pathOption(false),correctionMode(false);
	while ((c = getopt (argc, argv, "r:k:g:m:t:f:o:a:biqpc")) != -1){
	switch(c){
		case 'r':
			reads=optarg;
			break;
			case 'k':
				ka=stoi(optarg);
			break;
			case 'g':
				unitigs=(optarg);
			break;
			case 'm':
				errors=stoi(optarg);
			break;
			case 't':
				threads=stoi(optarg);
			break;
			case 'f':
				pathFile=(optarg);
			break;
			case 'o':
				noOverlapFile=(optarg);
			break;
			case 'a':
				notAlignedFile=(optarg);
			break;
			case 'b':
				brute=(true);
			break;
			case 'i':
				incomplete=(true);
			break;
			case 'q':
				fastq=(true);
			break;
			case 'c':
				correctionMode=true;
				// pathFile="corrected.fa";
			break;
		}
	}
	 if(reads!=""){
               	Aligner supervisor(unitigs,pathFile,noOverlapFile,notAlignedFile,ka,threads,errors,incomplete,fastq,pathOption,correctionMode);
               	supervisor.indexUnitigs();
				// supervisor.knowNeighbour();
               	supervisor.alignAll(!brute,reads);
       	}else{
				cout<<"-r read_file"<<endl
				<<"-k k_value (31)"<<endl
              	<<"-g unitig_file (unitig.dot)"<<endl
              	<<"-m n_missmatch (2)"<<endl
              	<<"-t n_thread (1)"<<endl
              	<<"-f path_file (paths)"<<endl
                <<"-o no_overlap_file (noOverlap.fa)"<<endl // Depreciated ??
              	<<"-a not_aligned_file (notAligned.fa)"<<endl
              	<<"-p to align on paths instead of walks"<<endl
              	<<"-q for fastq read file"<<endl;
                <<"-c for correction mode"<<endl;
        }
}
