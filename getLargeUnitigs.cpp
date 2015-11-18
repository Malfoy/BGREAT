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
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>



using namespace std;



void getBigUnitigs(const string input, const string output, int length){
	ifstream in(input.c_str(),ios::in);
	ofstream out(output.c_str(),ios::trunc);
	string line;

	out<<">N"<<endl;

	while(!in.eof()){
		getline(in,line);
		if(line.size()>length){
			line=line.substr(0,line.size()-1);
			transform(line.begin(),line.end(), line.begin(), ::toupper);
			out<<line;
		}
	}
	out<<endl;
}


int main(int argc, char ** argv)
{
	if(argc<4){
		cout<<"Need in, out and length..."<<endl;
	}
	if(argc==4){
		string arg1=argv[1];
		string arg2=argv[2];
		string arg3=argv[3];
		int l(atoi(arg3.c_str()));
		getBigUnitigs(arg1,arg2,l);
	}
}
