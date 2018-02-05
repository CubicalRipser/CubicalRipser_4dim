/* DenseCubicalGrids.cpp

Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.

This file is part of CubicalRipser_4dim.

CubicalRipser: C++ system for computation of Cubical persistence pairs
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

CubicalRipser is deeply depending on 'Ripser', software for Vietoris-Rips 
persitence pairs by Ulrich Bauer, 2015-2016.  We appreciate Ulrich very much.
We rearrange his codes of Ripser and add some new ideas for optimization on it 
and modify it for calculation of a Cubical filtration.

This part of CubicalRiper is a calculator of cubical persistence pairs for 
4 dimensional voxel data. The input data format conforms to that of DIPHA or of PERSEUS.
 See more descriptions in README.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include "DenseCubicalGrids.h"


using namespace std;

DenseCubicalGrids::DenseCubicalGrids(const std::string& filename, double _threshold, file_format _format)  {
	
	threshold = _threshold;
	format = _format;

	if(format == DIPHA){// ???.complex, DIPHA format
		std::ifstream reading_file; 

		ifstream fin( filename, ios::in | ios::binary ); 
		int64_t d;

		fin.read( ( char * ) &d, sizeof( int64_t ) ); // magic number
		//assert(d == 8067171840);
		fin.read( ( char * ) &d, sizeof( int64_t ) ); // type number
		//assert(d == 1);
		fin.read( ( char * ) &d, sizeof( int64_t ) ); //data num
		fin.read( ( char * ) &d, sizeof( int64_t ) ); // dim 
		dim = d;
		fin.read( ( char * ) &d, sizeof( int64_t ) );
		ax = d;
		fin.read( ( char * ) &d, sizeof( int64_t ) );
		ay = d;
		fin.read( ( char * ) &d, sizeof( int64_t ) );
		az = d;
		fin.read( ( char * ) &d, sizeof( int64_t ) );
		aw = d;
		//assert(0 < ax && ax < 510 && 0 < ay && ay < 510 && 0 < az && az < 510);
		cout << "ax : ay : az : aw = " << ax << " : " << ay << " : " << az << " : " << aw << endl;

		double dou;
		for(int w = 0; w < aw + 2; ++w){
			for(int z = 0; z < az + 2; ++z){
				for (int y = 0; y < ay + 2; ++y) {
					for (int x = 0; x < ax + 2; ++x) {
						if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az && 0 < w && w <= aw){
							if (!fin.eof()) {
								fin.read( ( char * ) &dou, sizeof( double ) );
								dense4[x][y][z][w] = dou;
							} else {
								cout << "file endof error " << endl;
							}
						}
						else {
							dense4[x][y][z][w] = threshold;
						}
					}
				}
			}
		}
		fin.close();
	}  else if(format == PERSEUS){// PERCEUS format

		std::ifstream reading_file; 
		reading_file.open(filename.c_str(), std::ios::in); 

		std::string reading_line_buffer; 
		std::getline(reading_file, reading_line_buffer); 
		dim = std::atoi(reading_line_buffer.c_str());
		std::getline(reading_file, reading_line_buffer); 
		ax = std::atoi(reading_line_buffer.c_str()); 
		std::getline(reading_file, reading_line_buffer); 
		ay = std::atoi(reading_line_buffer.c_str()); 
		std::getline(reading_file, reading_line_buffer); 
		az = std::atoi(reading_line_buffer.c_str());
		std::getline(reading_file, reading_line_buffer); 
		aw = std::atoi(reading_line_buffer.c_str());

		for(int w = 0; w < aw + 2; ++w){
			for(int z = 0; z < az + 2; ++z){
				for (int y = 0; y <ay + 2; ++y) { 
					for (int x = 0; x < ax + 2; ++x) { 
						if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az && 0 < w && w <= aw){ 
							if (!reading_file.eof()) { 
								std::getline(reading_file, reading_line_buffer); 
								dense4[x][y][z][w] = std::atoi(reading_line_buffer.c_str()); 
								if (dense4[x][y][z][w] == -1) { 
									dense4[x][y][z][w] = threshold; 
								} 
							} 
						}
						else { 
							dense4[x][y][z][w] = threshold; 
						} 
					} 
				}
			}
		}
	} 
}


double DenseCubicalGrids::getBirthday(int index, int dim){
		int cx = index & 0x7f;
		int cy = (index >> 7) & 0x7f;
		int cz = (index >> 14) & 0x7f;
		int cw = (index >> 21) & 0x7f;
		int cm = (index >> 28) & 0x0f;

		switch(dim){
			case 0:
				return dense4[cx][cy][cz][cw];
			case 1:
				switch(cm){
					case 0:
						return max(dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw]);
					case 1:
						return max(dense4[cx][cy][cz][cw], dense4[cx][cy + 1][cz][cw]);
					case 2:
						return max(dense4[cx][cy][cz][cw], dense4[cx][cy][cz + 1][cw]);
					default:
						return max(dense4[cx][cy][cz][cw], dense4[cx][cy][cz][cw + 1]);
					break;
				}
			case 2:
				switch(cm){
					case 0: // x - y (fix z, w)
						return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
						dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw]});
					case 1: // z - x (fix y, w)
						return max({dense4[cx][cy][cz][cw], dense4[cx][cy][cz + 1][cw], 
						dense4[cx + 1][cy][cz + 1][cw], dense4[cx + 1][cy][cz][cw]});
					case 2: // y - z (fix x, w)
						return max({dense4[cx][cy][cz][cw], dense4[cx][cy + 1][cz][cw], 
						dense4[cx][cy + 1][cz + 1][cw], dense4[cx][cy][cz + 1][cw]});
					case 3: // x - w
						return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
						dense4[cx + 1][cy][cz][cw + 1], dense4[cx][cy][cz][cw + 1]});
					case 4: // y - w
						return max({dense4[cx][cy][cz][cw], dense4[cx][cy + 1][cz][cw], 
						dense4[cx][cy + 1][cz][cw + 1], dense4[cx][cy][cz][cw + 1]});
					case 5: // z - w
						return max({dense4[cx][cy][cz][cw], dense4[cx][cy][cz + 1][cw], 
						dense4[cx][cy][cz + 1][cw + 1], dense4[cx][cy][cz][cw + 1]});
				}
			case 3:
				switch(cm){
					case 0: // x - y - z
						return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
									dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw],
									dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], 
									dense4[cx + 1][cy + 1][cz + 1][cw], dense4[cx][cy + 1][cz + 1][cw]
								});
					case 1: // x - y - w
						return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
									dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw],
									dense4[cx][cy][cz][cw + 1], dense4[cx + 1][cy][cz][cw + 1], 
									dense4[cx + 1][cy + 1][cz][cw + 1], dense4[cx][cy + 1][cz][cw + 1]
								});
					case 2: // x - z - w
						return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
									dense4[cx + 1][cy][cz][cw + 1], dense4[cx][cy][cz][cw + 1],
									dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], 
									dense4[cx + 1][cy][cz + 1][cw + 1], dense4[cx][cy][cz + 1][cw + 1]
								});
					case 3: // y - z - w
						return max({dense4[cx][cy][cz][cw], dense4[cx][cy][cz][cw + 1], 
									dense4[cx][cy + 1][cz][cw + 1], dense4[cx][cy + 1][cz][cw],
									dense4[cx][cy][cz + 1][cw], dense4[cx][cy][cz + 1][cw + 1], 
									dense4[cx][cy + 1][cz + 1][cw + 1], dense4[cx][cy + 1][cz + 1][cw]
								});
				}
			case 4:
				return max({dense4[cx][cy][cz][cw], dense4[cx + 1][cy][cz][cw], 
							dense4[cx + 1][cy + 1][cz][cw], dense4[cx][cy + 1][cz][cw],
							dense4[cx][cy][cz + 1][cw], dense4[cx + 1][cy][cz + 1][cw], 
							dense4[cx + 1][cy + 1][cz + 1][cw], dense4[cx][cy + 1][cz + 1][cw],
							dense4[cx][cy][cz][cw + 1], dense4[cx + 1][cy][cz][cw + 1], 
							dense4[cx + 1][cy + 1][cz][cw + 1], dense4[cx][cy + 1][cz][cw + 1],
							dense4[cx][cy][cz + 1][cw + 1], dense4[cx + 1][cy][cz + 1][cw + 1], 
							dense4[cx + 1][cy + 1][cz + 1][cw + 1], dense4[cx][cy + 1][cz + 1][cw + 1]
						});
		}
		return threshold;
}


