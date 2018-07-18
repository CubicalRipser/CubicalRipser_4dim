/* columns_to_reduce.cpp

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
#include <unordered_map>


#include "dense_cubical_grids.h"
#include "birthday_index.h"
#include "columns_to_reduce.h"

using namespace std;

ColumnsToReduce::ColumnsToReduce(DenseCubicalGrids* _dcg) {
		dim = 0;
		int ax = _dcg -> ax;
		int ay = _dcg -> ay;
		int az = _dcg -> az;
		int aw = _dcg -> aw;
		max_of_index = 128*128*128*(aw + 2);
		int index;
		double birthday;
		
		for(int w = aw; w > 0; --w){
			for(int z = az; z > 0; --z){
				for (int y = ay; y > 0; --y) {
					for (int x = ax; x > 0; --x) {
						birthday = _dcg -> dense4[x][y][z][w];
						index = x | (y << 7) | (z << 14) | (w << 21);
						if (birthday != _dcg -> threshold) {
							columns_to_reduce.push_back(BirthdayIndex(birthday, index, 0));
						}
					}
				}
			}
		}
		sort(columns_to_reduce.begin(), columns_to_reduce.end(), BirthdayIndexComparator());
}

int ColumnsToReduce::size() {
	return columns_to_reduce.size();
}
