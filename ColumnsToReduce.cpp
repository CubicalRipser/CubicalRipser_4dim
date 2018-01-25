//ColumnsToReduce.cpp


#include <algorithm>
#include <unordered_map>


#include "DenseCubicalGrids.h"
#include "BirthdayIndex.h"
#include "ColumnsToReduce.h"

using namespace std;


// vector<BirthdayIndex> columns_to_reduce;
// int dim;
// int max_of_index;

ColumnsToReduce::ColumnsToReduce(DenseCubicalGrids* _dcg) {
		dim = 0;
		int ax = _dcg -> ax;
		int ay = _dcg -> ay;
		int az = _dcg -> az;
		int aw = _dcg -> aw;
		max_of_index = 128*128*128*(aw + 2);
		int index;
		double birthday;
		
		// link_findのときは不要かも...
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
		sort(columns_to_reduce.rbegin(), columns_to_reduce.rend(), BirthdayIndexInverseComparator());
}

int ColumnsToReduce::size() {
	return columns_to_reduce.size();
}
