//SimplexCoboundaryEnumerator.cpp

#include <iostream>
#include <algorithm>
#include <string>

#include "DenseCubicalGrids.h"
#include "BirthdayIndex.h"
#include "ColumnsToReduce.h"
#include "SimplexCoboundaryEnumerator.h"


SimplexCoboundaryEnumerator::SimplexCoboundaryEnumerator(){
	nextCoface = BirthdayIndex(0, -1, 1);
}


void SimplexCoboundaryEnumerator::setSimplexCoboundaryEnumerator(BirthdayIndex _s, DenseCubicalGrids* _dcg) {
	simplex = _s;
	dcg = _dcg;
	dim = simplex.dim;
	birthtime = simplex.birthday;
	ax = _dcg -> ax;
	ay = _dcg -> ay;
	az = _dcg -> az;
	aw = _dcg -> aw;
	
	cx = simplex.index & 0x7f;
	cy = (simplex.index >> 7) & 0x7f;
	cz = (simplex.index >> 14) & 0x7f;
	cw = (simplex.index >> 21) & 0x7f;
	cm = (simplex.index >> 28) & 0x0f;

	threshold = _dcg->threshold;
	count = 0;
}


bool SimplexCoboundaryEnumerator::hasNextCoface() {
	int index = 0;
	double birthday = 0;
	switch (dim) {
	case 0: // dim0
			for (int i = count; i < 8; ++i) {
				switch (i){
					case 0: // w +
					index = (3 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | cx;
					birthday = max(birthtime, dcg -> dense4[cx][cy][cz][cw + 1]);
					break;

					case 1: // w -
					index = (3 << 28) | ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
					birthday = max(birthtime, dcg -> dense4[cx][cy][cz][cw - 1]);
					break;

					case 2: // z +
					index = (2 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | cx;
					birthday = max(birthtime, dcg -> dense4[cx][cy][cz + 1][cw]);
					break;

					case 3: // z -
					index = (2 << 28) | (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
					birthday = max(birthtime, dcg -> dense4[cx][cy][cz - 1][cw]);
					break;

					case 4: // y +
					index = (1 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | cx;
					birthday = max(birthtime, dcg -> dense4[cx][cy + 1][cz][cw]);
					break;

					case 5: // y -
					index = (1 << 28) | (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
					birthday = max(birthtime, dcg -> dense4[cx][cy - 1][cz][cw]);
					break;

					case 6: // x +
					index = (0 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | cx;
					birthday = max(birthtime, dcg -> dense4[cx + 1][cy][cz][cw]);
					break;

					case 7: // x -
					index = (0 << 28) | (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
					birthday = max(birthtime, dcg -> dense4[cx - 1][cy][cz][cw]);
					break;
				}
				if (birthday != threshold) {
					count = i + 1;
					nextCoface = BirthdayIndex(birthday, index, 1);
					return true;
				}
			}
			return false;
	case 1:// dim1
			switch (cm) {
				case 0: // dim1 type0 (x-axis -> )
				for(int i = count; i < 6; ++i){
					switch(i){
						case 0: // x - w +
						index = (3 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1]});
						break;

						case 1: // x - w -
						index = (3 << 28) | ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1]});
						break;

						case 2: // x - z +
						index = (1 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw]});
						break;

						case 3: // x - z -
						index = (1 << 28) | (cw << 21) |((cz - 1) << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw]});
						break;

						case 4: // x - y +
						index = (0 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw]});
						break;

						case 5: // x - y -
						index = (0 << 28) | (cw << 21) |(cz << 14) | ((cy - 1) << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw]});
						break;
					}

					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 2);
						return true;
					}
				}
				return false;

				case 1: // dim1 type1 (y-axis -> )
				for(int i = count; i < 6; ++i){
					switch(i){
						case 0: // y - w +
						index = (4 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx][cy + 1][cz][cw + 1]});
						break;

						case 1: // y - w -
						index = (4 << 28) | ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx][cy + 1][cz][cw - 1]});
						break;

						case 2: // y - z +
						index = (2 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw]});
						break;

						case 3: // y - z -
						index = (2 << 28) | (cw << 21) |((cz - 1) << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx][cy + 1][cz - 1][cw]});
						break;

						case 4: // y - x +
						index = (0 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw]});
						break;

						case 5: // y - x -
						index = (0 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | (cx - 1);
						birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 2);
						return true;
					}
				}
				return false;

				case 2: // dim1 type2 (z-axis -> )
				for(int i = count; i < 6; ++i){
					switch(i){
						case 0: // z - w +
						index = (5 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx][cy][cz + 1][cw + 1]});
						break;

						case 1: // z - w -
						index = (5 << 28) | ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx][cy][cz + 1][cw - 1]});
						break;

						case 2: // z - y +
						index = (2 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw]});
						break;

						case 3: // z - y -
						index = (2 << 28) | (cw << 21) |(cz << 14) | ((cy - 1) << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx][cy - 1][cz + 1][cw]});
						break;

						case 4: // z - x +
						index = (1 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw]});
						break;

						case 5: // z - x -
						index = (1 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | (cx - 1);
						birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy][cz + 1][cw]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 2);
						return true;
					}
				}
				return false;

				case 3: // dim1 type3 (w-axis -> )
				for(int i = count; i < 6; ++i){
					switch(i){
						case 0: // w - z +
						index = (5 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx][cy][cz + 1][cw + 1]});
						break;

						case 1: // w - z -
						index = (5 << 28) | (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx][cy][cz - 1][cw + 1]});
						break;

						case 2: // w - y +
						index = (4 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx][cy + 1][cz][cw + 1]});
						break;

						case 3: // w - y -
						index = (4 << 28) | (cw << 21) |(cz << 14) | ((cy - 1) << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx][cy - 1][cz][cw + 1]});
						break;

						case 4: // w - x +
						index = (3 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy][cz][cw + 1]});
						break;

						case 5: // w - x -
						index = (3 << 28) | (cw << 21) |(cz << 14) | (cy << 7) | (cx - 1);
						birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy][cz][cw + 1]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 2);
						return true;
					}
				}
				return false;
			}
			return false;
	case 2:// dim2
			switch (cm) {
				case 0: // dim2 type0 (fix x - y)
				for(int i = count; i < 4; ++i){
					switch(i){
						case 0: // w +
						index = (1 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1], 
							dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1]});
						break;

						case 1: // w -
						index = (1 << 28)| ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1], 
							dcg -> dense4[cx][cy + 1][cz][cw - 1],dcg -> dense4[cx + 1][cy + 1][cz][cw - 1]});
						break;

						case 2: // z +
						index = (0 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
							dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw]});
						break;

						case 3: // z -
						index = (0 << 28)| (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw], 
							dcg -> dense4[cx][cy + 1][cz - 1][cw],dcg -> dense4[cx + 1][cy + 1][cz - 1][cw]});
						break;

					}

					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 3);
						return true;
					}
				}
				return false;

				case 1: // dim2 type1 (fix x - z)
				for(int i = count; i < 4; ++i){
					switch(i){
						case 0: // w +
						index = (2 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1], 
							dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1]});
						break;

						case 1: // w -
						index = (2 << 28)| ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1], 
							dcg -> dense4[cx][cy][cz + 1][cw - 1],dcg -> dense4[cx + 1][cy][cz + 1][cw - 1]});
						break;

						case 2: // y +
						index = (0 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
							dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw]});
						break;

						case 3: // y -
						index = (0 << 28)| (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw], 
							dcg -> dense4[cx][cy - 1][cz + 1][cw],dcg -> dense4[cx + 1][cy - 1][cz + 1][cw]});
						break;

					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 3);
						return true;
					}
				}
				return false;

				case 2: // dim2 type3 (fix y - z)
				for(int i = count; i < 4; ++i){
					switch(i){
						case 0: // w +
						index = (3 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx][cy + 1][cz][cw + 1], 
							dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx][cy + 1][cz + 1][cw + 1]});
						break;

						case 1: // w -
						index = (3 << 28)| ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx][cy + 1][cz][cw - 1], 
							dcg -> dense4[cx][cy][cz + 1][cw - 1],dcg -> dense4[cx][cy + 1][cz + 1][cw - 1]});
						break;

						case 2: // x +
						index = (0 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
							dcg -> dense4[cx + 1][cy][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw]});
						break;

						case 3: // x -
						index = (0 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
						birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw], 
							dcg -> dense4[cx - 1][cy][cz + 1][cw],dcg -> dense4[cx - 1][cy + 1][cz + 1][cw]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 3);
						return true;
					}
				}
				return false;

				case 3: // dim2 type2 (fix x - w)
				for(int i = count; i < 4; ++i){
					switch(i){
						case 0: // z +
						index = (2 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
							dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1]});
						break;

						case 1: // z -
						index = (2 << 28)| (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw], 
							dcg -> dense4[cx][cy][cz - 1][cw + 1],dcg -> dense4[cx + 1][cy][cz - 1][cw + 1]});
						break;

						case 2: // y +
						index = (1 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
							dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1]});
						break;

						case 3: // y -
						index = (1 << 28)| (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw], 
							dcg -> dense4[cx][cy - 1][cz][cw + 1],dcg -> dense4[cx + 1][cy - 1][cz][cw + 1]});
						break;

					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 3);
						return true;
					}
				}
				return false;

				case 4: // dim2 type4 (fix y - w)
				for(int i = count; i < 4; ++i){
					switch(i){
						case 0: // z +
						index = (3 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw], 
							dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx][cy + 1][cz + 1][cw + 1]});
						break;

						case 1: // z -
						index = (3 << 28)| (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx][cy + 1][cz - 1][cw], 
							dcg -> dense4[cx][cy][cz - 1][cw + 1],dcg -> dense4[cx][cy + 1][cz - 1][cw + 1]});
						break;

						case 2: // x +
						index = (1 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
							dcg -> dense4[cx + 1][cy][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1]});
						break;

						case 3: // x -
						index = (1 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
						birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw], 
							dcg -> dense4[cx - 1][cy][cz][cw + 1],dcg -> dense4[cx - 1][cy + 1][cz][cw + 1]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 3);
						return true;
					}
				}
				return false;

				case 5: // dim2 type5 (fix z - w)
				for(int i = count; i < 4; ++i){
					switch(i){
						case 0: // y +
						index = (3 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx][cy + 1][cz + 1][cw], 
							dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx][cy + 1][cz + 1][cw + 1]});
						break;

						case 1: // y -
						index = (3 << 28)| (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx][cy - 1][cz + 1][cw], 
							dcg -> dense4[cx][cy - 1][cz][cw + 1],dcg -> dense4[cx][cy - 1][cz + 1][cw + 1]});
						break;

						case 2: // x +
						index = (2 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
							dcg -> dense4[cx + 1][cy][cz][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1]});
						break;

						case 3: // x -
						index = (2 << 28)| (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
						birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy][cz + 1][cw], 
							dcg -> dense4[cx - 1][cy][cz][cw + 1],dcg -> dense4[cx - 1][cy][cz + 1][cw + 1]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 3);
						return true;
					}
				}
				return false;
			}
			return false;
	case 3:// dim3
			switch (cm) {
				case 0: // dim3 type0 (x - y - z)
				for(int i = count; i < 2; ++i){
					switch(i){
						case 0: // w +
						index = (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw + 1], dcg -> dense4[cx + 1][cy][cz][cw + 1], 
										dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1],
										dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1],
										dcg -> dense4[cx][cy + 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
						break;

						case 1: // w -
						index = ((cw - 1) << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz][cw - 1], dcg -> dense4[cx + 1][cy][cz][cw - 1], 
										dcg -> dense4[cx][cy + 1][cz][cw - 1],dcg -> dense4[cx + 1][cy + 1][cz][cw - 1],
										dcg -> dense4[cx][cy][cz + 1][cw - 1],dcg -> dense4[cx + 1][cy][cz + 1][cw - 1],
										dcg -> dense4[cx][cy + 1][cz + 1][cw - 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw - 1]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 4);
						return true;
					}
				}
				return false;

				case 1: // dim3 type1 (x - y - w)
				for(int i = count; i < 2; ++i){
					switch(i){
						case 0: // z +
						index = (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz + 1][cw], dcg -> dense4[cx + 1][cy][cz + 1][cw], 
										dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw],
										dcg -> dense4[cx][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy][cz + 1][cw + 1],
										dcg -> dense4[cx][cy + 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
						break;

						case 1: // z -
						index = (cw << 21) | ((cz - 1) << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy][cz - 1][cw], dcg -> dense4[cx + 1][cy][cz - 1][cw], 
										dcg -> dense4[cx][cy + 1][cz - 1][cw],dcg -> dense4[cx + 1][cy + 1][cz - 1][cw],
										dcg -> dense4[cx][cy][cz - 1][cw + 1],dcg -> dense4[cx + 1][cy][cz - 1][cw + 1],
										dcg -> dense4[cx][cy + 1][cz - 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz - 1][cw + 1]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 4);
						return true;
					}
				}
				return false;

				case 2: // dim3 type2 (x - z - w)
				for(int i = count; i < 2; ++i){
					switch(i){
						case 0: // y +
						index = (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy + 1][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
										dcg -> dense4[cx][cy + 1][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw],
										dcg -> dense4[cx][cy + 1][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1],
										dcg -> dense4[cx][cy + 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
						break;

						case 1: // y -
						index = (cw << 21) | (cz << 14) | ((cy - 1) << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx][cy - 1][cz][cw], dcg -> dense4[cx + 1][cy - 1][cz][cw], 
										dcg -> dense4[cx][cy - 1][cz + 1][cw],dcg -> dense4[cx + 1][cy - 1][cz + 1][cw],
										dcg -> dense4[cx][cy - 1][cz][cw + 1],dcg -> dense4[cx + 1][cy - 1][cz][cw + 1],
										dcg -> dense4[cx][cy - 1][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy - 1][cz + 1][cw + 1]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 4);
						return true;
					}
				}
				return false;

				case 3: // dim3 type3 (y - z - w)
				for(int i = count; i < 2; ++i){
					switch(i){
						case 0: // x +
						index = (cw << 21) | (cz << 14) | (cy << 7) | cx;
						birthday = max({birthtime, dcg -> dense4[cx + 1][cy][cz][cw], dcg -> dense4[cx + 1][cy + 1][cz][cw], 
										dcg -> dense4[cx + 1][cy][cz + 1][cw],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw],
										dcg -> dense4[cx + 1][cy][cz][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz][cw + 1],
										dcg -> dense4[cx + 1][cy][cz + 1][cw + 1],dcg -> dense4[cx + 1][cy + 1][cz + 1][cw + 1]});
						break;

						case 1: // x -
						index = (cw << 21) | (cz << 14) | (cy << 7) | (cx - 1);
						birthday = max({birthtime, dcg -> dense4[cx - 1][cy][cz][cw], dcg -> dense4[cx - 1][cy + 1][cz][cw], 
										dcg -> dense4[cx - 1][cy][cz + 1][cw],dcg -> dense4[cx - 1][cy + 1][cz + 1][cw],
										dcg -> dense4[cx - 1][cy][cz][cw + 1],dcg -> dense4[cx - 1][cy + 1][cz][cw + 1],
										dcg -> dense4[cx - 1][cy][cz + 1][cw + 1],dcg -> dense4[cx - 1][cy + 1][cz + 1][cw + 1]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = BirthdayIndex(birthday, index, 4);
						return true;
					}
				}
				return false;
			}
			return false;
	}
	return false;
}

BirthdayIndex SimplexCoboundaryEnumerator::getNextCoface() {
	return nextCoface;
}
