//DenseCubicalGrids.h
#include <string>

enum file_format { DIPHA, PERSEUS };

class DenseCubicalGrids { // file_read

public:
	double threshold;
	int dim;
	int ax, ay, az, aw;
	double dense4[128][128][128][128];
	file_format format;

	DenseCubicalGrids(const std::string& filename, double _threshold, file_format _format) ; 
	double getBirthday(int index, int dim);
};
