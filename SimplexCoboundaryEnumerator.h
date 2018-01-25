//SimplexCoboundaryEnumerator.h

class SimplexCoboundaryEnumerator
{
public:
	BirthdayIndex simplex;
	DenseCubicalGrids* dcg;
	int dim;
	double birthtime;
	int ax, ay, az, aw;
	int cx, cy, cz, cw, cm;
	int count;
	BirthdayIndex nextCoface;
	double threshold;
	const int mode = 1; // 0 -> eight neighbourhoods 、 1 -> four neighbourhoods 

	SimplexCoboundaryEnumerator();
	void setSimplexCoboundaryEnumerator(BirthdayIndex _s, DenseCubicalGrids* _dcg);
	bool hasNextCoface() ;
	BirthdayIndex getNextCoface();
};
