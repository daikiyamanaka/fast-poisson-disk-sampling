#include <vector>
#include <Eigen/Core>
#include "TriangleMesh.hpp"

class ShapeOptimizer2D
{
public:
	ShapeOptimizer2D(std::vector<Eigen::Vector3d> points, std::vector<std::vector<int> > facets, std::vector<std::vector<int> > _neighbors);
	~ShapeOptimizer2D();

	std::vector<Eigen::Vector3d> getPoints();
	std::vector<std::vector<int> > getFacets();

	void optimize(Eigen::Vector3d center);
	void optimizeStepByStep(Eigen::Vector3d center);

private:
	void updateTrianglesInfo();
	inline static double area(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c){
		return 0.5*((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]));
  }
	std::vector<Eigen::Vector3d> points;
	std::vector<std::vector<int> > facets;
	std::vector<std::vector<int> > neighbors;	
	std::vector<Eigen::Vector3d> centroids;
	std::vector<double> areas;
	TriangleMesh mesh;
};
