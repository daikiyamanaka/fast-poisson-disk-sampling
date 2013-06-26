#include <Eigen/Core>
#include <cairo/cairo.h>
#include <vector>
#include <string>

class GeomRenderer
{
public:
	GeomRenderer();
	GeomRenderer(std::vector<Eigen::Vector3d> &points);
	GeomRenderer(std::vector<Eigen::Vector3d> &points, std::vector<std::vector<int> > &facets);	
	~GeomRenderer();

	void render(int width, int height);
	void renderPoint(Eigen::Vector3d p, Eigen::Vector3d color, double alpha);
	void renderEdge(Eigen::Vector3d p, Eigen::Vector3d n, Eigen::Vector3d color, double alpha);
	void save(std::string filename);
	//void save(std::string filename, int width, int height);

private:
	void init();
	std::vector<Eigen::Vector3d> points;
	std::vector<std::vector<int> > facets;

	// for geometry
	Eigen::Vector3d min, max;
	double ratio, origin_x, origin_y, scale_x, scale_y;
	cairo_surface_t *surface; 
	cairo_t *cr;

	// for rendering
	double point_size,line_size;
};
