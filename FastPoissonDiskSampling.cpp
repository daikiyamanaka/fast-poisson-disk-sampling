#include <iostream>
#include <string>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "PoissonDiskSampling.h"
#include "GeomRenderer.h"
#include "IO.h"

static const int width = 512;
static const int height = 512;

// functor to detemine distance between points
class UniformConverter : public PoissonDiskSampling::Converter{
public:
    UniformConverter(double len){
        len_ = len;
    }
	double operator()(double val){
        return len_;
	};
private:
    double len_;
};

int main(int argc, char ** argv){
	std::string outname;
	if(argc == 2){
		outname = std::string(argv[1]);
	}
	else{
		outname = "out.ply";
	}

	std::vector<Eigen::Vector3d> points;
	PoissonDiskSampling pds(width, height, 1.0);
    UniformConverter converter(10.0); // distance between points should 10.0
	pds.setConverter(converter);
	pds.sample(points);
	savePly("pds.ply", points);

	GeomRenderer gr(points);
	gr.render(width, height);
	gr.save("pds.png");
}
