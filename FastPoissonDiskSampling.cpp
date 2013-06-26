#include <iostream>
#include <string>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "Sampling.h"
#include "PoissonDiskSampling.h"
#include "triangle.h"
#include "GeomRenderer.h"
#include "ShapeOptimizer2D.h"
#include "IO.h"

static const int width = 128;
static const int height = 128;

// l = √3d((1+√(1-σ))/σ)

class MyConverter : public PoissonDiskSampling::Converter{
	public:	
		double operator()(double val){
			if(fabs(val) < thresh()){
				return max();
			}
			return sqrt(3)*0.1*(1+sqrt(1-val))/val;
		};
		static double thresh(){
			return 0.01;
		}
		static double max(){
			return 100;
		}
	};

int main(int argc, char ** argv){
	std::string outname;
	if(argc == 2){
		outname = std::string(argv[1]);
	}	
	else{
		outname = "out.ply";
	}	

	Eigen::VectorXd weight = Eigen::VectorXd(width*height);
	std::vector<double> w_vec(width*height);

	for(int y=0; y<height; y++){
		for(int x=0; x<width; x++){
			weight(x+y*width) = (double)x*(double)y/((double)width*(double)height);
			w_vec[x+y*width] = weight(x+y*width);
		}
	}

	//saveVector("density.raw", w_vec, width, height);

	std::vector<double > l_vec = w_vec;
	std::transform(l_vec.begin(), l_vec.end(), l_vec.begin(), MyConverter());

	//saveVector("radius.raw", l_vec, width, height);

	std::vector<Eigen::Vector3d> points;
	PoissonDiskSampling pds(width, height, 1.0);
	std::cout << "setDensityFunc" << std::endl;
	pds.setDensityFunc(w_vec);
	std::cout << "setConverter" << std::endl;		
	MyConverter converter;
	pds.setConverter(converter);
	std::cout << "sample" << std::endl;	
	pds.sample(points);
	savePly("pds.ply", points);

	GeomRenderer gr(points);
	gr.render(512, 512);
	gr.save("pds.png");

	/*
	std::vector<Eigen::Vector3d> points;
	SimpleSampling ss;
	ss.sample(weight, points, width, height, 1.0);

	triangulateio tri, out, vorout;
	out.pointlist = NULL;
	out.trianglelist = NULL;
	out.neighborlist = NULL;
	tri.numberofpoints = points.size();
	tri.pointlist = new REAL[tri.numberofpoints*2];
	for(int i=0; i<(int)points.size(); i++){
		tri.pointlist[i*2] = points[i][0];
		tri.pointlist[i*2+1] = points[i][1];
	}
	triangulate("zn", &tri, &out, &vorout);
	std::cout << "triangulated " << std::endl;
	std::vector<std::vector<int> > faces;//(out.numberoftriangles);
	std::vector<std::vector<int> > neighbors;
	points.resize(out.numberofpoints);
	for(int i=0; i<(int)points.size(); i++){
		points[i][0] = out.pointlist[i*2];
		points[i][1] = out.pointlist[i*2+1];	
		points[i][2] = 0;

	}	
	for(int i=0; i<out.numberoftriangles; i++){
		std::vector<int> face(3);
		face[0] = out.trianglelist[i*3];
		face[1] = out.trianglelist[i*3+1];
		face[2] = out.trianglelist[i*3+2];				
		//faces[i] = face;
		faces.push_back(face);
	}
	neighbors.resize(out.numberoftriangles);
	for(int i=0; i<(int)neighbors.size(); i++){
		std::vector<int> neighbor(3);
		neighbor[0] = out.neighborlist[i*3];		
		neighbor[1] = out.neighborlist[i*3+1];		
		neighbor[2] = out.neighborlist[i*3+2];						
		neighbors[i] = neighbor;
	}

	ShapeOptimizer2D so(points, faces, neighbors);
	Eigen::Vector3d ideal_center = Eigen::Vector3d((double)width/1.5, (double)height/1.5, 0.0);

	for(int i=0; i<100; i++){
			std::cout << "optimizing ... step " << i+1 << std::endl;
		so.optimizeStepByStep(ideal_center);
	}
	points = so.getPoints();
	faces = so.getFacets();
	GeomRenderer gr(points, faces);
	gr.render(640, 480);

	Eigen::Vector3d c(0.2, 0.2, 1.0);	
	gr.renderPoint(ideal_center, c, 1.0);
	gr.save("odt.png");
	*/
}
