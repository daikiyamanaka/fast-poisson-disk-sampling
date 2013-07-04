#include "ShapeOptimizer2D.h"
#include <iostream>

ShapeOptimizer2D::ShapeOptimizer2D(std::vector<Eigen::Vector3d> _points, std::vector<std::vector<int> > _facets, std::vector<std::vector<int> > _neighbors){
	points = _points;
	facets = _facets;
	//neighbors = _neighbors;
	//centroids.resize(facets.size());
	//areas.resize(facets.size());
	//updateTrianglesInfo();
	mesh.createFromFaceVertex(points, facets);
}

ShapeOptimizer2D::~ShapeOptimizer2D(){

}

std::vector<Eigen::Vector3d> ShapeOptimizer2D::getPoints(){
	return mesh.points;
}

std::vector<std::vector<int> > ShapeOptimizer2D::getFacets(){
	return mesh.facets;
}

void ShapeOptimizer2D::optimize(Eigen::Vector3d center){

}

void ShapeOptimizer2D::optimizeStepByStep(Eigen::Vector3d center){
	/*
	for(int i=0; i<facets.size(); i++){

		std::vector<int> neighbor = neighbors[i];
		std::cout << neighbor[0] <<" " << neighbor[1] << " " << neighbor[2] << std::endl;
		double T = 0;
		Eigen::Vector3d c = Eigen::Vector3d::Zero();
		for(int j=0; j<3; j++){
			if(neighbor[j] < 0){
				continue;
			}
			c += centroids[neighbor[j]];
			T += areas[neighbor[j]];
		}
		c /= T;
		points[i] = c;
	}	
	std::cout << "updateTrianglesInfo" << std::endl;
	updateTrianglesInfo();
	*/
	std::vector<float > f_weights(mesh.facets.size());
	std::fill(f_weights.begin(), f_weights.end(), 1.0);

	mesh.CODT(f_weights);
}

void ShapeOptimizer2D::updateTrianglesInfo(){	
	for(int i=0; i<(int)centroids.size(); i++){
		Eigen::Vector3d center = (points[facets[i][0]] + points[facets[i][1]] + points[facets[i][2]])/3;
		centroids[i] = center;
	}

	for(int i=0; i<(int)areas.size(); i++){		
		areas[i] = area(points[facets[i][0]], points[facets[i][1]], points[facets[i][2]]);
	}
}

