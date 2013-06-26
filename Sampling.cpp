#include "Sampling.h"
#include <ctime>
#include <iostream>
SimpleSampling::SimpleSampling(){
	srand(time(0));
}
SimpleSampling::~SimpleSampling(){
}
void SimpleSampling::sample(Eigen::VectorXd &weights, 
	                          std::vector<Eigen::Vector3d> &points, 
	                          int width, int height, double pitch)
{
 	//float uniform = ((double)rand()+1.0)/((double)RAND_MAX+2.0);
	points.resize(0);
 	for(int j=0; j<height; j++){
 		for(int i=0; i<width; i++){ 		
 		 	float uniform = ((double)rand()+1.0)/((double)RAND_MAX+2.0);
 		 	std::cout << uniform << " " << weights(i+j*width) << std::endl;
 		 	if(uniform < weights(i+j*width)){
 		 		points.push_back(Eigen::Vector3d(i*pitch, j*pitch, 0)); 		 		
 		 	}
 		}
	}
}
