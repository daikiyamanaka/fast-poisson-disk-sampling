#include "PoissonDiskSampling.h"

PoissonDiskSampling::PoissonDiskSampling(int width, int height, double pitch):width_(width), height_(height), pitch_(pitch)
{
	nodes.resize(width_*height_);

	for(int j=0; j<height_; j++){
		for(int i=0; i<width_; i++){
			nodes[i+j*width] = Node(i, j, 1.0);
		}
	}
	srand(time(0));
	f = Converter();
}

PoissonDiskSampling::~PoissonDiskSampling(){

}

void PoissonDiskSampling::sample(std::vector<Eigen::Vector3d> &points){
	points.resize(0);
	int count = 0;
	std::vector<Node>::iterator n_it;
	std::vector<int> neighbors;
	int max = 100000;
	while(count < max){
		Eigen::Vector3d point = genRandomPoint();
		int pid_x = point[0]/pitch_;
		int pid_y = point[1]/pitch_;
		int n_id = pid_x+pid_y*width_;
		n_it = nodes.begin()+ n_id;
		neighbors = calcNeighborIndex(*n_it, f(n_it->val));

		std::vector<int>::iterator id_it;
		bool free = true;
		for(id_it = neighbors.begin(); id_it != neighbors.end(); id_it++){
			if((nodes.begin()+*id_it)->visited == true){
				free = false;
			}
		}
		if(free){
			points.push_back(point);
			nodes.at(n_id).visited = true;
		}
		count ++;
	}
}

void PoissonDiskSampling::setDensityFunc(const std::vector<double> &dfunc){
	assert((int)dfunc.size() == width_*height_);
	std::vector<double>::const_iterator d_it = dfunc.begin();
	for(std::vector<Node>::iterator it = nodes.begin(); it != nodes.end(); ++it, ++d_it){
		it->val = *d_it;
	}
}

void PoissonDiskSampling::setConverter(boost::function<double(double)> _f){
	f = _f;
}

std::vector<int> PoissonDiskSampling::calcNeighborIndex(Node &node, double r){
	std::vector<int> neighbors;
	double dx = static_cast<double>(node.x_id);
	double dy = static_cast<double>(node.y_id);

	int min_x = (dx*pitch_-r)/pitch_;
	int max_x = (dx*pitch_+r)/pitch_;
	int min_y = (dy*pitch_-r)/pitch_;
	int max_y = (dy*pitch_+r)/pitch_;

	for(int j=min_y; j<=max_y; j++){
		for(int i=min_x; i<=max_x; i++){
			if(i<0 || i>=width_ || j<0 || j>=height_){
				continue;
			}
			neighbors.push_back(i+j*width_);
		}
	}
	return neighbors;
}

Eigen::Vector3d PoissonDiskSampling::genRandomPoint(){
	double x = ((double)rand()+1.0)/((double)RAND_MAX+2.0)*(double)width_;
	double y = ((double)rand()+1.0)/((double)RAND_MAX+2.0)*(double)height_;
	return Eigen::Vector3d(x, y, 0);
}
