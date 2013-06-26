#ifndef POISSONDISKSAMPLING_H_
#define POISSONDISKSAMPLING_H_

#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <Eigen/Sparse>
#include <boost/function.hpp>

class PoissonDiskSampling
{
public:
	class Node{
		public:
		Node(){};
		Node(int x, int y, double df):x_id(x), y_id(y), val(df){visited = false;};
		~Node(){};
		Node(const Node& rn){
			x_id = rn.x_id;
			y_id = rn.y_id;
			val = rn.val;
			visited = rn.visited;
		};
		int x_id, y_id;
		double val;
		bool visited;
	};

	class Converter{
		public:
			Converter(){};
			~Converter(){};		
			double operator()(double val){
				return val;
			};
		};

	PoissonDiskSampling(int width, int height, double pitch);
	~PoissonDiskSampling();

	void sample(std::vector<Eigen::Vector3d> &points);
	void setDensityFunc(const std::vector<double> &density_func);
	void setConverter(boost::function<double(double)> _f);

private:
	std::vector<int> calcNeighborIndex(Node &node, double r);
	Eigen::Vector3d genRandomPoint();		

	boost::function<double(double)> f;
	std::vector<Node> nodes;
	int width_, height_;
	double pitch_;

};
#endif
