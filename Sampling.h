#ifndef SAMPLING_H_
#define SAMPLING_H_

#include <Eigen/Sparse>

class SimpleSampling
{
public:
	SimpleSampling();
	~SimpleSampling();

 void sample(Eigen::VectorXd &weights, std::vector<Eigen::Vector3d> &points, int width, int height, double pitch);


	/* data */
};

#endif
