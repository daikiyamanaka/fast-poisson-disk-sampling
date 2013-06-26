/*
 * TriangleMesh.h
 *
 *  Created on: Jun 26, 2012
 *      Author: daikiyamanaka
 */

#ifndef TRIANGLEMESH_H_
#define TRIANGLEMESH_H_

#include <vector>
#include <Eigen/Dense>


/*
static Eigen::Vector3f centroid(Eigen::Vector3f &a, Eigen::Vector3f &b, Eigen::Vector3f &c){
    return (a+b+c)/3;
}
*/
class TriangleMesh {



public:
	TriangleMesh();
	virtual ~TriangleMesh();

    void clear();
	void read(std::string file_name);
	void write(std::string file_name);
    void readCTM(std::string filename);
    void createFromFaceVertex(std::vector<Eigen::Vector3d> &vertices, std::vector<std::vector<int> > &faces);
    void updateGeometry();

    // --- for Editing
	void add(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c);
    void ODT(std::vector<float> &p_weights);
    void CODT(std::vector<float> &p_weights);    

    // --- for IO
    void getCentroids(std::vector<Eigen::Vector3d> &_centroids);
    void getVertex(int index, Eigen::Vector3d &v);
    void setVertex(int index, Eigen::Vector3d &v);
    void getFace(int index, std::vector<int> &f);
    void setEdges(std::vector<std::pair<int, int> > &edges);
    void getEdges(std::vector<std::pair<int, int> > &edges);
    float getAverageLength();
    std::vector<int> getNeighbors(int index);
    int get_number_of_faces();
    int get_number_of_vertices();

    inline static double AREA2D(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c){
        return 0.5*((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]));
    }


    static Eigen::Vector3d centroid(Eigen::Vector3d &a, Eigen::Vector3d &b, Eigen::Vector3d &c);
    static Eigen::Vector3d circumcenter(Eigen::Vector3d &a, Eigen::Vector3d &b, Eigen::Vector3d &c);

    // point
	std::vector<Eigen::Vector3d> points;
    std::vector<std::vector<int> > v_f_map;
    std::vector<float> point_scalars;
    std::vector<bool> boundary_marker;

    // edge
    std::vector<std::pair<int, int> > edges;
    // face
    std::vector<Eigen::Vector3d> centroids;
    std::vector<Eigen::Vector3d> circumcenters;
	std::vector<std::vector<int> > facets;
    std::vector<std::vector<int> > neighbors;
    std::vector<float> areas;
    std::vector<float> scalars;

private:
    void normalizeScalar();
    void normalize();
    void initialize();

    void computeAverageLength();

    int num_of_vertices;
    int num_of_faces;

    float ave_length;

};

#endif /* TRIANGLEMESH_H_ */
