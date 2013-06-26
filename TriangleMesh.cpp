/*
 * TriangleMesh.cpp
 *
 *  Created on: Jun 26, 2012
 *      Author: daikiyamanaka
 */


#include <fstream>
#include <string>
#include <iostream>
#include "TriangleMesh.hpp"

TriangleMesh::TriangleMesh() {
	// TODO Auto-generated constructor stub
    num_of_vertices = 0;
    num_of_faces = 0;
    ave_length = -1;
}


TriangleMesh::~TriangleMesh() {
    // TODO Auto-generated destructor stub

}

void TriangleMesh::clear(){
    points.clear();
    centroids.clear();
    facets.clear();
    neighbors.clear();
    v_f_map.clear();
    areas.clear();
    scalars.clear();
    point_scalars.clear();
}

void TriangleMesh::readCTM(std::string filename){

    std::cout << filename.c_str() << std::endl;

    FILE* in = fopen(filename.c_str(), "r");
    fscanf(in, "%d", &num_of_vertices);
    //ver = new float[verN][2];
    points.resize(num_of_vertices);

    fscanf(in, "%d", &num_of_faces);
    facets.resize(num_of_faces);
    scalars.resize(num_of_faces);

    std::cout << "num_of_vertices : " << num_of_vertices << std::endl;
    std::cout << "num_of_faces : " << num_of_faces << std::endl;

    for(int i=0; i<num_of_vertices; i++){
        float x, y;
      fscanf(in, "%g %g", &x, &y);
      points[i][0] = x;
      points[i][1] = y;
      points[i][2] = 0;
    }

    for(int i=0; i<num_of_faces; i++){
        int p0, p1, p2;
        fscanf(in, "%d %d %d", &p0, &p1, &p2);
        std::vector<int> face(3);
        face[0] = p0;
        face[1] = p1;
        face[2] = p2;
        facets[i] = face;
        //std::cout << face << std::endl;
    }

    for(int i=0; i<num_of_faces; i++){
        float value;
        fscanf(in, "%g", &value);
        scalars[i] = value;
    }

    fclose(in);

    normalizeScalar();
    normalize();
    initialize();

    //std::cout << "ctm file readed" << std::endl;

}

void TriangleMesh::createFromFaceVertex(std::vector<Eigen::Vector3d> &vertices, std::vector<std::vector<int> > &faces){
    clear();

    num_of_vertices = vertices.size();
    num_of_faces = faces.size();

    points.resize(num_of_vertices);
    facets.resize(num_of_faces);

    for(int i=0; i<num_of_vertices; i++){
        points[i][0] = vertices[i][0];
        points[i][1] = vertices[i][1];
        points[i][2] = vertices[i][2];
    }
    for(int i=0; i<num_of_faces; i++){
        std::vector<int> face(3);
        face[0] = faces[i][0];
        face[1] = faces[i][1];
        face[2] = faces[i][2];
        facets[i] = face;
        //std::cout << face[0] << " " << face[1] << " " << face[2] << std::endl;
        //facets.push_back(face);
    }
    scalars.resize(num_of_faces);
    std::fill(scalars.begin(), scalars.end(), 0);
    initialize();

    std::cout << "initial mesh crated from vertex and face list" <<  std::endl;
}

void TriangleMesh::updateGeometry(){

    for(int i=0; i<num_of_faces; i++){
        std::vector<int> face = facets[i];
        centroids[i] = centroid(points[face[0]], points[face[1]], points[face[2]]);
        circumcenters[i] = circumcenter(points[face[0]], points[face[1]], points[face[2]]);
        areas[i] = AREA2D(points[face[0]], points[face[1]], points[face[2]]);
    }

}

void TriangleMesh::add(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c){

	std::vector<int> face(3);

	if(std::find(points.begin(), points.end(), a) == points.end()){
		points.push_back(a);
	}
	if(std::find(points.begin(), points.end(), b) == points.end()){
		points.push_back(b);
	}
	if(std::find(points.begin(), points.end(), c) == points.end()){
		points.push_back(c);
	}

	face[0] = std::find(points.begin(), points.end(), a) - points.begin();
	face[1] = std::find(points.begin(), points.end(), b) - points.begin();
	face[2] = std::find(points.begin(), points.end(), c) - points.begin();

	facets.push_back(face);

}

void TriangleMesh::ODT(std::vector<float> &p_weights){
    for(int i=0; i<num_of_vertices; i++){
        if(boundary_marker[i]){
            continue;
        }
        std::vector<int> ring_faces = v_f_map[i];
        float weight = 0;
        Eigen::Vector3d weighted_p = Eigen::Vector3d::Zero();
        for(int j=0; j<ring_faces.size(); j++){
            int n_index = ring_faces[j];
            weight += p_weights[n_index]*areas[n_index];
            weighted_p += p_weights[n_index]*areas[n_index]*circumcenters[n_index];
        }
        if(ring_faces.size() != 0){
            points[i] = weighted_p/weight;
        }
    }
    updateGeometry();
}

void TriangleMesh::CODT(std::vector<float> &f_weights){
    for(int i=0; i<num_of_vertices; i++){
        if(boundary_marker[i]){
            continue;
        }
        std::vector<int> ring_faces = v_f_map[i];
        float weight = 0;
        Eigen::Vector3d weighted_p = Eigen::Vector3d::Zero();
        for(int j=0; j<ring_faces.size(); j++){
            int n_index = ring_faces[j];
            weight += f_weights[n_index]*areas[n_index];
            weighted_p += f_weights[n_index]*areas[n_index]*centroids[n_index];
        }
        if(ring_faces.size() != 0){
            points[i] = weighted_p/weight;
        }
    }
    updateGeometry();
}

void TriangleMesh::getCentroids(std::vector<Eigen::Vector3d> &_centroids){
    _centroids = centroids;
}

void TriangleMesh::getVertex(int index, Eigen::Vector3d &v){
    v = points[index];
}

void TriangleMesh::setVertex(int index, Eigen::Vector3d &v){
    points[index] = v;
}

void TriangleMesh::getFace(int index, std::vector<int> &f){
    f = facets[index];
}

void TriangleMesh::setEdges(std::vector<std::pair<int, int> > &_edges){
    edges = _edges;
    computeAverageLength();
}

void TriangleMesh::getEdges(std::vector<std::pair<int, int> > &_edges){
    _edges = edges;
}

float TriangleMesh::getAverageLength(){
    return ave_length;
}

std::vector<int> TriangleMesh::getNeighbors(int index){
    return neighbors[index];
}

int TriangleMesh::get_number_of_faces(){
//	return facets.size();
    return num_of_faces;
}

int TriangleMesh::get_number_of_vertices(){
    return num_of_vertices;
}

void TriangleMesh::normalizeScalar(){
    float max, min;
    max = min = scalars[0];
    for(int i=0; i<num_of_faces; i++){
        if(scalars[i] > max){
            max = scalars[i];
        }
        else if(scalars[i] < min){
            min = scalars[i];
        }
    }
    for(int i=0; i<num_of_faces; i++){
        scalars[i] = (scalars[i] - min)/(max - min);
    }
}

void TriangleMesh::normalize(){
    float max[3], min[3], mid[3];
    for(int i=0; i<3; i++){
        max[i] = min[i] = points[0][i];
        for(int j=1; j<num_of_vertices; j++){
            float v = points[j][i];
            if(max[i] < v){
                max[i] = v;
            }
            else if(min[i] > v){
                min[i] = v;
            }
        }
        mid[i] = 0.5*(max[i] + min[i]);
    }
    float size = max[0] - min[0];

    for(int i=0; i<num_of_vertices; i++){
        points[i][0] = (points[i][0] - mid[0])/size;
        points[i][1] = (points[i][1] - mid[1])/size;
        points[i][2] = (points[i][2] - mid[2])/size;
    }
}

void TriangleMesh::initialize(){

    std::cout << "initialize" << std::endl;

    int *ringN = new int[num_of_vertices];
    for(int i=0; i<num_of_vertices; i++)
        ringN[i]= 0;
    for(int i=0; i<num_of_faces; i++)
        for(int j=0; j<3; j++)
            ringN[facets[i][j]]++;

    std::cout << "initialize 1" << std::endl;

    int **ring = new int*[num_of_vertices];
    for(int i=0; i<num_of_vertices; i++){
       ring[i] = new int[ringN[i]];
       ringN[i] = 0;
    }

    for(int i=0; i<num_of_faces; i++)
        for(int j=0; j<3; j++)
            ring[facets[i][j]][ringN[facets[i][j]]++] = i;

    std::cout << "initialize 2" << std::endl;
    //adT = new int[triN][3];
    neighbors.resize(num_of_faces);
    for(int i=0; i<num_of_faces; i++){
        std::vector<int> face;
        face.resize(3);
        neighbors[i] = face;
    }

    for(int i=0; i<num_of_faces; i++)
        neighbors[i][0] = neighbors[i][1] = neighbors[i][2] = -1;
    for(int i=0; i<num_of_vertices; i++){
        int rN = ringN[i];
        int* r = ring[i];
        for(int j=0; j<rN; j++){
            std::vector<int> t = facets[r[j]];
            int m = 2;
            if(t[1] == i)
                m = 0;
            else if(t[2] == i)
                m = 1;
            if(neighbors[r[j]][m] != -1)
                continue;
            int v = t[(m+2)%3];
            for(int k=0; k<rN; k++){
                if(k == j)
                    continue;
            std::vector<int> tk = facets[r[k]];
            for(int l=0; l<3; l++)
                if(tk[l] == v && tk[(l+1)%3] == i){
                    neighbors[r[j]][m] = r[k];
                    neighbors[r[k]][(l+2)%3] = r[j];
                    break;
                }
            }
        }
    }

    boundary_marker.resize(num_of_vertices);
    std::fill(boundary_marker.begin(), boundary_marker.end(), false);

    // boundary
    for(int i=0; i<num_of_faces; i++){
        for(int j=0; j<3; j++){
            if(neighbors[i][j] == -1){
                boundary_marker[facets[i][(j+1)%3]] = true;
                boundary_marker[facets[i][(j+2)%3]] = true;
            }
        }
    }


    //circumcenter
    circumcenters.resize(num_of_faces);
    for(int i=0; i<num_of_faces; i++){
        circumcenters[i] = circumcenter(points[facets[i][0]], points[facets[i][1]], points[facets[i][2]]);
    }


    // area
    areas.resize(num_of_faces);
    // vertex face map
    v_f_map.resize(num_of_vertices);

    for(int i=0; i<num_of_faces; i++){
        std::vector<int> face = facets[i];
        areas[i] = AREA2D(points[face[0]], points[face[1]], points[face[2]]);
        for(int j=0; j<3; j++){
            if(face[j] != -1){
                v_f_map[face[j]].push_back(i);
            }
        }
    }

    /*
    point_scalars.resize(num_of_vertices);
    // --- compute point scalars
    for(int i=0; i<num_of_vertices; i++){
        float sum = 0;
        float weight_sum = 0;
        for(int j=0; j<ringN[i]; j++){
            int n_tri = ring[i][j];
            std::vector<int> face = facets[n_tri];
            Eigen::Vector3f c = centroid(points[face[0]], points[face[1]], points[face[2]]);
            float weight = (points[i] - c).norm();
            weight = 1/weight;
            sum += weight*scalars[n_tri];
            weight_sum += weight;
        }
        point_scalars[i] = sum/weight_sum;
    }
    */

    centroids.resize(num_of_faces);
    for(int i=0; i<num_of_faces; i++){
        centroids[i] = centroid(points[facets[i][0]],
                                points[facets[i][1]],
                                points[facets[i][2]]);
    }

    for(int i=0; i<num_of_vertices; i++)
      delete[] ring[i];
    delete[] ring;
    delete[] ringN;
}

Eigen::Vector3d TriangleMesh::centroid(Eigen::Vector3d &a, Eigen::Vector3d &b, Eigen::Vector3d &c){
    return (a+b+c)/3.0;
}

Eigen::Vector3d TriangleMesh::circumcenter(Eigen::Vector3d &a, Eigen::Vector3d &b, Eigen::Vector3d &c){
    double xba, yba, xca, yca;
    double balength, calength;
    double denominator;
    double xcirca, ycirca;

    Eigen::Vector3d center;

    /* Use coordinates relative to point `a' of the triangle. */
    xba = b[0] - a[0];
    yba = b[1] - a[1];
    xca = c[0] - a[0];
    yca = c[1] - a[1];
    /* Squares of lengths of the edges incident to `a'. */
    balength = xba * xba + yba * yba;
    calength = xca * xca + yca * yca;

    /* Calculate the denominator of the formulae. */
  //#ifdef EXACT
    /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
    /*   to ensure a correctly signed (and reasonably accurate) result, */
    /*   avoiding any possibility of division by zero.                  */
    //denominator = 0.5 / orient2d(b, c, a);
  //#else
    /* Take your chances with floating-point roundoff. */
    denominator = 0.5 / (xba * yca - yba * xca);
  //#endif

    /* Calculate offset (from `a') of circumcenter. */
    xcirca = (yca * balength - yba * calength) * denominator;
    ycirca = (xba * calength - xca * balength) * denominator;
    center[0] = xcirca + a[0];
    center[1] = ycirca + a[1];
    center[2] = 0;

    return center;

}

void TriangleMesh::computeAverageLength(){

    if(ave_length > 0){
        return;
    }

    float length = 0;
    for(int i=0; i<edges.size(); i++){
        int p0 = edges[i].first;
        int p1 = edges[i].second;
        length += sqrt(pow(points[p0][0] - points[p1][0], 2) + pow(points[p0][1] - points[p1][1], 2));
    }
    ave_length = length/((float) edges.size());
}
