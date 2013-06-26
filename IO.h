#include <iostream>
#include <vector>
#include <fstream>
#include <Eigen/Core>


void savePly(const std::string &filename, std::vector<Eigen::Vector3d> points){
    //std::cout << "CTMesh.save to " << filename << std::endl;
    std::ofstream ofs;
    ofs.open((filename).c_str());
    // --- header ---
    ofs << "ply" << std::endl;
    ofs << "format ascii 1.0" << std::endl;
    ofs << "comment VCGLIB generated" << std::endl;
    ofs << "element vertex " <<  points.size() << std::endl;
    ofs << "property float x" << std::endl;
    ofs << "property float y" << std::endl;
    ofs << "property float z" << std::endl;
    ofs << "element face " << 0 << std::endl;
    ofs << "property list uchar int vertex_indices" << std::endl;
    ofs << "end_header" << std::endl;

    for(int i=0; i<(int)points.size(); i++){
        ofs << points[i][0] << " " << points[i][1] << " " << points[i][2] << std::endl;
    }
    ofs.close();
}

void savePly(const std::string &filename, std::vector<Eigen::Vector3d> &points, std::vector<std::vector<int> > &faces){
    std::cout << "CTMesh.save to " << filename << std::endl;
    std::ofstream ofs;
    ofs.open((filename).c_str());
    // --- header ---
    ofs << "ply" << std::endl;
    ofs << "format ascii 1.0" << std::endl;
    ofs << "comment VCGLIB generated" << std::endl;
    ofs << "element vertex " <<  points.size() << std::endl;
    ofs << "property float x" << std::endl;
    ofs << "property float y" << std::endl;
    ofs << "property float z" << std::endl;
    ofs << "element face " << faces.size() << std::endl;
    ofs << "property list uchar int vertex_indices" << std::endl;
    ofs << "end_header" << std::endl;

    std::cout << points.size() << std::endl;

    for(int i=0; i<(int)points.size(); i++){
        ofs << points[i][0] << " " << points[i][1] << " " << 0 << std::endl;
    }
    std::cout << faces.size() << std::endl;
    for(int i=0; i<(int)faces.size(); i++){
    	ofs << "3 " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << std::endl;
    }
    ofs.close();
}

void saveVector(std::string name, Eigen::VectorXd &vec, int w, int h){
    // --- output ---
    std::ofstream ofs(name.c_str(), std::ios::out| std::ios::binary);
    float *buff = new float[w];
    for(int y=0; y<h; y++){
        for(int x=0; x<w; x++){
            buff[x] = (float)vec(x+y*w);
        }
        ofs.write((char *)buff, w*sizeof(float));
    }
    delete(buff);
    ofs.close();    
}
void saveVector(std::string name, std::vector<double> &vec, int w, int h){
    // --- output ---
    std::ofstream ofs(name.c_str(), std::ios::out| std::ios::binary);
    float *buff = new float[w];
    for(int y=0; y<h; y++){
        for(int x=0; x<w; x++){
            buff[x] = (float)vec[x+y*w];
        }
        ofs.write((char *)buff, w*sizeof(float));
    }
    delete(buff);
    ofs.close();    
}
