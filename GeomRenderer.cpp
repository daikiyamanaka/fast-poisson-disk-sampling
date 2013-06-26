#include "GeomRenderer.h"

GeomRenderer::GeomRenderer(){
	init();
}
GeomRenderer::GeomRenderer(std::vector<Eigen::Vector3d> &_points){
	points = _points;
	init();
}
GeomRenderer::GeomRenderer(std::vector<Eigen::Vector3d> &_points, std::vector<std::vector<int> > &_facets){
	points = _points;
	facets = _facets;
	init();
}
GeomRenderer::~GeomRenderer(){
  cairo_destroy (cr);
  cairo_surface_destroy (surface);
}

void GeomRenderer::render(int width, int height){
	if(points.size() == 0)
		return;
	min = max = points[0];
	for(int i=0; i<(int)points.size(); i++){
		Eigen::Vector3d point = points[i];
		for(int j=0; j<3; j++){
			if(point[j] > max[j]){
				max[j] = point[j];
			}
			else if(point[j] < min[j]){
				min[j] = point[j];
			}
		}
	}
	ratio = 0.1;

	origin_x = width*ratio/2.0;
	origin_y = height*ratio/2.0;

	scale_x = ((double)width*(1-ratio))/(max[0] - min[0]);
	scale_y = ((double)height*(1-ratio))/(max[1] - min[1]);		

	surface = 
	  cairo_image_surface_create (CAIRO_FORMAT_ARGB32, width, height);
  cr = cairo_create (surface);

  cairo_set_source_rgba (cr, 0.1, 0.1, 0.1, 0.6);
  cairo_set_line_width(cr, line_size);
  for(int i=0; i<facets.size(); i++){
  	for(int j=0; j<facets[i].size(); j++){
  		Eigen::Vector3d p_point, n_point;
  		p_point = points[facets[i][j]];
  		if(j == facets[i].size()-1){
  			n_point = points[facets[i][0]];
  		}
  		else{
  			n_point = points[facets[i][j+1]];
  		}
  		cairo_move_to(cr, (p_point[0]- min[0])*scale_x + origin_x, (p_point[1] - min[1])*scale_y + origin_y);
  		cairo_line_to(cr, (n_point[0]- min[0])*scale_x + origin_x, (n_point[1] - min[1])*scale_y + origin_y);  		
  		cairo_stroke(cr);
  	}
  }

  cairo_set_source_rgba (cr, 1, 0.2, 0.2, 1.0);
  //cairo_set_line_width (cr, 6.0);

  for(int i=0; i<(int)points.size(); i++){  	
  	double x = (points[i][0] - min[0])*scale_x + origin_x;
  	double y = (points[i][1] - min[1])*scale_y + origin_y;  	

  	cairo_arc (cr, x, y, point_size, 0, 2*M_PI);
  	cairo_fill (cr);
  }
}

void GeomRenderer::renderPoint(Eigen::Vector3d point, Eigen::Vector3d color, double alpha){
	  cairo_set_source_rgba (cr, color[0], color[1], color[2], alpha);
  	double x = (point[0] - min[0])*scale_x + origin_x;
  	double y = (point[1] - min[1])*scale_y + origin_y;  	
  	cairo_arc (cr, x, y, point_size, 0, 2*M_PI);
  	cairo_fill (cr);	
}

void GeomRenderer::renderEdge(Eigen::Vector3d p, Eigen::Vector3d n, Eigen::Vector3d color, double alpha){
	  cairo_set_source_rgba (cr, color[0], color[1], color[2], alpha);	  
	  cairo_set_line_width(cr, line_size);
		cairo_move_to(cr, (p[0]- min[0])*scale_x + origin_x, (p[1] - min[1])*scale_y + origin_y);
		cairo_line_to(cr, (n[0]- min[0])*scale_x + origin_x, (n[1] - min[1])*scale_y + origin_y);  		
		cairo_stroke(cr);
}

void GeomRenderer::save(std::string filename){
	  cairo_surface_write_to_png (surface, filename.c_str());
}

void GeomRenderer::init(){
	surface = NULL;
	cr = NULL;
	point_size = 3.0;
	line_size = 3.0;
}

/*
void GeomRenderer::save(std::string filename, int width, int height){
	if(points.size() == 0)
		return;
	Eigen::Vector3d min, max;
	min = max = points[0];
	for(int i=0; i<(int)points.size(); i++){
		Eigen::Vector3d point = points[i];
		for(int j=0; j<3; j++){
			if(point[j] > max[j]){
				max[j] = point[j];
			}
			else if(point[j] < min[j]){
				min[j] = point[j];
			}
		}
	}


	double ratio = 0.1;

	double origin_x = width*ratio/2.0;
	double origin_y = height*ratio/2.0;

	double scale_x = ((double)width*(1-ratio))/(max[0] - min[0]);
	double scale_y = ((double)height*(1-ratio))/(max[1] - min[1]);	

	cairo_surface_t *surface = 
	  cairo_image_surface_create (CAIRO_FORMAT_ARGB32, width, height);
  cairo_t *cr = cairo_create (surface);

  cairo_set_source_rgba (cr, 0.1, 0.1, 0.1, 0.6);
  cairo_set_line_width(cr, 3.0);
  for(int i=0; i<facets.size(); i++){
  	for(int j=0; j<facets[i].size(); j++){
  		Eigen::Vector3d p_point, n_point;
  		p_point = points[facets[i][j]];
  		if(j == facets[i].size()-1){
  			n_point = points[facets[i][0]];
  		}
  		else{
  			n_point = points[facets[i][j+1]];
  		}
  		cairo_move_to(cr, (p_point[0]- min[0])*scale_x + origin_x, (p_point[1] - min[1])*scale_y + origin_y);
  		cairo_line_to(cr, (n_point[0]- min[0])*scale_x + origin_x, (n_point[1] - min[1])*scale_y + origin_y);  		
  		cairo_stroke(cr);
  	}
  }

  cairo_set_source_rgba (cr, 1, 0.2, 0.2, 1.0);
  cairo_set_line_width (cr, 6.0);

  for(int i=0; i<(int)points.size(); i++){  	
  	double x = (points[i][0] - min[0])*scale_x + origin_x;
  	double y = (points[i][1] - min[1])*scale_y + origin_y;  	

  	cairo_arc (cr, x, y, 3.0, 0, 2*M_PI);
  	cairo_fill (cr);
  }

  cairo_destroy (cr);
  cairo_surface_write_to_png (surface, filename.c_str());
  cairo_surface_destroy (surface);
}
*/