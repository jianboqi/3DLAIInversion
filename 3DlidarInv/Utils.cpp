#include "Utils.h"
#include <fstream>
#include <iomanip>


void output_obj(delaunator::Delaunator d, std::string file_path) {
	std::ofstream out;
	out.open(file_path, std::ios::out);
	for (std::size_t i = 0; i < d.coords.size()/2; i++) {
		out << "v " << std::setprecision(15)<< d.coords[2 * i] << " " << d.coords[2 * i + 1] << " 0" << std::endl;
	}
	/*for (std::size_t i = 0; i < d.triangles.size(); i += 3) {
		out << "v " << std::setprecision(15) << std::setprecision(15) << d.coords[2 * d.triangles[i]] << " " << d.coords[2 * d.triangles[i] + 1] << " 0" << std::endl;
		out << "v " << std::setprecision(15) << std::setprecision(15) << d.coords[2 * d.triangles[i+1]] << " " << d.coords[2 * d.triangles[i+1] + 1] << " 0" << std::endl;
		out << "v " << std::setprecision(15) << std::setprecision(15) << d.coords[2 * d.triangles[i+2]] << " " << d.coords[2 * d.triangles[i+2] + 1] << " 0" << std::endl;
	}*/
	for (std::size_t i = 0; i < d.triangles.size(); i += 3) {
		out << "f " << d.triangles[i] + 1 << " " << d.triangles[i + 1] + 1 << " " << d.triangles[i + 2] + 1 << std::endl;
	}
	out.close();
}