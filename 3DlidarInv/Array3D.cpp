#include "Array3D.h"
#include <fstream>

Array3D::Array3D() {
	this->xsize = 0;
	this->ysize = 0;
	this->zsize = 0;
	data = NULL;
}

Array3D::Array3D(int xsize, int ysize, int zsize) {
	this->xsize = xsize;
	this->ysize = ysize;
	this->zsize = zsize;
	data = new double **[xsize];
	for (int x = 0; x < xsize; x++) {
		data[x] = new double*[ysize];
		for (int y = 0; y < ysize; y++) {
			data[x][y] = new double[zsize];
			for (int z = 0; z < zsize; z++) {
				data[x][y][z] = 0;
			}
		}
	}
}

void Array3D::resize(int xsize, int ysize, int zsize) {
	this->xsize = xsize;
	this->ysize = ysize;
	this->zsize = zsize;
	data = new double** [xsize];
	for (int x = 0; x < xsize; x++) {
		data[x] = new double* [ysize];
		for (int y = 0; y < ysize; y++) {
			data[x][y] = new double[zsize];
		}
	}
}

void Array3D::put(int x, int y, int z, double value) {
	data[x][y][z] = value;
}

double Array3D::get(int x, int y, int z) {
	return data[x][y][z];
}

void Array3D::save_as_txt(std::string txt_path, Eigen::Vector3d resolution) {
	std::ofstream outf(txt_path);
	for (int x = 0; x < xsize; x++) {
		for (int y = 0; y < ysize; y++) {
			for (int z = 0; z < zsize; z++) {
				if(data[x][y][z] > 0.001 && data[x][y][z]<9.999)
					outf << (x + 0.5) * (resolution[0]) << " " << (y + 0.5) * (resolution[1]) << " " <<
						(z + 0.5) * (resolution[2]) << " " << data[x][y][z] << std::endl;
			}
		}
	}
	outf.close();
}


Array3D::~Array3D() { 
	//用完数组后，用delete将内存释放
	for (int x = 0; x < xsize; x++) {
		for (int y = 0; y < ysize; y++) {
			delete[] data[x][y];
		}
	}
	for (int x = 0; x < xsize; x++) {
		delete[] data[x];
	}
	delete[] data;
}