#ifndef ARRAY3D_H_
#define ARRAY3D_H_
#include <Eigen/Core>
/*
*
* x从左到右
* y从下往上
* z高度
*/

class Array3D {
public:
	int xsize, ysize, zsize;
	double ***data;
public:
	Array3D();
	Array3D(int xsize, int yszie, int zsize);

	void resize(int xsize, int yszie, int zsize);

	//Put data into array
	void put(int x, int y, int z, double value);
	double get(int x, int y, int z);

	void save_as_txt(std::string txt_path, Eigen::Vector3d resolution);

	double operator[](Eigen::Vector3i coord) const { return data[coord.x()][coord.y()][coord.z()]; }
	double& operator[](Eigen::Vector3i coord) { return data[coord.x()][coord.y()][coord.z()]; }

	~Array3D();

};


#endif // !ARRAY3D_H_
