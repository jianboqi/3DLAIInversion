#ifndef _GRIDESTIMATOR_H_
#define _GRIDESTIMATOR_H_
#include "PADEstimator.h"
#include <Eigen/Core>
#include "Array3D.h"

#define RAY_O_OFFSET 1000

class GridPADEstimator:public PADEstimator{
public:
	Eigen::Vector3d resolution;

	Array3D pad_3d; //PAD
	Array3D inc_3d; //��������
	Array3D out_3d; //��������

	Array3D num_rays_3d; //�����Ĺ�������
	Array3D path_lengths_3d; //������·�������ۼ�

	Eigen::Vector3d min_xyz;
	Eigen::Vector3d max_xyz;
	Eigen::Vector3i matrix_dims; //���г߶�
	//int pad_xsize, pad_ysize, pad_zsize;

public:
	GridPADEstimator(double voxel_size_x, double voxel_size_y, double voxel_size_z) {
		this->resolution[0] = voxel_size_x;
		this->resolution[1] = voxel_size_y;
		this->resolution[2] = voxel_size_z;
	}

	GridPADEstimator(double size) {
		this->resolution[0] = size;
		this->resolution[1] = size;
		this->resolution[2] = size;
	}

	void set_resolution(double r) {
		this->resolution[0] = r;
		this->resolution[1] = r;
		this->resolution[2] = r;
	}

	//��ʼ��
	void init(std::vector<std::string> & input_files);

	

	void virtual pad_inverse(std::vector<std::string>& input_files);
	void virtual save_to_file(std::string out_file_path, OutFileFormat out_format = OutFileFormat::TXT);

	/// <summary>
	/// ����LAI
	/// </summary>
	void lai_inverse();

private:
	//ִ�й���׷��
	void ray_traversal_per_flight();

	//������voxel���㷨
	std::vector<Eigen::Vector3i> voxel_traversal(Eigen::Vector3d ray_start, Eigen::Vector3d ray_end, std::vector<double>& traversed_lengths);
};



#endif // !_GRIDESTIMATOR_H_

