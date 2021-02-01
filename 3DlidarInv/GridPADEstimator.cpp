#include "GridPADEstimator.h"
#include <lasreader.hpp>
#include <cmath>


void GridPADEstimator::pad_inverse(std::vector<std::string>& input_files) {
	std::cout << " - Estimator: GridPADEstimator" << std::endl;
    this->init(input_files);
	//计算每个输入航线文件中的pulse相关信息，比如类型、方向、最邻近纯地面脉冲
	std::cout << " - Retrieving pulses information..." << std::endl;
	for (int i = 0; i < input_files.size(); i++) {
		std::cout << " - Processing file: " << input_files[i] << std::endl;
		this->pulse_info_per_flight(input_files[i]);
		//执行光线追踪
		std::cout << " - Do ray traversal... "<< std::endl;
		this->ray_traversal_per_flight();
	}
	std::cout << " - LAI inversion... " << std::endl;
	this->lai_inverse();
	std::cout << " - Done." << std::endl;
}

void GridPADEstimator::save_to_file(std::string out_file_path, OutFileFormat out_format) {
	if (out_format == OutFileFormat::TXT) {
		pad_3d.save_as_txt(out_file_path, resolution);
	}
}

//反演LAI

void GridPADEstimator::lai_inverse() {
	for (int x = 0; x < matrix_dims[0]; x++) {
		for (int y = 0; y < matrix_dims[1]; y++) {
			for (int z = 0; z < matrix_dims[2]; z++) {
				double inc_energy = inc_3d.data[x][y][z];
				if (inc_energy != 0) {
					double trans = out_3d.data[x][y][z] / inc_energy;
					double avg_path_length = path_lengths_3d.data[x][y][z] / num_rays_3d.data[x][y][z];
					if(avg_path_length != 0)
						pad_3d.data[x][y][z] = -2 * std::log(trans) / avg_path_length;
				}

			}
		}
	}
}


void GridPADEstimator::init(std::vector<std::string> & input_files) {
	constexpr double inf = std::numeric_limits<double>::infinity();
	min_xyz = Eigen::Vector3d(inf, inf, inf);
	max_xyz = Eigen::Vector3d(-inf, -inf, -inf);
	if (m_file_type == "las") {
		for (int i = 0; i < input_files.size(); i++) {
			LASreadOpener lasreadopener;
			lasreadopener.set_file_name(input_files[i].c_str());
			LASreader* lasreader = lasreadopener.open();
			if (lasreader == 0)
			{
				fprintf(stderr, "ERROR: could not open lasreader\n");
				exit(0);
			}
			double min_x = lasreader->get_min_x();
			double min_y = lasreader->get_min_y();
			double min_z = lasreader->get_min_z();
			double max_x = lasreader->get_max_x();
			double max_y = lasreader->get_max_y();
			double max_z = lasreader->get_max_z();
			if (min_x < min_xyz[0]) { min_xyz[0] = min_x; }
			if (min_y < min_xyz[1]) { min_xyz[1] = min_y; }
			if (min_z < min_xyz[2]) { min_xyz[2] = min_z; }
			if (max_x > max_xyz[0]) { max_xyz[0] = max_x; }
			if (max_y > max_xyz[1]) { max_xyz[1] = max_y; }
			if (max_z > max_xyz[2]) { max_xyz[2] = max_z; }
			lasreader->close();
		}
	}
	int xsize = (int)ceil((max_xyz[0] - min_xyz[0]) / this->resolution[0]);
	int ysize = (int)ceil((max_xyz[1] - min_xyz[1]) / this->resolution[1]);
	int zsize = (int)ceil((max_xyz[2] - min_xyz[2]) / this->resolution[2]);
	pad_3d.resize(xsize, ysize, zsize);
	inc_3d.resize(xsize, ysize, zsize);
	out_3d.resize(xsize, ysize, zsize);
	num_rays_3d.resize(xsize, ysize, zsize);
	path_lengths_3d.resize(xsize, ysize, zsize);
	matrix_dims = Eigen::Vector3i(xsize, ysize, zsize);
}

void GridPADEstimator::ray_traversal_per_flight() {
	for (std::size_t i = 0; i < this->m_pulses.size(); i++) {
		Pulse* pulse = this->m_pulses[i];
		if (pulse->has_valid_nearest_point) {
			Eigen::Vector3d ray_dir = pulse->dir;
			//将点进行归一化，即减去原点
			Eigen::Vector3d ray_new_o = pulse->get_first_echo()->pos - RAY_O_OFFSET * ray_dir - min_xyz;
			//计算光线与矩阵顶部的交点
			double t0 = (max_xyz[2] - min_xyz[2] - ray_new_o[2]) / ray_dir[2];
			Eigen::Vector3d ray_o_top = ray_new_o + t0 * ray_dir;
			Eigen::Vector3d ray_new_end = pulse->get_last_echo()->pos - min_xyz;
			std::vector<double> traversed_lengths;
			std::vector<Eigen::Vector3i> visited_voxels = this->voxel_traversal(ray_o_top, ray_new_end, traversed_lengths);

			// Calculate which voxel contains a actual echo
			std::vector<Eigen::Vector3i> _echo_voxel_coords = pulse->get_voxel_coords_of_echos(min_xyz, resolution);
			std::vector<double> _echo_inc_intensites = pulse->get_echo_inc_intensities();
			// Remove multiple echoes in the same voxel
			std::vector<Eigen::Vector3i> echo_voxel_coords{ _echo_voxel_coords [0]};
			std::vector<double> echo_inc_intensites{ _echo_inc_intensites[0]};
			for (int k = 1; k < _echo_voxel_coords.size(); k++) {
				if (!(_echo_voxel_coords[k] == _echo_voxel_coords[k - 1])) {
					echo_voxel_coords.push_back(_echo_voxel_coords[k]);
					echo_inc_intensites.push_back(_echo_inc_intensites[k]);
				}
			}
			int correct_return_num = echo_voxel_coords.size();
			int cur_return_cell_index = 0;

			if (pulse->get_pulse_type() == PULSE_TYPE::VEG_GROUND) {
				for (int j = 0; j < visited_voxels.size() - 1; j++) {
					Eigen::Vector3i visited_voxel_coord = visited_voxels[j];
					if (visited_voxel_coord[0] >= 0 && visited_voxel_coord[0] < matrix_dims[0] &&
						visited_voxel_coord[1] >= 0 && visited_voxel_coord[1] < matrix_dims[1] &&
						visited_voxel_coord[2] >= 0 && visited_voxel_coord[2] < matrix_dims[2]) {
						//If the current visiting voxel contains an acutal echo
						if (visited_voxel_coord == echo_voxel_coords[cur_return_cell_index]) {
							inc_3d[visited_voxel_coord] += echo_inc_intensites[cur_return_cell_index];
							out_3d[visited_voxel_coord] += echo_inc_intensites[cur_return_cell_index + 1];
							cur_return_cell_index++;
						}
						else {
							inc_3d[visited_voxel_coord] += echo_inc_intensites[cur_return_cell_index];
							out_3d[visited_voxel_coord] += echo_inc_intensites[cur_return_cell_index];
						}
						path_lengths_3d[visited_voxel_coord] += traversed_lengths[j];
						num_rays_3d[visited_voxel_coord] += 1;
					}
				}
			}

			if (pulse->get_pulse_type() == PULSE_TYPE::PURE_VEG) {
				for (int j = 0; j < visited_voxels.size(); j++) {
					Eigen::Vector3i visited_voxel_coord = visited_voxels[j];
					if (visited_voxel_coord[0] >= 0 && visited_voxel_coord[0] < matrix_dims[0] &&
						visited_voxel_coord[1] >= 0 && visited_voxel_coord[1] < matrix_dims[1] &&
						visited_voxel_coord[2] >= 0 && visited_voxel_coord[2] < matrix_dims[2]) {
						//If the current visiting voxel contains an acutal echo
						if (visited_voxel_coord == echo_voxel_coords[cur_return_cell_index]) {
							inc_3d[visited_voxel_coord] += echo_inc_intensites[cur_return_cell_index];

							if (cur_return_cell_index < correct_return_num - 1) {
								out_3d[visited_voxel_coord] += echo_inc_intensites[cur_return_cell_index + 1];
							}
							cur_return_cell_index++;
						}
						else {
							inc_3d[visited_voxel_coord] += echo_inc_intensites[cur_return_cell_index];
							out_3d[visited_voxel_coord] += echo_inc_intensites[cur_return_cell_index];
						}
						path_lengths_3d[visited_voxel_coord] += traversed_lengths[j];
						num_rays_3d[visited_voxel_coord] += 1;
					}
				}
			}
		}
	}
}

std::vector<Eigen::Vector3i> GridPADEstimator::voxel_traversal(Eigen::Vector3d ray_start, Eigen::Vector3d ray_end, std::vector<double> &traversed_lengths) {
	std::vector<Eigen::Vector3i> visited_voxels;
	double res_x = this->resolution[0];
	double res_y = this->resolution[1];
	double res_z = this->resolution[2];
	// This id of the first/current voxel hit by the ray.
  // Using floor (round down) is actually very important,
  // the implicit int-casting will round up for negative numbers.
	Eigen::Vector3i current_voxel((int)std::floor(ray_start[0] / res_x),
		(int)std::floor(ray_start[1] / res_y),
		(int)std::floor(ray_start[2] / res_z));

	// The id of the last voxel hit by the ray.
	// TODO: what happens if the end point is on a border?
	Eigen::Vector3i last_voxel((int)std::floor(ray_end[0] / res_x),
		(int)std::floor(ray_end[1] / res_y),
		(int)std::floor(ray_end[2] / res_z));

	Eigen::Vector3d ray = ray_end - ray_start;
	ray.normalize();

	// In which direction the voxel ids are incremented.
	int stepX = (ray[0] >= 0) ? 1 : -1; // correct
	int stepY = (ray[1] >= 0) ? 1 : -1; // correct
	int stepZ = (ray[2] >= 0) ? 1 : -1; // correct

	// Distance along the ray to the next voxel border from the current position (tMaxX, tMaxY, tMaxZ).
	double next_voxel_boundary_x = (ray[0] > 0) ? (current_voxel[0] + 1)*res_x: current_voxel[0] * res_x; // correct
	double next_voxel_boundary_y = (ray[1] > 0) ? (current_voxel[1] + 1)*res_y: current_voxel[1] * res_y; // correct
	double next_voxel_boundary_z = (ray[2] > 0) ? (current_voxel[2] + 1)*res_z: current_voxel[2] * res_z; // correct

	// tMaxX, tMaxY, tMaxZ -- distance until next intersection with voxel-border
  // the value of t at which the ray crosses the first vertical voxel boundary
	double tMaxX = (ray[0] != 0) ? (next_voxel_boundary_x - ray_start[0]) / ray[0] : DBL_MAX; //
	double tMaxY = (ray[1] != 0) ? (next_voxel_boundary_y - ray_start[1]) / ray[1] : DBL_MAX; //
	double tMaxZ = (ray[2] != 0) ? (next_voxel_boundary_z - ray_start[2]) / ray[2] : DBL_MAX; //

	double t = std::min(tMaxX, std::min(tMaxY, tMaxZ));
	traversed_lengths.push_back(t);
	double pre_t = t;

	// tDeltaX, tDeltaY, tDeltaZ --
  // how far along the ray we must move for the horizontal component to equal the width of a voxel
  // the direction in which we traverse the grid
  // can only be FLT_MAX if we never go in that direction
	double tDeltaX = (ray[0] != 0) ? res_x / ray[0] * stepX : DBL_MAX;
	double tDeltaY = (ray[1] != 0) ? res_y / ray[1] * stepY : DBL_MAX;
	double tDeltaZ = (ray[2] != 0) ? res_z / ray[2] * stepZ : DBL_MAX;

	visited_voxels.push_back(current_voxel);

	while (last_voxel != current_voxel) {
		if (tMaxX < tMaxY) {
			if (tMaxX < tMaxZ) {
				current_voxel[0] += stepX;
				tMaxX += tDeltaX;
			}
			else {
				current_voxel[2] += stepZ;
				tMaxZ += tDeltaZ;
			}
		}
		else {
			if (tMaxY < tMaxZ) {
				current_voxel[1] += stepY;
				tMaxY += tDeltaY;
			}
			else {
				current_voxel[2] += stepZ;
				tMaxZ += tDeltaZ;
			}
		}
		visited_voxels.push_back(current_voxel);
		t = std::min(tMaxX, std::min(tMaxY, tMaxZ));
		traversed_lengths.push_back(t - pre_t);
		pre_t = t;
	}
	return visited_voxels;

}