#ifndef _PULSE_H_
#define _PULSE_H_

#include <vector>
#include "Echo.h"
#include<algorithm>

enum class PULSE_TYPE { VEG_GROUND, PURE_VEG, PURE_GROUND};

bool echo_cmp(Echo* echo1, Echo* echo2);

class Pulse {
public:
	std::vector<Echo*> m_echo_list;
	PULSE_TYPE m_pulse_type;
	Eigen::Vector3d dir;
	//脉冲是否具有最邻近的点，在数据边缘附近，脉冲经过光线追踪后与地面交点可能超出边界
	bool has_valid_nearest_point = true;
	//double nearest_intensity;

	void add_echo(Echo* echo) {
		this->m_echo_list.push_back(echo);
	}

	Echo* get_last_echo() {
		return m_echo_list[m_echo_list.size() - 1];
	}

	Echo* get_first_echo() {
		return m_echo_list[0];
	}

	std::vector<Echo*> get_echo_list() {
		return m_echo_list;
	}

	//按return number降序排序
	void sort_according_return_num() {
		std::sort(m_echo_list.begin(), m_echo_list.end(), echo_cmp);
	}

	double get_gps_time() {
		return m_echo_list[0]->gps_time;
	}

	char get_scan_angle_rank() {
		return m_echo_list[0]->scan_angle_rank;
	}

	PULSE_TYPE get_pulse_type() {
		return this->m_pulse_type;
	}

	double get_cumulated_intensity_from(int start_from) {
		double tmp=0;
		for (int i = start_from; i < m_echo_list.size(); i++) {
			tmp += m_echo_list[i]->intensity;
		}
		return tmp;
	}

	//不包含to
	double get_cumulated_intensity_to(int to) {
		double tmp=0;
		for (int i = 0; i < to; i++) {
			tmp += m_echo_list[i]->intensity;
		}
		return tmp;
	}

	//echo.classification are always 1 or 2
	PULSE_TYPE determine_pulse_type();

	//whether pulse is invalid
	bool is_valid();

	/// <summary>
	/// Get the voxel coordinates of echo
	/// </summary>
	/// <param name="xyz_min"></param>
	/// <param name="resolution"></param>
	/// <returns></returns>
	std::vector<Eigen::Vector3i> get_voxel_coords_of_echos(Eigen::Vector3d xyz_min, Eigen::Vector3d resolution) {
		std::vector<Eigen::Vector3i> voxel_coords;
		for (int i = 0; i < m_echo_list.size(); i++) {
			Echo* echo = m_echo_list[i];
			voxel_coords.push_back(Eigen::Vector3i((int)std::floor((echo->pos.x()- xyz_min.x()) / resolution.x()),
				(int)std::floor((echo->pos.y() - xyz_min.y()) / resolution.y()),
				(int)std::floor((echo->pos.z() - xyz_min.z()) / resolution.z())));
		}
		return voxel_coords;
	}

	/// <summary>
	/// Get the incident intensity at each echo
	/// </summary>
	/// <returns></returns>
	std::vector<double> get_echo_inc_intensities() {
		std::vector<double> intensities;
		for (int i = 0; i < m_echo_list.size(); i++) {
			Echo* echo = m_echo_list[i];
			intensities.push_back(echo->inc_intensity);
		}
		return intensities;
	}

};



#endif // !_PULSE_H_

