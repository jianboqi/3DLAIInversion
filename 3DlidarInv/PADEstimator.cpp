#include "PADEstimator.h"
#include <lasreader.hpp>
#include "nanort.h"
#include <fstream>
#include "GridPADEstimator.h"

//#define _DEBUG_MODE
//#define _DEBUG_MODE_OUT_OBJ
#ifdef _DEBUG_MODE
#define DEBUG_DIR "D:\\TMP\\"
#endif // _DEBUG_MODE



void PADEstimator::pad_inverse(std::vector<std::string> & input_files) {}
void PADEstimator::save_to_file(std::string out_file_path, OutFileFormat out_format) {}

void PADEstimator::pulse_info_per_flight(std::string input_file) {
	this->parse_pulses_from_file(input_file);
	this->determine_pulse_type();
	this->determine_pulse_dir();
	this->determine_nearest_pure_ground_intensity();
}

void PADEstimator::parse_pulses_from_file(std::string input_file) {
	m_pulses.clear();
	if (m_file_type == "las") {
		this->parse_pulses_from_point_cloud_las(input_file);
	}
	else {

	}
}

/// <summary>
/// 计算每个脉冲的入射能量
/// </summary>
void PADEstimator::determine_nearest_pure_ground_intensity() {
	//对每个植被-地面脉冲，找一个最邻近的地面脉冲

	//构建所有纯地面脉冲的kdtree
	std::vector<std::vector<double>> pt2d;
	std::vector<double> pure_ground_intensity_list;

	for (int i = 0; i < m_pulses.size(); i++) {
		if (m_pulses[i]->get_pulse_type() == PULSE_TYPE::PURE_GROUND) {
			pt2d.push_back(std::vector<double>{m_pulses[i]->m_echo_list[0]->pos.x(),
				m_pulses[i]->m_echo_list[0]->pos.y()});
			pure_ground_intensity_list.push_back(m_pulses[i]->get_last_echo()->intensity);
		}
	}

	typedef KDTreeVectorOfVectorsAdaptor< std::vector<std::vector<double>>, double >  my_kd_tree_t;
	my_kd_tree_t mat_index(2, pt2d, 10);
	mat_index.index->buildIndex();
	// do a knn search
	const size_t num_results = 1;
	std::vector<size_t>   ret_indexes(num_results);
	std::vector<double> out_dists_sqr(num_results);

	std::size_t num_pulses_without_np = 0;
	std::size_t number_of_pulses_invalid_intensity = 0;

	for (int i = 0; i < m_pulses.size(); i++) {
		//对于纯植被脉冲，首先需要计算其与地面的交点，再计算最邻近脉冲
		if (m_pulses[i]->get_pulse_type() == PULSE_TYPE::VEG_GROUND) {
				std::vector<double> query_pt{ m_pulses[i]->get_last_echo()->pos[0], m_pulses[i]->get_last_echo()->pos[1] };
				mat_index.index->knnSearch(&query_pt[0], 1, &ret_indexes[0], &out_dists_sqr[0]);
				double qstar = pure_ground_intensity_list[ret_indexes[0]];
				double qg = m_pulses[i]->get_last_echo()->intensity;
				if (qstar <= qg) {
					m_pulses[i]->has_valid_nearest_point = false;
					number_of_pulses_invalid_intensity++;
					for (int j = 0; j < m_pulses[i]->get_echo_list().size(); j++) {
						m_pulses[i]->get_echo_list()[j]->inc_intensity = m_pulses[i]->get_cumulated_intensity_from(j);
					}
				}
				else {
					for (int j = 0; j < m_pulses[i]->get_echo_list().size(); j++) {
						m_pulses[i]->get_echo_list()[j]->inc_intensity = (qg * m_pulses[i]->get_cumulated_intensity_to(j) +
							qstar * m_pulses[i]->get_cumulated_intensity_from(j)) / (qstar - qg);
					}
				}
		}
		else if (m_pulses[i]->get_pulse_type() == PULSE_TYPE::PURE_VEG) {
			for (int j = 0; j < m_pulses[i]->get_echo_list().size(); j++) {
				m_pulses[i]->get_echo_list()[j]->inc_intensity = m_pulses[i]->get_cumulated_intensity_from(j);
			}
		}
	}
	std::cout << "   - Number of pulses without nearest pure-ground points: " << num_pulses_without_np << " -> "
		<< std::setprecision(4) << (100.0 * num_pulses_without_np) / m_pulses.size() << "%" << std::endl;

	std::cout << "   - Number of pulses with invalid nearest pure-ground points: " << number_of_pulses_invalid_intensity << " -> "
		<< std::setprecision(4) << (100.0 * number_of_pulses_invalid_intensity) / m_pulses.size() << "%" << std::endl;

}

/// <summary>
/// 计算每个脉冲的入射能量
/// 此版本使用了三角网化，求交点，对于计算冠层总体透过率，需要计算纯
/// 植被脉冲与地面的交点
/// </summary>
void PADEstimator::determine_nearest_pure_ground_intensity_with_nanort() {
	//对每个植被-地面脉冲，找一个最邻近的地面脉冲

	//构建所有纯地面脉冲的kdtree
	std::vector<std::vector<double>> pt2d;
	std::vector<double> pure_ground_intensity_list;
	std::vector<double> xy_coords;
	std::vector<double> xyz_coords;

#ifdef _DEBUG_MODE
	std::ofstream ofs;
	ofs.open((std::string)DEBUG_DIR + "pure_ground.txt", std::ios::out);
#endif // _DEBUG_MODE

	for (int i = 0; i < m_pulses.size(); i++) {
		if (m_pulses[i]->get_pulse_type() == PULSE_TYPE::PURE_GROUND) {
			pt2d.push_back(std::vector<double>{m_pulses[i]->m_echo_list[0]->pos.x(),
				m_pulses[i]->m_echo_list[0]->pos.y()});
			pure_ground_intensity_list.push_back(m_pulses[i]->get_last_echo()->intensity);
#ifdef _DEBUG_MODE
				ofs << std::setprecision(15) << m_pulses[i]->m_echo_list[0]->pos.x() << " " << std::setprecision(15) << m_pulses[i]->m_echo_list[0]->pos.y() << std::endl;
#endif // _DEBUG_MODE
		}

		if (m_pulses[i]->get_pulse_type() == PULSE_TYPE::VEG_GROUND ||
			m_pulses[i]->get_pulse_type() == PULSE_TYPE::PURE_GROUND) {
			double x = m_pulses[i]->get_last_echo()->pos.x();
			double y = m_pulses[i]->get_last_echo()->pos.y();
			double z = m_pulses[i]->get_last_echo()->pos.z();
			xy_coords.push_back(x);
			xy_coords.push_back(y);
			xyz_coords.push_back(x);
			xyz_coords.push_back(y);
			xyz_coords.push_back(z);
		}
	}
#ifdef _DEBUG_MODE
	ofs.close();
#endif // _DEBUG_MODE

	
	typedef KDTreeVectorOfVectorsAdaptor< std::vector<std::vector<double>>, double >  my_kd_tree_t;
	my_kd_tree_t mat_index(2, pt2d, 10);
	mat_index.index->buildIndex();

	//用三角化的方式，首先做三角网,然后用nanort计算交点
	delaunator::Delaunator d(xy_coords);
#ifdef _DEBUG_MODE_OUT_OBJ
	output_obj(d, "D:\\terrain_mesh.obj");
#endif // _DEBUG_MODE_OUT_OBJ

	const std::vector<unsigned int> triange_indices(d.triangles.begin(), d.triangles.end());
	nanort::TriangleMesh<double> triangle_mesh(xyz_coords.data(), triange_indices.data(), /* stride */sizeof(double) * 3);
	nanort::TriangleSAHPred<double> triangle_pred(xyz_coords.data(), triange_indices.data(), /* stride */sizeof(double) * 3);
	nanort::TriangleIntersector<double, nanort::TriangleIntersection<double> > triangle_intersecter(xyz_coords.data(), triange_indices.data(), /* stride */sizeof(double) * 3);
	nanort::BVHAccel<double> accel;
	nanort::BVHBuildOptions<double> build_options; // Use default option
	std::cout << " - Building data structure for ray-terrain intersection..."<<std::endl;
	bool ret = accel.Build((unsigned int)(d.triangles.size()/3.0), triangle_mesh, triangle_pred, build_options);
	std::cout << " - Finding nearest ground points for pure vegetation pulses..." << std::endl;
	if (ret) {
		// do a knn search
		const size_t num_results = 1;
		std::vector<size_t>   ret_indexes(num_results);
		std::vector<double> out_dists_sqr(num_results);

#ifdef _DEBUG_MODE
		std::ofstream ofs1;
		ofs1.open((std::string)DEBUG_DIR + "terrain_intersected.txt", std::ios::out);
		std::ofstream ofs2;
		ofs2.open((std::string)DEBUG_DIR + "found_nearest.txt", std::ios::out);
		std::ofstream ofs3;
		ofs3.open((std::string)DEBUG_DIR + "Pulse_without_np.txt");
#endif // _DEBUG_MODE

		std::size_t num_pulses_without_np = 0;
		std::size_t number_of_pulses_invalid_intensity = 0;

		for (int i = 0; i < m_pulses.size(); i++) {
			//对于纯植被脉冲，首先需要计算其与地面的交点，再计算最邻近脉冲
			if (m_pulses[i]->get_pulse_type() == PULSE_TYPE::VEG_GROUND) {			
				nanort::TriangleIntersection<double> isect;
				nanort::Ray<double> ray;
				ray.org[0] = m_pulses[i]->get_first_echo()->pos[0];
				ray.org[1] = m_pulses[i]->get_first_echo()->pos[1];
				ray.org[2] = m_pulses[i]->get_first_echo()->pos[2];
				ray.dir[0] = m_pulses[i]->dir[0];
				ray.dir[1] = m_pulses[i]->dir[1];
				ray.dir[2] = m_pulses[i]->dir[2];
				ray.min_t = 0.0f;
				ray.max_t = 1.0e+30f;
				bool hit = accel.Traverse(ray, triangle_intersecter, &isect);
				if (hit) {
					std::vector<double> query_pt{ ray.org[0]+ ray.dir[0]* isect.t,
					ray.org[1] + ray.dir[1] * isect.t };
					std::cout << "....................................................." << std::endl;
					std::cout << " - Ground point-hit: " <<std::setprecision(15)<< query_pt[0] << " " << std::setprecision(15) << query_pt[1] << std::endl;
					std::cout << " - Ground point-lat: " << std::setprecision(15) << m_pulses[i]->get_last_echo()->pos[0] << " " << std::setprecision(15) << m_pulses[i]->get_last_echo()->pos[1] << std::endl;
					mat_index.index->knnSearch(&query_pt[0], 1, &ret_indexes[0], &out_dists_sqr[0]);
					double qstar = pure_ground_intensity_list[ret_indexes[0]];
					double qg = m_pulses[i]->get_last_echo()->intensity;
					if (qstar <= qg) {
						m_pulses[i]->has_valid_nearest_point = false;
						number_of_pulses_invalid_intensity++;
						for (int j = 0; j < m_pulses[i]->get_echo_list().size(); j++) {
							m_pulses[i]->get_echo_list()[j]->inc_intensity = m_pulses[i]->get_cumulated_intensity_from(j);
						}
					}
					else {
						for (int j = 0; j < m_pulses[i]->get_echo_list().size(); j++) {
							m_pulses[i]->get_echo_list()[j]->inc_intensity = (qg * m_pulses[i]->get_cumulated_intensity_to(j) +
								qstar * m_pulses[i]->get_cumulated_intensity_from(j))/(qstar - qg);
						}
					}
					 

#ifdef _DEBUG_MODE
					ofs1 << std::setprecision(15) << query_pt[0] << " " << std::setprecision(15) << query_pt[1] << std::endl;
					ofs2 << std::setprecision(15) << pt2d[ret_indexes[0]][0] << " " << std::setprecision(15) << pt2d[ret_indexes[0]][1] << std::endl;
#endif
					//std::cout << "ret_index[" << 0 << "]=" << ret_indexes[0] << " out_dist_sqr=" << std::sqrt(out_dists_sqr[0]) 
					//	<<" NP: " << std::setprecision(15) << pt2d[ret_indexes[0]][0]<<" " << std::setprecision(15) << pt2d[ret_indexes[0]][1] << std::endl;
				}
				else {
#ifdef _DEBUG_MODE
					ofs3 << std::setprecision(15) << ray.org[0] << " " << std::setprecision(15) << ray.org[1] << std::endl;
#endif
					m_pulses[i]->has_valid_nearest_point = false;
					num_pulses_without_np++;
				}
			}
			else if (m_pulses[i]->get_pulse_type() == PULSE_TYPE::PURE_VEG) {
				for (int j = 0; j < m_pulses[i]->get_echo_list().size(); j++) {
					m_pulses[i]->get_echo_list()[j]->inc_intensity = m_pulses[i]->get_cumulated_intensity_from(j);
				}
			}
		}
		std::cout << "   - Number of pulses without nearest pure-ground points: " << num_pulses_without_np << " -> "
			<< std::setprecision(4) << (100.0*num_pulses_without_np) / m_pulses.size() << "%" << std::endl;

		std::cout << "   - Number of pulses with invalid nearest pure-ground points: " << number_of_pulses_invalid_intensity << " -> "
			<< std::setprecision(4) << (100.0*number_of_pulses_invalid_intensity) / m_pulses.size() << "%" << std::endl;;

#ifdef _DEBUG_MODE
		ofs1.close();
		ofs2.close();
		ofs3.close();
#endif // _DEBUG_MODE

	}
	else {
		std::cout << " - Error: Building data structure for ray-terrain intersection failed." << std::endl;
		exit(0);
	}
	
}


void PADEstimator::determine_pulse_dir() {
	std::map<char, std::vector< Eigen::Vector3d>> scan_angle_directions;
	std::map<char, Eigen::Vector3d> scan_angle_direction_dic;
	for (int i = 0; i < m_pulses.size(); i++) {
		std::vector<Echo*> point_list = m_pulses[i]->m_echo_list;
		if (point_list.size() > 1) {  //具有多个回波点的脉冲可直接计算方向
			m_pulses[i]->dir = point_list[point_list.size() - 1]->pos
				- point_list[0]->pos;
			m_pulses[i]->dir.normalize();

			char scan_angle_rank = m_pulses[i]->get_scan_angle_rank();
			if (scan_angle_directions.count(scan_angle_rank == 0)) {
				scan_angle_directions[scan_angle_rank] = std::vector< Eigen::Vector3d>{ m_pulses[i]->dir };
			}
			else {
				scan_angle_directions[scan_angle_rank].push_back(m_pulses[i]->dir);
			}
		}
	}
	std::map<char, std::vector< Eigen::Vector3d>>::iterator iter;
	for (iter = scan_angle_directions.begin(); iter != scan_angle_directions.end(); iter++) {
		Eigen::Vector3d avg_dir(0, 0, 0);
		for (int i = 0; i < iter->second.size(); i++) {
			avg_dir += iter->second[i];
		}
		avg_dir /= iter->second.size();
		avg_dir.normalize();
		scan_angle_direction_dic[iter->first] = avg_dir;
	}

	//确定其他pulse的方向
	int pulse_with_no_dir = 0;
	for (int i = 0; i < m_pulses.size(); i++) {
		std::vector<Echo*> point_list = m_pulses[i]->m_echo_list;
		char scan_angle_rank = m_pulses[i]->get_scan_angle_rank();
		if (point_list.size() <= 1) {  //只有一个回波点
			if (scan_angle_direction_dic.count(scan_angle_rank) > 0) {
				m_pulses[i]->dir = scan_angle_direction_dic[scan_angle_rank];
			}
			else {
				pulse_with_no_dir++;
				m_pulses[i]->dir = Eigen::Vector3d(0, 0, -1);
			}
		}
	}
	std::cout << "   - Pulse with no direction: " << pulse_with_no_dir << " -> "
		<< std::setprecision(4) << (100.0*pulse_with_no_dir) / m_pulses.size() << "%" << std::endl;

}


void PADEstimator::determine_pulse_type() {
	long int num_pv = 0, num_pg = 0, num_vg = 0;
	for (int i = 0; i < m_pulses.size(); i++) {
		PULSE_TYPE pulse_type = m_pulses[i]->determine_pulse_type();
		switch (pulse_type)
		{
		case PULSE_TYPE::VEG_GROUND:
			num_vg++;
			break;
		case PULSE_TYPE::PURE_VEG:
			num_pv++;
			break;
		case PULSE_TYPE::PURE_GROUND:
			num_pg++;
			break;
		default:
			break;
		}
	}
	std::cout << "     - Pure-Vegetation Pulse: " << num_pv << " -> "
		<< std::setprecision(4) << (100.0*num_pv) / m_pulses.size() << "%" << std::endl;
	std::cout << "     - Pure-Ground Pulse: " << num_pg << " -> "
		<< std::setprecision(4) << (100.0*num_pg) / m_pulses.size() << "%" << std::endl;
	std::cout << "     - Vegetation-Ground Pulse: " << num_vg << " -> "
		<< std::setprecision(4) << (100.0*num_vg) / m_pulses.size() << "%" << std::endl;
}


void PADEstimator::parse_pulses_from_point_cloud_las(std::string las_file) {
	long int num_error_pulses = 0;
	LASreadOpener lasreadopener;
	lasreadopener.set_file_name(las_file.c_str());
	LASreader* lasreader = lasreadopener.open();
	if (lasreader == 0)
	{
		fprintf(stderr, "ERROR: could not open lasreader\n");
	}
	std::map<double, Pulse*> tmp_pulses;
	long int tot_points = 0;
	long int error_points = 0;
	while (lasreader->read_point())
	{
		Echo *echo = new Echo();
		echo->pos[0] = lasreader->point.quantizer->get_x(lasreader->point.X);
		echo->pos[1] = lasreader->point.quantizer->get_y(lasreader->point.Y);
		echo->pos[2] = lasreader->point.quantizer->get_z(lasreader->point.Z);
		echo->gps_time = lasreader->point.get_gps_time();
		echo->return_number = lasreader->point.get_return_number();
		echo->number_of_returns = lasreader->point.get_number_of_returns();
		echo->classification = lasreader->point.get_classification();
		echo->scan_angle_rank = lasreader->point.get_scan_angle_rank();
		echo->user_data = lasreader->point.get_user_data();
		echo->point_source_ID = lasreader->point.get_point_source_ID();
		echo->intensity = lasreader->point.get_intensity();

		double gps_time = echo->gps_time;
		if (tmp_pulses.count(gps_time) == 0) {
			tmp_pulses[gps_time] = new Pulse();
			tmp_pulses[gps_time]->add_echo(echo);
		}
		else {
			tmp_pulses[gps_time]->add_echo(echo);
		}
		tot_points++;
	}
	std::map<double, Pulse*>::iterator iter;
	for (iter = tmp_pulses.begin(); iter != tmp_pulses.end(); iter++) {
		if (iter->second->m_echo_list.size() > 1) {
			iter->second->sort_according_return_num();
		}
		if (iter->second->is_valid()) {
			m_pulses.push_back(iter->second);
		}
		else {
			num_error_pulses++;
			error_points += iter->second->m_echo_list.size();
		}
	}
	std::cout << " - Total number of points: " << tot_points << std::endl;
	std::cout << "   - Number of points from incomplete pulses: " << error_points
		<< " -> " << std::setprecision(4) << (100.0*error_points) / tot_points << "%" << std::endl;
	std::cout << " - Total number of pulses: " << tmp_pulses.size() << std::endl;
	std::cout << "   - Number of error pulses: " << num_error_pulses << " -> "
		<< std::setprecision(4) << (100.0*num_error_pulses) / tmp_pulses.size() << "%" << std::endl;
	std::cout << "   - Number of valid pulses: " << m_pulses.size() << " -> "
		<< std::setprecision(4) << (100.0* m_pulses.size()) / tmp_pulses.size() << "%" << std::endl;

}