#ifndef _PADESTIMATOR_H_
#define _PADESTIMATOR_H_

#include "Pulse.h"
#include <vector>
#include <map>
#include <string>
#include "Utils.h"
#include "nanoflann.hpp"
#include "KDTreeVectorOfVectorsAdaptor.h"
#include "delaunator.hpp"
#include<stdexcept>
#include <iostream>
#include <iomanip>

enum class EstimatorType{
	GRID, 
	ALPHA_SHAPE
};

enum class OutFileFormat {
	TXT,
	NPY
};


class PADEstimator {
public:
	std::vector<Pulse*> m_pulses; //临时储存每个las文件（或者txt）中的pulse
	
	std::string m_file_type;
	//long int
public:

	void set_file_type(std::string ft) {this->m_file_type = ft;}

	// invert for each flight
	//该函数由子类继承
	void virtual pad_inverse(std::vector<std::string> & input_files);
	void virtual save_to_file(std::string out_file_path, OutFileFormat out_format = OutFileFormat::TXT);

	//parse flight information
	void pulse_info_per_flight(std::string input_file);

private:

	//recover pulses from point cloud
	void parse_pulses_from_file(std::string input_file);

	//recover pulses from las point cloud
	void parse_pulses_from_point_cloud_las(std::string las_file);

	//determine pulse direction
	void determine_pulse_dir();

	//determine pulse type: pure-ground, pure vegetation...
	void determine_pulse_type();

	//determine pulse echo incident energy
	void determine_nearest_pure_ground_intensity();

	void determine_nearest_pure_ground_intensity_with_nanort();

};

#endif // !_PADESTIMATOR_H_

