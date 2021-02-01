// Author: Jianbo Qi
// Date: 2020/7/15
// This is for handling lidar data for inverting 3D PAD

#include "argparse.hpp"
#include <string>
#include <iostream>
#include <lasreader.hpp>
#include <iomanip>
#include "PADEstimator.h"
#include "GridPADEstimator.h"
#include "Echo.h"

int main(int argc, const char** argv) {
	ArgumentParser parser;
	parser.addArgument("-i", "--input", '+', false);
	parser.addArgument("-o", "--output", 1, false);
	parser.addArgument("--res_x", 1);
	parser.addArgument("--res_y", 1);
	parser.addArgument("--res_z", 1);
	parser.addArgument("-f", "--filetype", 1);
	parser.parse(argc, argv);

	/*处理参数*/
	std::vector<std::string> input_files = parser.retrieve<std::vector<std::string> >("input");
	std::string out_file = parser.retrieve<std::string>("output");

	double res_x=1, res_y=1, res_z=1;
	if (parser.exists("res_x")) {
		res_x = atof(parser.retrieve<std::string>("res_x").c_str());
	}

	if (parser.exists("res_y")) {
		res_y = atof(parser.retrieve<std::string>("res_y").c_str());
	}

	if (parser.exists("res_z")) {
		res_z = atof(parser.retrieve<std::string>("res_z").c_str());
	}

	std::cout << " - Resolutions: x=" << res_x <<" y="<<res_y<<" z="<<res_z<< std::endl;

	std::string file_type = "las";
	if (parser.exists("filetype")) {
		file_type = parser.retrieve<std::string>("filetype");
	}
	
	std::cout << " - File Type: " << file_type << std::endl;

	EstimatorType estimator_type = EstimatorType::GRID;

	/*开始计算*/
	PADEstimator* pad_estimator=NULL;
	if (estimator_type == EstimatorType::GRID) {
		pad_estimator = new GridPADEstimator(res_x, res_y, res_z);
		pad_estimator->set_file_type(file_type);
	}
	pad_estimator->pad_inverse(input_files);
	pad_estimator->save_to_file(out_file, OutFileFormat::TXT);
	
	return 0;
}