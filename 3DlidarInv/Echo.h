#ifndef _ECHO_H_
#define _ECHO_H_

#include <Eigen/Core>

class Echo {
public:
	Eigen::Vector3d pos;
	double gps_time;
	unsigned char return_number;
	unsigned char number_of_returns;
	unsigned char classification;
	char scan_angle_rank;
	unsigned char user_data;
	unsigned char point_source_ID;
	double intensity;

	double inc_intensity; //根据最邻近纯地面intensity计算的incident能量
};

#endif // !_ECHO_H_

