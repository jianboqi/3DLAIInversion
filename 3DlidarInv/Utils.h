#ifndef _UTILS_H_
#define _UTILS_H_
#include <algorithm>
#include "delaunator.hpp"
#include <vector>
#include <string>

//output triangle mesh to obj file
void output_obj(delaunator::Delaunator d, std::string file_path);

template<typename InputIterator, typename ValueType>
InputIterator closest(InputIterator first, InputIterator last, ValueType value)
{
	return std::min_element(first, last, [&](ValueType x, ValueType y)
	{
		return std::abs(x - value) < std::abs(y - value);
	});
}

template<typename T>
class Ray {
	T org[3];        // [in] must set
	T dir[3];        // [in] must set
	T min_t;         // [in] must set
	T max_t;         // [in] must set
	unsigned int type;  // optional. ray type.
};

class BVHTraceOptions {
	// Trace rays only in face ids range. faceIdsRange[0] < faceIdsRange[1]
	// default: 0 to 0x3FFFFFFF(2G faces)
	unsigned int prim_ids_range[2];
	bool cull_back_face; // default: false
};

/// <summary>
/// 对输入的文件进行解析，如何含有*，则进行通配符匹配
/// </summary>
/// <param name="input_files"></param>
/// <returns></returns>
//std::vector<std::string> parse_input_files(std::vector<std::string> input_files) {
//
//}

#endif // !_UTILS_H_

