#pragma once
#include <Eigen/Eigen>

class Keypoint {
public:
	unsigned int ID;
	Eigen::RowVector2d Position;
};
