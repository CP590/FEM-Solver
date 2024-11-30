#pragma once
#include <Eigen/Eigen>

class Node {
public:
	unsigned int ID;
	Eigen::RowVector2d Position;
	Eigen::RowVector2i Constraints;
	Eigen::RowVector2d Forces;
	Eigen::RowVector2d Displacements;
};
