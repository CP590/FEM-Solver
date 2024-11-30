#pragma once
#include <Eigen/Eigen>

class Element {
public:
	unsigned int ID;
	Eigen::RowVectorXi Nodes;
	Eigen::MatrixXd Strains;
	Eigen::MatrixXd Stresses;
};
