#pragma once
#include <Eigen/Eigen>

class Area {
public:
    unsigned int ID;
    Eigen::RowVector4i Lines; // Lines which make up area boundary
    std::vector<unsigned int> Nodes;
};
