#pragma once
#include <Eigen/Eigen>

class Line {
public:
    unsigned int ID;
    Eigen::Matrix<unsigned int, 1, 2> Endpoints;
    unsigned int Divisions; // Number of elements to make along line
    float SpacingRatio = 1.0;
    std::vector<unsigned int> Nodes; // Store node references
    int Owner = -1;
    std::vector<unsigned int> Neighbours;
    std::vector<unsigned int> OverlappingLines;

};
