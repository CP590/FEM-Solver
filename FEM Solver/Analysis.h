#pragma once
#include <Eigen/Eigen>
#include <Tuple>
#include <vector>
#include "Material.h"

void InitialiseFMatrix();

void InitialiseuMatrix();

void InitialiseKMatrix();

void ApplyConstraintsToLine(unsigned int id, bool ConstrainX, bool ConstrainY);

void ApplyConstraintsToNode(unsigned int id, bool ConstrainX, bool ConstrainY);

unsigned short int ReturnConstraintValue(bool Constrain);

void PrintNodeConstraints();

void PrintNodeConstraints(unsigned int id);

void ApplyXForcesToLine(unsigned int id, double Force);

void DeleteXForcesOnLine(unsigned int id);

void ApplyXForcesToNode(unsigned int id, double Force);

void DeleteXForcesOnNode(unsigned int id);

void ApplyYForcesToLine(unsigned int id, double Force);

void DeleteYForcesOnLine(unsigned int id);

void ApplyYForcesToNode(unsigned int id, double Force);

void DeleteYForcesOnNode(unsigned int id);

void PrintNodeForces();

void PrintNodeForces(unsigned int id);

void SetMaterial(double E, double nu, double thickness);

void SolveSystem();

void TransferForcesToFMatrix();

void GenerateKMatrix();

// Assumed only isotropic for now
Eigen::MatrixXd GenerateDMatrix();

Eigen::MatrixXd GenerateDMatrixPlaneStress();

Eigen::MatrixXd GenerateDMatrixPlaneStrain();

Eigen::MatrixXd GenerateDMatrixIsotropic();

void GeneratekMatrix(Eigen::MatrixXd* D, Eigen::MatrixXd* k, int id);

void GenerateBMatrix(int id, int j, Eigen::MatrixXd* IntegrationPoints, Eigen::MatrixXd* B, double* Jdet);

void GenerateShapeFunctions(Eigen::MatrixXd* IntegrationPoints, Eigen::MatrixXd* ShapeFunctions);

// Returns integration points matrix depending on element type and evaluation points (Gauss or nodes)
Eigen::MatrixXd ReturnIntegrationPoints(int Locations);

void ReturnGaussIntegrationPoints(Eigen::MatrixXd* IntegrationPoints);

void ReturnExtrapolationIntegrationPoints(Eigen::MatrixXd* IntegrationPoints);

void ReturnNodeIntegrationPoints(Eigen::MatrixXd* IntegrationPoints);

void PopulateKMatrix(int id, Eigen::MatrixXd* k);

std::tuple<Eigen::VectorXd, std::vector<int>> GenerateReducedSystem();

std::vector<int> TransferConstraints();

int SumConstraints();

void ReduceMatrices(std::vector<int>* DeletedIndices, int* FinalMatrixSize);

std::vector<int> GenerateKeptIndices(std::vector<int>* DeletedIndices);

void PopulateDisplacements(Eigen::VectorXd* uReduced, std::vector<int>* KeptIndices);

void PrintNodeDisplacements();

void PrintNodeDisplacements(unsigned int id);

void PrintElementDisplacements();

void PrintElementDisplacements(unsigned int id);

void GenerateNodeStrainsAndStresses();

void GenerateExtrapolatedStrainsAndStresses();

void GetElementDisplacements(unsigned int id, Eigen::VectorXd* Displacements);

void GenerateElementNodeStrainsAndStresses(unsigned int id, Eigen::MatrixXd* NodeIntegrationPoints, Eigen::MatrixXd* D, Eigen::VectorXd* Displacements, int MatrixSize, Eigen::VectorXd* Strains, Eigen::VectorXd* Stresses);

void GenerateElementExtrapolatedStrainsAndStresses(unsigned int id, Eigen::MatrixXd* GaussIntegrationPoints, Eigen::MatrixXd* ShapeFunctions, Eigen::MatrixXd* D, Eigen::VectorXd* Displacements, int MatrixSize, Eigen::MatrixXd* GaussStrains, Eigen::MatrixXd* GaussStresses, Eigen::MatrixXd* ExtrapolatedStrains, Eigen::MatrixXd* ExtrapolatedStresses);

void PrintElementStrains();

void PrintElementStrains(unsigned int id);

void PrintElementStresses();

void PrintElementStresses(unsigned int id);
