#pragma once
#include <iostream>
#include <Eigen/Eigen>
#include <deque>
#include <vector>
#include "Keypoint.h"
#include "Line.h"
#include "Area.h"
#include "Node.h"
#include "Element.h"

namespace Geometry {

	std::deque<Keypoint> CreateKeypointsArray();
	std::deque<Line> CreateLinesArray();
	std::deque<Area> CreateAreasArray();
	std::deque<Node> CreateNodesArray();
	std::deque<Element> CreateElementsArray();

	std::vector<unsigned int> CreateDeletedKeypointIDs();
	std::vector<unsigned int> CreateDeletedLineIDs();
	std::vector<unsigned int> CreateDeletedAreaIDs();
	std::vector<unsigned int> CreateDeletedNodeIDs();
	std::vector<unsigned int> CreateDeletedElementIDs();

	std::vector<Keypoint*> CreateKeypointPointers();
	std::vector<Line*> CreateLinePointers();
	std::vector<Area*> CreateAreaPointers();
	std::vector<Node*> CreateNodePointers();
	std::vector<Element*> CreateElementPointers();

	//// KEYPOINTS ////
	void CreateKeypoint(double x, double y);

	bool KeypointExistsInArray(double x, double y);

	bool KeypointExists(int id, double x, double y);

	unsigned int AssignKeypointID();

	void UpdateNextNewKeypointID();

	void AssignKeypointPointer(unsigned int id);

	void DeleteKeypoint(unsigned int id);

	void AddDeletedKeypointID(unsigned int id);

	void AdjustKeypointPointers(unsigned int id, int DistanceToStart, int DistanceToEnd);

	void PrintKeypoints();

	void PrintKeypoints(unsigned int id);

	//// LINES ////
	void CreateLine(unsigned int Point0, unsigned int Point1);

	bool LineExistsInArray(unsigned int Keypoint0, unsigned int Keypoint1);

	bool LineExists(int id, unsigned int Keypoint0, unsigned int Keypoint1);

	unsigned int AssignLineID();

	void UpdateNextNewLineID();

	void AssignLinePointer(unsigned int id);

	void DeleteLine(unsigned int id);

	void AddDeletedLineID(unsigned int id);

	void AdjustLinePointers(unsigned int id, int DistanceToStart, int DistanceToEnd);

	void AssignLineDivision(unsigned int id, unsigned int div);

	void AssignLineSpacingRatio(unsigned int id, float SpacingRatio);

	void PrintLines();

	void PrintLines(unsigned int id);

	void PrintLineNodes();

	void PrintLineNodes(unsigned int id);

	//// AREAS ////
	void Create4SidedArea(unsigned int Line0, unsigned int Line1, unsigned int Line2, unsigned int Line3);

	bool AreaExistsInArray(std::vector<int>* PotentialLines);

	std::vector<int> PackageExistingLinesOfArea(int id);

	// Checks if area of n lines only consist of n unique keypoints
	bool AreaKeypointsUnique(std::vector<int>* PotentialLines);

	std::vector<int> PackageKeypointsOfPotentialLines(std::vector<int>* PotentialLines);

	Eigen::RowVectorXi OrderAreaLines(std::vector<int>* AreaLines);

	std::tuple<std::vector<double>, std::vector<double>> CalculateAveragePositionsOfAreaLines(std::vector<int>* AreaLines);

	unsigned int AssignAreaID();

	void UpdateNextNewAreaID();

	void AssignAreaPointer(unsigned int id);

	void AssignLineOwnersAndNeighbours(unsigned int id);

	void AssignAreaNeighbours(unsigned int NewArea);

	void AreasShareLine(unsigned int Area0, unsigned int Area1);

	void AreaLinesOverlap(unsigned int Area0, unsigned int Area1);

	bool LinesCollinear(unsigned int Line0, unsigned int Line1);

	void CheckIfLinesShareKeypoint(Eigen::VectorXi* SharedKeypointIndices, unsigned int Line0, unsigned int Line1);

	int ReturnOtherEndpoint(unsigned int id);

	std::tuple<std::vector<double>, std::vector<double>> PackageKeypointPositions(unsigned int Line);

	bool SharedKeypointInRangeOfLine(std::vector<double> KeypointXPositions, std::vector<double> KeypointYPositions, unsigned int Keypoint);

	bool KeypointInRangeOfLine(std::vector<double> KeypointXPositions, std::vector<double> KeypointYPositions, unsigned int Keypoint);

	bool SharedXPosition(unsigned int Keypoint, std::vector<double>* KeypointXPositions);

	bool SharedYPosition(unsigned int Keypoint, std::vector<double>* KeypointYPositions);

	void DeleteArea(unsigned int id);

	void AddDeletedAreaID(unsigned int id);

	void AdjustAreaPointers(unsigned int id, int DistanceToStart, int DistanceToEnd);

	void PrintAreaLines();

	void PrintAreaLines(unsigned int id);

	void PrintAreaNodes();

	void PrintAreaNodes(unsigned int id);

	void GenerateAreaMesh(unsigned int id);

	bool AreaLineDivisionsValid(unsigned int id);

	Eigen::MatrixXi PackageAreaLineEndpoints(unsigned int id);

	std::vector<Eigen::MatrixXd> CreateAreaLineNodes(unsigned int id, Eigen::MatrixXi* AreaLines);

	Eigen::RowVectorXd PackageAreaLineEndpointPositions(Eigen::RowVectorXi AreaLine);

	bool AreaLineGradientFinite(Eigen::RowVectorXd* EndpointPositions);

	Eigen::RowVectorXd CalculateConstantLineParameters(Eigen::RowVectorXd* EndpointPositions);

	Eigen::VectorXd GenerateSpacingRatioIntervals(unsigned int div, float SpacingRatio);

	bool ComputingNSLines(int i);

	void CreateNSLineNodesFinite(Eigen::RowVectorXd* EndpointPositions, Eigen::RowVectorXd* ConstantLineParameters, Eigen::VectorXd* SpacingRatioIntervals, Eigen::MatrixXd* BoundaryNS, int i);

	void CreateEWLineNodesFinite(Eigen::RowVectorXd* EndpointPositions, Eigen::RowVectorXd* ConstantLineParameters, Eigen::VectorXd* SpacingRatioIntervals, Eigen::MatrixXd* BoundaryEW, int i);

	void CreateNSLineNodesInfinite(Eigen::RowVectorXd* EndpointPositions, Eigen::VectorXd* SpacingRatioIntervals, Eigen::MatrixXd* BoundaryNS, int i);

	void CreateEWLineNodesInfinite(Eigen::RowVectorXd* EndpointPositions, Eigen::VectorXd* SpacingRatioIntervals, Eigen::MatrixXd* BoundaryEW, int i);

	// Resize the nodal matrices of an area
	std::vector<Eigen::MatrixXd> InitialiseAreaNodes(unsigned int id);

	void InputAreaBoundaryNodes(std::vector<Eigen::MatrixXd>* AreaLineNodes, std::vector<Eigen::MatrixXd>* NodePositions);

	// Constructs initial geometry matrix using Transfinite Interpolation
	void InitialiseAreaGeometry(std::vector<Eigen::MatrixXd>* NodePositions);

	void CreateAreaNodes(unsigned int id, std::vector<Eigen::MatrixXd>* NodePositions);

	void CreateNode(double x, double y);

	unsigned int AssignNodeID();

	void UpdateNextNewNodeID();

	void AssignNodePointer(unsigned int id);

	void UpdateNodeCount();

	void PrintNodes();

	void PrintNodes(unsigned int id);

	void InitialiseElements();

	unsigned int ReturnStrainMatrixSize();

	void CreateAreaElements(unsigned int id);

	bool ElementsInitialised();

	void CreateElement(Eigen::RowVectorXi* ElementNodes);

	unsigned int AssignElementID();

	void UpdateNextNewElementID();

	void AssignElementPointer(unsigned int id);

	void UpdateElementCount();

	void CheckMeshQuality();

	void ReturnElementLineVectors(Eigen::MatrixXd* LineVectors, unsigned int id);

	void ReturnElementLineLengths(Eigen::MatrixXd* LineVectors, Eigen::VectorXd* LineLengths);

	void ReturnElementArea(Eigen::VectorXd* LineLengths, Eigen::VectorXd* LineAngles, Eigen::VectorXd* Areas, unsigned int id);

	void CheckAspectRatio(Eigen::VectorXd* AspectRatios, Eigen::VectorXd* LineLengths, unsigned int id);

	void ReturnLineAngles(Eigen::MatrixXd* LineVectors, Eigen::VectorXd* LineLengths, Eigen::VectorXd* LineAngles);

	void CheckSkewness(Eigen::VectorXd* Skewnesses, Eigen::VectorXd* LineLengths, Eigen::VectorXd* Areas, unsigned int id);

	void PrintElements();

	void PrintElements(unsigned int id);
}
