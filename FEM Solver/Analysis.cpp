#include "Geometry.h"
#include "Analysis.h"

extern unsigned int Case;
extern unsigned int ElementType;
extern std::deque<Keypoint> KeypointsArray;
extern std::deque<Line> LinesArray;
extern std::deque<Area> AreasArray;
extern std::deque<Node> NodesArray;
extern std::deque<Element> ElementsArray;
extern std::vector<Keypoint*> KeypointPointers;
extern std::vector<Line*> LinePointers;
extern std::vector<Area*> AreaPointers;
extern std::vector<Node*> NodePointers;
extern std::vector<Element*> ElementPointers;
extern unsigned int NodeCount;
extern unsigned int ElementCount;

Material AnalysisMaterial = Material::Material();

Eigen::VectorXd F;
Eigen::VectorXd u;
Eigen::MatrixXd K;

void InitialiseFMatrix() {
	F = Eigen::VectorXd::Zero(2 * NodeCount);
}

void InitialiseuMatrix() {
	u = Eigen::VectorXd::Zero(2 * NodeCount);
}

void InitialiseKMatrix() {
	K = Eigen::MatrixXd::Zero(2 * NodeCount, 2 * NodeCount);
}

void ApplyConstraintsToLine(unsigned int id, bool ConstrainX, bool ConstrainY) {
	for (int i = 0; i < (*LinePointers[id]).Nodes.size(); i++) {
		(*NodePointers[(*LinePointers[id]).Nodes[i]]).Constraints(0) = ReturnConstraintValue(ConstrainX);
		(*NodePointers[(*LinePointers[id]).Nodes[i]]).Constraints(1) = ReturnConstraintValue(ConstrainY);
	}
}

void ApplyConstraintsToNode(unsigned int id, bool ConstrainX, bool ConstrainY) {
	for (int i = 0; i < (*LinePointers[id]).Nodes.size(); i++) {
		(*NodePointers[id]).Constraints(0) = ReturnConstraintValue(ConstrainX);
		(*NodePointers[id]).Constraints(1) = ReturnConstraintValue(ConstrainY);
	}
}

unsigned short int ReturnConstraintValue(bool Constrain) {
	if (Constrain) {
		return 1; // Might switch around
	}
	return 0;
}

void PrintNodeConstraints() {
	std::cout << "NODE CONSTRAINTS" << std::endl;
	for (unsigned int i = 0; i < NodesArray.size(); i++) {
		std::cout << "Node " << (*NodePointers[i]).ID << ":" << std::endl;
		std::cout << NodesArray[(*NodePointers[i]).ID].Constraints(0) << ", " << NodesArray[(*NodePointers[i]).ID].Constraints(1) << std::endl;
	}
	std::cout << " " << std::endl;
}

void PrintNodeConstraints(unsigned int id) {
	std::cout << "NODE CONSTRAINTS" << std::endl;
	std::cout << "Node " << id << ":" << std::endl;
	std::cout << (*NodePointers[id]).Constraints(0) << ", " << (*NodePointers[id]).Constraints(1) << std::endl;
	std::cout << " " << std::endl;
}

void ApplyXForcesToLine(unsigned int id, double Force) {
	double NodeForce = Force / (*LinePointers[id]).Nodes.size();
	for (int i = 0; i < (*LinePointers[id]).Nodes.size(); i++) {
		(*NodePointers[(*LinePointers[id]).Nodes[i]]).Forces(0) += NodeForce;
	}
}

void DeleteXForcesOnLine(unsigned int id) {
	for (int i = 0; i < (*LinePointers[id]).Nodes.size(); i++) {
		(*NodePointers[(*LinePointers[id]).Nodes[i]]).Forces(0) = 0.0;
	}
}

void ApplyXForcesToNode(unsigned int id, double Force) {
	(*NodePointers[id]).Forces(0) += Force;
}

void DeleteXForcesOnNode(unsigned int id) {
	(*NodePointers[id]).Forces(0) = 0.0;
}

void ApplyYForcesToLine(unsigned int id, double Force) {
	double NodeForce = Force / (*LinePointers[id]).Nodes.size();
	for (int i = 0; i < (*LinePointers[id]).Nodes.size(); i++) {
		(*NodePointers[(*LinePointers[id]).Nodes[i]]).Forces(1) += NodeForce;
	}
}

void DeleteYForcesOnLine(unsigned int id) {
	for (int i = 0; i < (*LinePointers[id]).Nodes.size(); i++) {
		(*NodePointers[(*LinePointers[id]).Nodes[i]]).Forces(1) = 0.0;
	}
}

void ApplyYForcesToNode(unsigned int id, double Force) {
	(*NodePointers[id]).Forces(1) += Force;
}

void DeleteYForcesOnNode(unsigned int id) {
	(*NodePointers[id]).Forces(1) = 0.0;
}

void PrintNodeForces() {
	std::cout << "NODE FORCES" << std::endl;
	for (unsigned int i = 0; i < NodesArray.size(); i++) {
		std::cout << "Node " << NodesArray[i].ID << ":" << std::endl;
		std::cout << NodesArray[i].Forces(0) << ", " << NodesArray[i].Forces(1) << std::endl;
	}
	std::cout << " " << std::endl;
}

void PrintNodeForces(unsigned int id) {
	std::cout << "NODE FORCES" << std::endl;
	std::cout << "Node " << id << ":" << std::endl;
	std::cout << (*NodePointers[id]).Forces(0) << ", " << (*NodePointers[id]).Forces(1) << std::endl;
	std::cout << " " << std::endl;
}

void SetMaterial(double E, double nu, double thickness) { // Material properties dependent on case! Should work for isotropic
	AnalysisMaterial.E = E;
	AnalysisMaterial.nu = nu;
	AnalysisMaterial.thickness = thickness;
}

void SolveSystem() {
	std::cout << "Node count: " << NodeCount << std::endl;
	std::cout << "Element count: " << ElementCount << std::endl;
	TransferForcesToFMatrix();
	GenerateKMatrix();

	Eigen::VectorXd uReduced;
	std::vector<int> KeptIndices;
	std::tie(uReduced, KeptIndices) = GenerateReducedSystem();
	std::cout << "Reduced node count: " << uReduced.size() / 2.0 << std::endl;
	std::cout << " " << std::endl;

	uReduced = K.lu().solve(F);
	PopulateDisplacements(&uReduced, &KeptIndices);
}

void TransferForcesToFMatrix() {
	InitialiseFMatrix();
	for (int i = 0; i < NodeCount; i++) {
		F(2 * i) = NodesArray[i].Forces(0); // Does using Array matter? Is connectivity covered enough?
		F(2 * i + 1) = NodesArray[i].Forces(1);
	}
}

void GenerateKMatrix() {
	InitialiseKMatrix();
	InitialiseuMatrix();
	Eigen::MatrixXd D = GenerateDMatrix();
	Eigen::MatrixXd k = Eigen::MatrixXd::Zero(2 * ElementType, 2 * ElementType);
	K = Eigen::MatrixXd::Zero(2 * NodeCount, 2 * NodeCount);

	for (int i = 0; i < ElementCount; i++) {
		GeneratekMatrix(&D, &k, i);
		PopulateKMatrix(i, &k);
		k.setZero();
	}

	K *= AnalysisMaterial.thickness;
}

// Assumed only isotropic for now
Eigen::MatrixXd GenerateDMatrix() {
	Eigen::MatrixXd D;
	switch (Case) {
	case 0:
		D = GenerateDMatrixPlaneStress();
		return D;
	case 1:
		D = GenerateDMatrixPlaneStrain();
		return D;
	case 2:
		D = GenerateDMatrixIsotropic();
		return D;
	default:
		D = Eigen::MatrixXd::Zero(1, 1); // Invalid case specified, should be error and returns nothing but can be overwritten
		return D;
	}
}

Eigen::MatrixXd GenerateDMatrixPlaneStress() {
	unsigned int DMatrixSize = Geometry::ReturnStrainMatrixSize();
	Eigen::MatrixXd D(DMatrixSize, DMatrixSize);
	double E = AnalysisMaterial.E;
	double nu = AnalysisMaterial.nu;
	D << 1.0, nu, 0.0,
		nu, 1.0, 0.0,
		0.0, 0.0, 0.5 * (1.0 - nu);
	D *= E / (1 - pow(nu, 2));
	std::cout << "Plane stress isotropic" << std::endl;
	return D;
}

Eigen::MatrixXd GenerateDMatrixPlaneStrain() {
	unsigned int DMatrixSize = Geometry::ReturnStrainMatrixSize();
	Eigen::MatrixXd D(DMatrixSize, DMatrixSize);
	double E = AnalysisMaterial.E;
	double nu = AnalysisMaterial.nu;
	D << 1.0 - nu, nu, 0.0,
		nu, 1.0 - nu, 0.0,
		0.0, 0.0, 0.5 * (1.0 - 2.0 * nu);
	D *= E / ((1.0 + nu) * (1.0 - 2.0 * nu));
	std::cout << "Plane strain isotropic" << std::endl;
	return D;
}

Eigen::MatrixXd GenerateDMatrixIsotropic() {
	unsigned int DMatrixSize = Geometry::ReturnStrainMatrixSize();
	Eigen::MatrixXd D(DMatrixSize, DMatrixSize);
	double E = AnalysisMaterial.E;
	double nu = AnalysisMaterial.nu;
	// Lamé constants
	double mu = E / (2 * (1 + nu));
	double lamda = E * nu / ((1 + nu) * (1 - 2 * nu));
	D << 2.0 * mu + lamda, lamda, 0.0, lamda, 0.0, 0.0,
		lamda, 2.0 * mu + lamda, 0.0, lamda, 0.0, 0.0,
		0.0, 0.0, mu, 0.0, 0.0, 0.0,
		lamda, lamda, 0.0, 2.0 * mu + lamda, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, mu, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, mu;
	std::cout << "3D Isotropic" << std::endl;
	return D;
}

void GeneratekMatrix(Eigen::MatrixXd* D, Eigen::MatrixXd* k, int id) {
	Eigen::MatrixXd IntegrationPoints = ReturnIntegrationPoints(0);
	for (int j = 0; j < ElementType; j++) {
		Eigen::MatrixXd B;
		double Jdet;
		GenerateBMatrix(id, j, &IntegrationPoints, &B, &Jdet);
		Eigen::MatrixXd k_temp = Jdet * B.transpose() * (*D) * B; // Preallocate k_temp size?
		*k += k_temp;
	}
}

void GenerateBMatrix(int id, int j, Eigen::MatrixXd* IntegrationPoints, Eigen::MatrixXd* B, double* Jdet) { // id, j not good
	if (Case == 0 || Case == 1) { // Also waiting for classes
		Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2, 2);
		Eigen::MatrixXd dN_st = Eigen::MatrixXd::Zero(2, ElementType);
		Eigen::MatrixXd dN_xy = Eigen::MatrixXd::Zero(2, ElementType);
		*B = Eigen::MatrixXd::Zero(3, 2 * ElementType);

		dN_st << (*IntegrationPoints)(j, 1) - 1.0, 1.0 - (*IntegrationPoints)(j, 1), -(*IntegrationPoints)(j, 1) - 1.0, (*IntegrationPoints)(j, 1) + 1.0,
			(*IntegrationPoints)(j, 0) - 1.0, -(*IntegrationPoints)(j, 0) - 1.0, 1.0 - (*IntegrationPoints)(j, 0), (*IntegrationPoints)(j, 0) + 1.0;
		dN_st /= 4.0;

		Eigen::VectorXd x(ElementType), y(ElementType);
		for (int i = 0; i < ElementType; i++) {
			x(i) = NodesArray[ElementsArray[id].Nodes(i)].Position(0);
			y(i) = NodesArray[ElementsArray[id].Nodes(i)].Position(1);
		}

		J.row(0) = dN_st * x;
		J.row(1) = dN_st * y;

		for (int i = 0; i < ElementType; i++) {
			dN_xy.col(i) = J.inverse() * dN_st.col(i);
			(*B)(0, 2 * i) = (*B)(2, 2 * i + 1) = dN_xy(0, i);
			(*B)(1, 2 * i + 1) = (*B)(2, 2 * i) = dN_xy(1, i);
		}
		*Jdet = J.determinant();
	}
}

void GenerateShapeFunctions(Eigen::MatrixXd* IntegrationPoints, Eigen::MatrixXd* ShapeFunctions) {
	for (int i = 0; i < ElementType; i++) {
		(*ShapeFunctions).row(i) << 1 - (*IntegrationPoints)(i, 0) - (*IntegrationPoints)(i, 1) - (*IntegrationPoints)(i, 0) * (*IntegrationPoints)(i, 1),
			1 + (*IntegrationPoints)(i, 0) - (*IntegrationPoints)(i, 1) - (*IntegrationPoints)(i, 0) * (*IntegrationPoints)(i, 1),
			1 - (*IntegrationPoints)(i, 0) + (*IntegrationPoints)(i, 1) - (*IntegrationPoints)(i, 0) * (*IntegrationPoints)(i, 1),
			1 + (*IntegrationPoints)(i, 0) + (*IntegrationPoints)(i, 1) + (*IntegrationPoints)(i, 0) * (*IntegrationPoints)(i, 1);
	}
}

// Returns integration points matrix depending on element type and evaluation points (Gauss, extrapolated or nodes)
Eigen::MatrixXd ReturnIntegrationPoints(int Locations) {
	switch (ElementType) { // Horrible, will stay until dif element types added, same function names under different classes
	case 4:
		Eigen::MatrixXd IntegrationPoints(4, 2);
		switch (Locations) {
		case 0:
			ReturnGaussIntegrationPoints(&IntegrationPoints);
			return IntegrationPoints;
		case 1:
			ReturnExtrapolationIntegrationPoints(&IntegrationPoints);
			return IntegrationPoints;
		case 2:
			ReturnNodeIntegrationPoints(&IntegrationPoints);
			return IntegrationPoints;

		}
	}
}

void ReturnGaussIntegrationPoints(Eigen::MatrixXd* IntegrationPoints) { // 4 Points only!
	*IntegrationPoints << -sqrt(3), -sqrt(3),
		sqrt(3), -sqrt(3),
		-sqrt(3), sqrt(3),
		sqrt(3), sqrt(3);
	*IntegrationPoints /= 3.0;
}

void ReturnExtrapolationIntegrationPoints(Eigen::MatrixXd* IntegrationPoints) {
	*IntegrationPoints << -sqrt(3), -sqrt(3),
		sqrt(3), -sqrt(3),
		-sqrt(3), sqrt(3),
		sqrt(3), sqrt(3);
}

void ReturnNodeIntegrationPoints(Eigen::MatrixXd* IntegrationPoints) { // 4 Points only!
	*IntegrationPoints << -1, -1,
		1, -1,
		-1, 1,
		1, 1;
}

void PopulateKMatrix(int id, Eigen::MatrixXd* k) {
	for (int i = 0; i < ElementType; i++) {
		for (int j = 0; j < ElementType; j++) {
			K.block(2 * ElementsArray[id].Nodes(i), 2 * ElementsArray[id].Nodes(j), 2, 2) += (*k).block(2 * i, 2 * j, 2, 2);
		}
	}
}

std::tuple<Eigen::VectorXd, std::vector<int>> GenerateReducedSystem() {
	std::vector<int> DeletedIndices = TransferConstraints();
	int FinalMatrixSize;
	ReduceMatrices(&DeletedIndices, &FinalMatrixSize);
	Eigen::VectorXd uReduced(FinalMatrixSize);
	std::vector<int> KeptIndices = GenerateKeptIndices(&DeletedIndices);

	return std::make_tuple(uReduced, KeptIndices);
}

std::vector<int> TransferConstraints() {
	int ConstraintsCount = SumConstraints();
	std::vector<int> DeletedIndices;
	DeletedIndices.reserve(ConstraintsCount);
	for (int i = 0; i < NodeCount; i++) {
		if (NodesArray[i].Constraints(0) == 1) {
			K.row(2 * i).setZero();
			K.col(2 * i).setZero();
			F(2 * i) *= 0;
			DeletedIndices.emplace_back(2 * i);
		}
		if (NodesArray[i].Constraints(1) == 1) {
			K.row(2 * i + 1).setZero();
			K.col(2 * i + 1).setZero();
			F(2 * i + 1) *= 0;
			DeletedIndices.emplace_back(2 * i + 1);
		}
	}
	return DeletedIndices;
}

int SumConstraints() {
	int ConstraintsCount = 0;
	for (int i = 0; i < NodeCount; i++) {
		if (NodesArray[i].Constraints(0) == 1) {
			ConstraintsCount++;
		}
		if (NodesArray[i].Constraints(1) == 1) {
			ConstraintsCount++;
		}
	}
	return ConstraintsCount;
}

void ReduceMatrices(std::vector<int>* DeletedIndices, int* FinalMatrixSize) {
	int rows = K.rows();
	int columns = K.cols();
	int Constraints = SumConstraints();
	int DeletedIndicesSize = (*DeletedIndices).size(); // Quicker than constantly calculating size?

	for (int i = 0; i < DeletedIndicesSize; i++) {
		K.block((*DeletedIndices)[i] - i, 0, rows - 1 - ((*DeletedIndices)[i] - i), columns) = K.block((*DeletedIndices)[i] - i + 1, 0, rows - 1 - ((*DeletedIndices)[i] - i), columns);
		F.block((*DeletedIndices)[i] - i, 0, rows - 1 - ((*DeletedIndices)[i] - i), 1) = F.block((*DeletedIndices)[i] - i + 1, 0, rows - 1 - ((*DeletedIndices)[i] - i), 1);
		K.conservativeResize(rows - 1, columns);
		F.conservativeResize(rows - 1, 1);
		rows = K.rows();
	}

	for (int i = 0; i < DeletedIndicesSize; i++) {
		K.block(0, (*DeletedIndices)[i] - i, rows, columns - 1 - ((*DeletedIndices)[i] - i)) = K.block(0, (*DeletedIndices)[i] - i + 1, rows, columns - 1 - ((*DeletedIndices)[i] - i));
		K.conservativeResize(rows, columns - 1);
		columns = K.cols();
	}
	*FinalMatrixSize = rows;
}

std::vector<int> GenerateKeptIndices(std::vector<int>* DeletedIndices) { // Might not be necessary for populating displacements since they could be init to 0
	std::vector<int> KeptIndices;
	KeptIndices.reserve(2 * NodeCount);
	for (int i = 0; i < 2 * NodeCount; i++) {
		KeptIndices.emplace_back(i);
	}

	int DeletedIndicesSize = (*DeletedIndices).size();
	for (int i = 0; i < DeletedIndicesSize; i++) {
		KeptIndices.erase(KeptIndices.begin() + (*DeletedIndices)[i] - i);
	}
	return KeptIndices;
}

void PopulateDisplacements(Eigen::VectorXd* uReduced, std::vector<int>* KeptIndices) {
	int KeptIndicesSize = (*KeptIndices).size();
	for (int i = 0; i < KeptIndicesSize; i++) {
		u((*KeptIndices)[i]) = (*uReduced)(i);
		(*NodePointers[(*KeptIndices)[i] / 2]).Displacements((*KeptIndices)[i] % 2) = (*uReduced)(i);
	}
}

void PrintNodeDisplacements() {
	std::cout << "NODE DISPLACEMENTS" << std::endl;
	for (unsigned int i = 0; i < NodesArray.size(); i++) {
		std::cout << "Node " << NodesArray[i].ID << ":" << std::endl;
		std::cout << NodesArray[i].Displacements(0) << ", " << NodesArray[i].Displacements(1) << std::endl;
	}
	std::cout << " " << std::endl;
}

void PrintNodeDisplacements(unsigned int id) {
	std::cout << "NODE DISPLACEMENTS" << std::endl;
	std::cout << "Node " << id << ":" << std::endl;
	std::cout << (*NodePointers[id]).Displacements(0) << ", " << (*NodePointers[id]).Displacements(1) << std::endl;
	std::cout << " " << std::endl;
}

void PrintElementDisplacements() {
	std::cout << "ELEMENT DISPLACEMENTS" << std::endl;
	for (unsigned int i = 0; i < ElementsArray.size(); i++) {
		std::cout << "Element " << ElementsArray[i].ID << ":" << std::endl;
		for (int j = 0; j < ElementType; j++) {
			std::cout << "Node " << ElementsArray[i].Nodes(j) << ":" << std::endl;
			std::cout << (*NodePointers[ElementsArray[i].Nodes(j)]).Displacements(0) <<
				", " << (*NodePointers[ElementsArray[i].Nodes(j)]).Displacements(1) << std::endl;
		}
	}
	std::cout << " " << std::endl;
}

void PrintElementDisplacements(unsigned int id) {
	std::cout << "ELEMENT DISPLACEMENTS" << std::endl;
	std::cout << "Element " << id << ":" << std::endl;
	for (int j = 0; j < ElementType; j++) {
		std::cout << "Node " << (*ElementPointers[id]).Nodes(j) << ":" << std::endl;
		std::cout << (*NodePointers[(*ElementPointers[id]).Nodes(j)]).Displacements(0) << ", " << (*NodePointers[(*ElementPointers[id]).Nodes(j)]).Displacements(1) << std::endl;
	}
	std::cout << "" << std::endl;
}

void GenerateNodeStrainsAndStresses() {
	Eigen::MatrixXd NodeIntegrationPoints = ReturnIntegrationPoints(2);
	Eigen::MatrixXd D = GenerateDMatrix();
	Eigen::VectorXd Displacements(2 * ElementType);
	int MatrixSize = Geometry::ReturnStrainMatrixSize();
	Eigen::VectorXd Strains(MatrixSize);
	Eigen::VectorXd Stresses(MatrixSize);//

	for (int i = 0; i < ElementCount; i++) {
		GetElementDisplacements(i, &Displacements);
		GenerateElementNodeStrainsAndStresses(i, &NodeIntegrationPoints, &D, &Displacements, MatrixSize, &Strains, &Stresses);
	}
}

void GenerateExtrapolatedStrainsAndStresses() {
	Eigen::MatrixXd GaussIntegrationPoints = ReturnIntegrationPoints(0);
	Eigen::MatrixXd ExtrapolationIntegrationPoints = ReturnIntegrationPoints(1);
	Eigen::MatrixXd D = GenerateDMatrix();
	Eigen::VectorXd Displacements(2 * ElementType);
	int MatrixSize = Geometry::ReturnStrainMatrixSize();
	Eigen::MatrixXd GaussStrains(ElementType, MatrixSize);
	Eigen::MatrixXd GaussStresses(ElementType, MatrixSize);
	Eigen::MatrixXd ExtrapolatedStrains(ElementType, MatrixSize);
	Eigen::MatrixXd ExtrapolatedStresses(ElementType, MatrixSize);
	Eigen::MatrixXd ShapeFunctions(ElementType, ElementType);
	GenerateShapeFunctions(&ExtrapolationIntegrationPoints, &ShapeFunctions);

	for (int i = 0; i < ElementCount; i++) {
		GetElementDisplacements(i, &Displacements);
		GenerateElementExtrapolatedStrainsAndStresses(i, &GaussIntegrationPoints, &ShapeFunctions, &D, &Displacements, MatrixSize, &GaussStrains, &GaussStresses, &ExtrapolatedStrains, &ExtrapolatedStresses);
	}
}

void GetElementDisplacements(unsigned int id, Eigen::VectorXd* Displacements) {
	for (int i = 0; i < ElementType; i++) {
		(*Displacements)(2 * i) = (*NodePointers[(*ElementPointers[id]).Nodes(i)]).Displacements(0);
		(*Displacements)(2 * i + 1) = (*NodePointers[(*ElementPointers[id]).Nodes(i)]).Displacements(1);
	}
}

void GenerateElementNodeStrainsAndStresses(unsigned int id, Eigen::MatrixXd* NodeIntegrationPoints, Eigen::MatrixXd* D, Eigen::VectorXd* Displacements, int MatrixSize, Eigen::VectorXd* Strains, Eigen::VectorXd* Stresses) {
	for (int j = 0; j < ElementType; j++) {
		Eigen::MatrixXd B;
		double Jdet;
		GenerateBMatrix(id, j, NodeIntegrationPoints, &B, &Jdet);
		*Strains = B * (*Displacements); // ex, gamxy, sigx, tauxy correct BUT ey, sigy incorrect
		*Stresses = *D * (*Strains);

		(*ElementPointers[id]).Strains.row(j) = *Strains;
		(*ElementPointers[id]).Stresses.block(j, 0, 1, MatrixSize) = (*Stresses).transpose(); // Might be able to change to avoid needing to transpose, or not since stresses must be cvec
		(*ElementPointers[id]).Stresses(j, MatrixSize) = pow(pow((*ElementPointers[id]).Stresses(j, 0), 2)
			- (*ElementPointers[id]).Stresses(j, 0) * (*ElementPointers[id]).Stresses(j, 1)
			+ pow((*ElementPointers[id]).Stresses(j, 1), 2) + 3 * pow((*ElementPointers[id]).Stresses(j, 2), 2), 0.5);
	}
}

void GenerateElementExtrapolatedStrainsAndStresses(unsigned int id, Eigen::MatrixXd* GaussIntegrationPoints, Eigen::MatrixXd* ShapeFunctions, Eigen::MatrixXd* D, Eigen::VectorXd* Displacements, int MatrixSize, Eigen::MatrixXd* GaussStrains, Eigen::MatrixXd* GaussStresses, Eigen::MatrixXd* ExtrapolatedStrains, Eigen::MatrixXd* ExtrapolatedStresses) {
	for (int j = 0; j < ElementType; j++) {
		Eigen::MatrixXd B;
		double Jdet;
		GenerateBMatrix(id, j, GaussIntegrationPoints, &B, &Jdet);
		(*GaussStrains).row(j) = B * (*Displacements); // ex, gamxy, sigx, tauxy correct BUT ey, sigy incorrect
		(*GaussStresses).row(j) = *D * ((*GaussStrains).row(j)).transpose();
	}

	for (int j = 0; j < ElementType; j++) {
		(*ExtrapolatedStrains)(j, 0) = (*ShapeFunctions).row(j) * (*GaussStrains).col(0);
		(*ExtrapolatedStrains)(j, 1) = (*ShapeFunctions).row(j) * (*GaussStrains).col(1);
		(*ExtrapolatedStrains)(j, 2) = (*ShapeFunctions).row(j) * (*GaussStrains).col(2); // Depends on case
		(*ExtrapolatedStresses)(j, 0) = (*ShapeFunctions).row(j) * (*GaussStresses).col(0);
		(*ExtrapolatedStresses)(j, 1) = (*ShapeFunctions).row(j) * (*GaussStresses).col(1);
		(*ExtrapolatedStresses)(j, 2) = (*ShapeFunctions).row(j) * (*GaussStresses).col(2);

		(*ElementPointers[id]).Strains.row(j) = (*ExtrapolatedStrains).row(j);
		(*ElementPointers[id]).Stresses.block(j, 0, 1, MatrixSize) = (*ExtrapolatedStresses).row(j);
		(*ElementPointers[id]).Stresses(j, MatrixSize) = pow(pow((*ElementPointers[id]).Stresses(j, 0), 2)
			- (*ElementPointers[id]).Stresses(j, 0) * (*ElementPointers[id]).Stresses(j, 1)
			+ pow((*ElementPointers[id]).Stresses(j, 1), 2) + 3 * pow((*ElementPointers[id]).Stresses(j, 2), 2), 0.5);
	}
}

void PrintElementStrains() {
	std::cout << "ELEMENT STRAINS" << std::endl;
	for (unsigned int i = 0; i < ElementsArray.size(); i++) {
		std::cout << "Element " << ElementsArray[i].ID << ":" << std::endl;
		for (int j = 0; j < ElementType; j++) {
			std::cout << "Node " << ElementsArray[i].Nodes(j) << ":" << std::endl;
			std::cout << ElementsArray[i].Strains.row(j) << std::endl;
		}
	}
	std::cout << "" << std::endl;
}

void PrintElementStrains(unsigned int id) {
	std::cout << "ELEMENT STRAINS" << std::endl;
	std::cout << "Element " << id << ":" << std::endl;
	for (int i = 0; i < ElementType; i++) {
		std::cout << "Node " << (*ElementPointers[id]).Nodes(i) << ":" << std::endl;
		std::cout << (*ElementPointers[id]).Strains.row(i) << std::endl;
	}
	std::cout << "" << std::endl;
}

void PrintElementStresses() {
	std::cout << "ELEMENT STRESSES" << std::endl;
	for (unsigned int i = 0; i < ElementsArray.size(); i++) {
		std::cout << "Element " << ElementsArray[i].ID << ":" << std::endl;
		for (int j = 0; j < ElementType; j++) {
			std::cout << "Node " << ElementsArray[i].Nodes(j) << ":" << std::endl;
			std::cout << ElementsArray[i].Stresses.row(j) << std::endl;
		}
	}
	std::cout << "" << std::endl;
}

void PrintElementStresses(unsigned int id) {
	std::cout << "ELEMENT STRESSES" << std::endl;
	std::cout << "Element " << id << ":" << std::endl;
	for (int i = 0; i < ElementType; i++) {
		std::cout << "Node " << (*ElementPointers[id]).Nodes(i) << ":" << std::endl;
		std::cout << (*ElementPointers[id]).Stresses.row(i) << std::endl;
	}
	std::cout << "" << std::endl;
}
