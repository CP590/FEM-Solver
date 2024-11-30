#include "Geometry.h"

extern unsigned int ElementType;
extern unsigned int Case;
unsigned int NodeCount;
unsigned int ElementCount;
unsigned int NextNewKeypointID = 0;
unsigned int NextNewLineID = 0;
unsigned int NextNewAreaID = 0;
unsigned int NextNewNodeID = 0;
unsigned int NextNewElementID = 0;

std::deque<Keypoint> Geometry::CreateKeypointsArray() {
    std::deque<Keypoint> KeypointsArray;
    return KeypointsArray;
}
std::deque<Line> Geometry::CreateLinesArray() {
    std::deque<Line> LinesArray;
    return LinesArray;
}
std::deque<Area> Geometry::CreateAreasArray() {
    std::deque<Area> AreasArray;
    return AreasArray;
}
std::deque<Node> Geometry::CreateNodesArray() {
    std::deque<Node> NodesArray;
    return NodesArray;
}
std::deque<Element> Geometry::CreateElementsArray() {
    std::deque<Element> ElementsArray;
    return ElementsArray;
}

std::vector<unsigned int> Geometry::CreateDeletedKeypointIDs() {
    std::vector<unsigned int> DeletedKeypointIDs;
    DeletedKeypointIDs.reserve(10);
    return DeletedKeypointIDs;
}
std::vector<unsigned int> Geometry::CreateDeletedLineIDs() {
    std::vector<unsigned int> DeletedLineIDs;
    DeletedLineIDs.reserve(10);
    return DeletedLineIDs;
}
std::vector<unsigned int> Geometry::CreateDeletedAreaIDs() {
    std::vector<unsigned int> DeletedAreaIDs;
    DeletedAreaIDs.reserve(10);
    return DeletedAreaIDs;
}
std::vector<unsigned int> Geometry::CreateDeletedNodeIDs() {
    std::vector<unsigned int> DeletedNodeIDs;
    DeletedNodeIDs.reserve(10);
    return DeletedNodeIDs;
}
std::vector<unsigned int> Geometry::CreateDeletedElementIDs() {
    std::vector<unsigned int> DeletedElementIDs;
    DeletedElementIDs.reserve(10);
    return DeletedElementIDs;
}

std::vector<Keypoint*> Geometry::CreateKeypointPointers() {
    std::vector<Keypoint*> KeypointPointers;
    KeypointPointers.reserve(10);
    return KeypointPointers;
}
std::vector<Line*> Geometry::CreateLinePointers() {
    std::vector<Line*> LinePointers;
    LinePointers.reserve(10);
    return LinePointers;
}
std::vector<Area*> Geometry::CreateAreaPointers() {
    std::vector<Area*> AreaPointers;
    AreaPointers.reserve(10);
    return AreaPointers;
}
std::vector<Node*> Geometry::CreateNodePointers() {
    std::vector<Node*> NodePointers;
    NodePointers.reserve(10);
    return NodePointers;
}
std::vector<Element*> Geometry::CreateElementPointers() {
    std::vector<Element*> ElementPointers;
    ElementPointers.reserve(10);
    return ElementPointers;
}

std::deque<Keypoint> KeypointsArray = Geometry::CreateKeypointsArray();
std::deque<Line> LinesArray = Geometry::CreateLinesArray();
std::deque<Area> AreasArray = Geometry::CreateAreasArray();
std::deque<Node> NodesArray = Geometry::CreateNodesArray();
std::deque<Element> ElementsArray = Geometry::CreateElementsArray();

std::vector<unsigned int> DeletedKeypointIDs = Geometry::CreateDeletedKeypointIDs();
std::vector<unsigned int> DeletedLineIDs = Geometry::CreateDeletedLineIDs();
std::vector<unsigned int> DeletedAreaIDs = Geometry::CreateDeletedAreaIDs();
std::vector<unsigned int> DeletedNodeIDs = Geometry::CreateDeletedNodeIDs();
std::vector<unsigned int> DeletedElementIDs = Geometry::CreateDeletedElementIDs();

std::vector<Keypoint*> KeypointPointers = Geometry::CreateKeypointPointers();
std::vector<Line*> LinePointers = Geometry::CreateLinePointers();
std::vector<Area*> AreaPointers = Geometry::CreateAreaPointers();
std::vector<Node*> NodePointers = Geometry::CreateNodePointers();
std::vector<Element*> ElementPointers = Geometry::CreateElementPointers();

Keypoint KP;
Line L;
Area A;
Node N;
Element E;

void Geometry::CreateKeypoint(double x, double y) {
    if (!Geometry::KeypointExistsInArray(x, y)) {
        KP.ID = Geometry::AssignKeypointID();
        KP.Position << x, y;
        KeypointsArray.emplace_back(KP);
        Geometry::AssignKeypointPointer(KP.ID);
    }
}

bool Geometry::KeypointExistsInArray(double x, double y) {
    for (int i = 0; i < KeypointsArray.size(); i++) {
        if (Geometry::KeypointExists(i, x, y)) {
            return true;
        }
    }
    return false;
}

bool Geometry::KeypointExists(int id, double x, double y) {
    if (KeypointsArray[id].Position(0) == x && KeypointsArray[id].Position(1) == y) {
        return true;
    }
    return false;
}

unsigned int Geometry::AssignKeypointID() {
    unsigned int KeypointID;
    if (DeletedKeypointIDs.size() > 0) {
        KeypointID = DeletedKeypointIDs[0];
        DeletedKeypointIDs.erase(DeletedKeypointIDs.begin());
        return KeypointID;
    }
    else {
        KeypointID = NextNewKeypointID;
        Geometry::UpdateNextNewKeypointID();
        return KeypointID;
    }
}

void Geometry::UpdateNextNewKeypointID() {
    if (KeypointsArray.size() != 0) {
        if (DeletedKeypointIDs.size() == 0) {
            NextNewKeypointID++;
        }
    }
    else {
        NextNewKeypointID++;
    }
}

void Geometry::AssignKeypointPointer(unsigned int id) {
    if (id == NextNewKeypointID - 1 && KeypointPointers.size() == id) {
        KeypointPointers.emplace_back(&(KeypointsArray.back()));
    }
    else if (id == NextNewKeypointID - 1 && KeypointPointers.back() == nullptr) {
        KeypointPointers[id] = &(KeypointsArray.back());
    }
    else {
        KeypointPointers[id] = &(KeypointsArray.back());
    }
}

void Geometry::DeleteKeypoint(unsigned int id) { // Should check if keypoint with ID exists
    Geometry::AddDeletedKeypointID(id); // Sort?
    std::deque<Keypoint>::iterator Location = std::find_if(KeypointsArray.begin(), KeypointsArray.end(), [id](Keypoint& i) {return i.ID == id; });
    if (Location != KeypointsArray.end()) {
        int DistanceToStart = Location - KeypointsArray.begin();
        int DistanceToEnd = KeypointsArray.end() - Location - 1;
        KeypointsArray.erase(Location);
        Geometry::AdjustKeypointPointers(id, DistanceToStart, DistanceToEnd);
    }
}

void Geometry::AddDeletedKeypointID(unsigned int id) {
    DeletedKeypointIDs.emplace_back(id);
}

void Geometry::AdjustKeypointPointers(unsigned int id, int DistanceToStart, int DistanceToEnd) {
    if (DistanceToStart != 0 && DistanceToEnd != 0) {
        int Size = KeypointPointers.size() - 1;
        if (DistanceToStart < DistanceToEnd) {
            //delete(KeypointPointers[0]); delete ONLY for heap allocated, will have to stay on stack
            KeypointPointers[0] = nullptr;
            for (int i = 0; i < DistanceToStart; i++) {
                KeypointPointers[i] = KeypointPointers[i + 1];
            }
            KeypointPointers[id] = nullptr;
        }
        else {
            //delete(KeypointPointers.back());
            KeypointPointers.back() = nullptr;
            for (int i = 0; i < DistanceToEnd; i++) {
                KeypointPointers[Size - i] = KeypointPointers[Size - (i + 1)];
            }
            KeypointPointers[id] = nullptr;
        }
    }
    else {
        KeypointPointers[id] = nullptr;
    }
}

// Print keypoints with positions
void Geometry::PrintKeypoints() {
    std::cout << "KEYPOINT POSITIONS" << std::endl;
    for (unsigned int i = 0; i < KeypointsArray.size(); i++) {
        std::cout << "Keypoint " << KeypointsArray[i].ID << ":" << std::endl; // Cannot ID nullptr! switch "i" to KeypointsArray[i].ID
        std::cout << KeypointsArray[i].Position(0) << ", " << KeypointsArray[i].Position(1) << std::endl;
    }
    std::cout << " " << std::endl;
}

void Geometry::PrintKeypoints(unsigned int id) {
    std::cout << "KEYPOINT POSITION" << std::endl;
    std::cout << "Keypoint " << id << ":" << std::endl;
    std::cout << (*KeypointPointers[id]).Position(0) << ", " << (*KeypointPointers[id]).Position(1) << std::endl;
    std::cout << " " << std::endl;
}

void Geometry::CreateLine(unsigned int Point0, unsigned int Point1) {
    if (!Geometry::LineExistsInArray(Point0, Point1)) {
        L.ID = Geometry::AssignLineID();
        L.Endpoints << (int)Point0, (int)Point1;
        LinesArray.emplace_back(L);
        Geometry::AssignLinePointer(L.ID);
    }
}

bool Geometry::LineExistsInArray(unsigned int Keypoint0, unsigned int Keypoint1) {
    for (unsigned int i = 0; i < LinesArray.size(); i++) {
        if (Geometry::LineExists(i, Keypoint0, Keypoint1)) {
            return true;
        }
    }
    return false;
}

bool Geometry::LineExists(int id, unsigned int Keypoint0, unsigned int Keypoint1) {
    if ((LinesArray[id].Endpoints(0) == Keypoint0 && LinesArray[id].Endpoints(1) == Keypoint1) ||
        (LinesArray[id].Endpoints(0) == Keypoint1 && LinesArray[id].Endpoints(1) == Keypoint0)) {
        return true;
    }
    return false;
}

unsigned int Geometry::AssignLineID() {
    unsigned int LineID;
    if (DeletedLineIDs.size() > 0) {
        LineID = DeletedLineIDs[0];
        DeletedLineIDs.erase(DeletedLineIDs.begin());
        return LineID;
    }
    else {
        LineID = NextNewLineID;
        Geometry::UpdateNextNewLineID();
        return LineID;
    }
}

void Geometry::UpdateNextNewLineID() {
    if (LinesArray.size() != 0) {
        if (DeletedLineIDs.size() == 0) {
            NextNewLineID++;
        }
    }
    else {
        NextNewLineID++;
    }
}

void Geometry::AssignLinePointer(unsigned int id) {
    if (id == NextNewLineID - 1) {
        LinePointers.emplace_back(&(LinesArray.back()));
    }
    else {
        LinePointers[id] = &(LinesArray.back());
    }
}

void Geometry::DeleteLine(unsigned int id) {
    Geometry::AddDeletedLineID(id);
    std::deque<Line>::iterator Location = std::find_if(LinesArray.begin(), LinesArray.end(), [id](Line& i) {return i.ID == id; });
    if (Location != LinesArray.end()) {
        int DistanceToStart = Location - LinesArray.begin();
        int DistanceToEnd = LinesArray.end() - Location - 1;
        LinesArray.erase(Location); // What if already part of an area?
        Geometry::AdjustLinePointers(id, DistanceToStart, DistanceToEnd);
    }
}

// Do these need their own function?
void Geometry::AddDeletedLineID(unsigned int id) {
    DeletedLineIDs.emplace_back(id);
}

void Geometry::AdjustLinePointers(unsigned int id, int DistanceToStart, int DistanceToEnd) {
    if (DistanceToStart != 0 && DistanceToEnd != 0) {
        int Size = LinePointers.size() - 1;
        if (DistanceToStart < DistanceToEnd) {
            delete(LinePointers[0]);
            for (int i = 0; i < DistanceToStart; i++) {
                LinePointers[i] = LinePointers[i + 1];
            }
            LinePointers[id] = nullptr;
        }
        else {
            delete(LinePointers.back());
            for (int i = 0; i < DistanceToEnd; i++) {
                LinePointers[Size - i] = LinePointers[Size - (i + 1)];
            }
            LinePointers[id] = nullptr;
        }
    }
    else {
        LinePointers[id] = nullptr;
    }
}

void Geometry::AssignLineDivision(unsigned int id, unsigned int div) {
    (*LinePointers[id]).Divisions = div;
    (*LinePointers[id]).Nodes.reserve(div);
}

void Geometry::AssignLineSpacingRatio(unsigned int id, float SpacingRatio) { // Could combine with AssignLineDivision
    (*LinePointers[id]).SpacingRatio = SpacingRatio;
}

void Geometry::PrintLines() {
    std::cout << "LINE KEYPOINTS" << std::endl;
    for (unsigned int i = 0; i < LinesArray.size(); i++) {
        std::cout << "Line " << LinesArray[i].ID << ":" << std::endl;
        std::cout << LinesArray[i].Endpoints(0) << ", " << LinesArray[i].Endpoints(1) << std::endl;
    }
    std::cout << " " << std::endl;
}

void Geometry::PrintLines(unsigned int id) {
    std::cout << "LINE KEYPOINTS" << std::endl;
    std::cout << "Line " << id << ":" << std::endl;
    std::cout << (*LinePointers[id]).Endpoints(0) << ", " << (*LinePointers[id]).Endpoints(1) << std::endl;
    std::cout << " " << std::endl;
}

void Geometry::PrintLineNodes() {
    std::cout << "LINE NODES" << std::endl;
    for (unsigned int i = 0; i < LinesArray.size(); i++) {
        std::cout << "Line " << (*LinePointers[i]).ID << ":" << std::endl;
        std::string NodeList;
        for (unsigned int j = 0; j < (*LinePointers[(*LinePointers[i]).ID]).Nodes.size(); j++) {
            int CurrentNode = (*LinePointers[(*LinePointers[i]).ID]).Nodes[j]; // Should be Nodes[j].ID
            std::string CurrentNodeString = std::to_string(CurrentNode);
            NodeList.append(CurrentNodeString);
            NodeList.append(", ");
        }
        NodeList.pop_back();
        NodeList.pop_back();
        std::cout << NodeList << std::endl;
    }
    std::cout << " " << std::endl;
}

void Geometry::PrintLineNodes(unsigned int id) {
    std::cout << "LINE NODES" << std::endl;
    std::cout << "Line " << id << ":" << std::endl;
    std::string NodeList;
    for (unsigned int j = 0; j < (*LinePointers[id]).Nodes.size(); j++) {
        int CurrentNode = (*LinePointers[id]).Nodes[j]; // Same as func above
        std::string CurrentNodeString = std::to_string(CurrentNode);
        NodeList.append(CurrentNodeString);
        NodeList.append(", ");
    }
    NodeList.pop_back();
    NodeList.pop_back();
    std::cout << NodeList << std::endl;
    std::cout << " " << std::endl;
}

void Geometry::Create4SidedArea(unsigned int Line0, unsigned int Line1, unsigned int Line2, unsigned int Line3) {
    std::vector<int> AreaLines = { (int)Line0, (int)Line1, (int)Line2, (int)Line3 };
    if (!AreaExistsInArray(&AreaLines) && AreaKeypointsUnique(&AreaLines)) {
        Eigen::RowVectorXi OrderedAreaLines = OrderAreaLines(&AreaLines);
        A.ID = Geometry::AssignAreaID();
        A.Lines << OrderedAreaLines;
        AreasArray.emplace_back(A);
        Geometry::AssignAreaPointer(A.ID);
        Geometry::AssignLineOwnersAndNeighbours(A.ID);
        Geometry::AssignAreaNeighbours(A.ID);
    }
}

bool Geometry::AreaExistsInArray(std::vector<int>* PotentialLines) {
    do {
        for (int i = 0; i < AreasArray.size(); i++) {
            std::vector<int> ExistingLines = PackageExistingLinesOfArea(i);
            if (*PotentialLines == ExistingLines) {
                return true;
            }
        }
    } while (std::next_permutation((*PotentialLines).begin(), (*PotentialLines).end()));
    return false;
}

std::vector<int> Geometry::PackageExistingLinesOfArea(int id) {
    std::vector<int> ExistingLines;
    for (int i = 0; i < A.Lines.size(); i++) {
        ExistingLines.push_back({ (*AreaPointers[id]).Lines(i) });
    }
    return ExistingLines;
}

// Checks if area of n lines only consist of n unique keypoints
bool Geometry::AreaKeypointsUnique(std::vector<int>* PotentialLines) {
    std::vector<int> Keypoints = Geometry::PackageKeypointsOfPotentialLines(PotentialLines);

    std::sort(Keypoints.begin(), Keypoints.end());
    Keypoints.erase(std::unique(Keypoints.begin(), Keypoints.end()), Keypoints.end());

    if (Keypoints.size() == (*PotentialLines).size()) {
        return true;
    }
    return false;
}

std::vector<int> Geometry::PackageKeypointsOfPotentialLines(std::vector<int>* PotentialLines) { // Could be used for any set of lines
    std::vector<int> Keypoints;
    for (int i = 0; i < (*PotentialLines).size(); i++) {
        Keypoints.push_back({ (int)(*LinePointers[(*PotentialLines)[i]]).Endpoints(0) });
        Keypoints.push_back({ (int)(*LinePointers[(*PotentialLines)[i]]).Endpoints(1) });
    }
    return Keypoints;
}

Eigen::RowVectorXi Geometry::OrderAreaLines(std::vector<int>* AreaLines) {
    std::vector<double> AverageXPosition;
    std::vector<double> AverageYPosition;
    std::tie(AverageXPosition, AverageYPosition) = CalculateAveragePositionsOfAreaLines(AreaLines);

    int IndexOfMaxX = std::max_element(AverageXPosition.begin(), AverageXPosition.end()) - AverageXPosition.begin(); // Could sort AvgXPos, Ypos and use index 0 and .back()
    int IndexOfMinX = std::min_element(AverageXPosition.begin(), AverageXPosition.end()) - AverageXPosition.begin();
    int IndexOfMaxY = std::max_element(AverageYPosition.begin(), AverageYPosition.end()) - AverageYPosition.begin();
    int IndexOfMinY = std::min_element(AverageYPosition.begin(), AverageYPosition.end()) - AverageYPosition.begin();

    Eigen::RowVectorXi OrderedAreaLines((*AreaLines).size());
    OrderedAreaLines << (*AreaLines)[IndexOfMinY], (*AreaLines)[IndexOfMaxY], (*AreaLines)[IndexOfMinX], (*AreaLines)[IndexOfMaxX];
    return OrderedAreaLines;
}

// EMPLACE_BACK, RESERVE ELEMENTTYPE
std::tuple<std::vector<double>, std::vector<double>> Geometry::CalculateAveragePositionsOfAreaLines(std::vector<int>* AreaLines) {
    std::vector<double> AverageXPosition;
    std::vector<double> AverageYPosition;
    for (int i = 0; i < (*AreaLines).size(); i++) {
        AverageXPosition.push_back({ 0.5 * ((*KeypointPointers[(*LinePointers[(*AreaLines)[i]]).Endpoints(0)]).Position(0) + (*KeypointPointers[(*LinePointers[(*AreaLines)[i]]).Endpoints(1)]).Position(0)) });
        AverageYPosition.push_back({ 0.5 * ((*KeypointPointers[(*LinePointers[(*AreaLines)[i]]).Endpoints(0)]).Position(1) + (*KeypointPointers[(*LinePointers[(*AreaLines)[i]]).Endpoints(1)]).Position(1)) });
    }
    return std::make_tuple(AverageXPosition, AverageYPosition);
}

unsigned int Geometry::AssignAreaID() {
    unsigned int AreaID;
    if (DeletedAreaIDs.size() > 0) {
        AreaID = DeletedAreaIDs[0];
        DeletedAreaIDs.erase(DeletedAreaIDs.begin());
        return AreaID;
    }
    else {
        AreaID = NextNewAreaID;
        Geometry::UpdateNextNewAreaID();
        return AreaID;
    }
}

void Geometry::UpdateNextNewAreaID() {
    if (AreasArray.size() != 0) {
        if (DeletedAreaIDs.size() == 0) {
            NextNewAreaID++;
        }
    }
    else {
        NextNewAreaID++;
    }
}

void Geometry::AssignAreaPointer(unsigned int id) {
    if (id == NextNewAreaID - 1) {
        AreaPointers.emplace_back(&(AreasArray.back()));
    }
    else {
        AreaPointers[id] = &(AreasArray.back());
    }
}

void Geometry::AssignLineOwnersAndNeighbours(unsigned int id) {
    for (int i = 0; i < (*AreaPointers[id]).Lines.size(); i++) {
        if ((*LinePointers[(*AreaPointers[id]).Lines(i)]).Owner == -1) {
            (*LinePointers[(*AreaPointers[id]).Lines(i)]).Owner = id;
        }
        else {
            (*LinePointers[(*AreaPointers[id]).Lines(i)]).Neighbours.emplace_back((*LinePointers[(*AreaPointers[id]).Lines(i)]).Owner);
        }
    }
}

// Overlapping lines need to be deleted if an area is deleted
void Geometry::AssignAreaNeighbours(unsigned int NewArea) {
    for (int i = 0; i < AreasArray.size() - 1; i++) {
        //Geometry::AreasShareLine(AreasArray[i].ID, NewArea);
        Geometry::AreaLinesOverlap(AreasArray[i].ID, NewArea);
    }
}

// Bool variant can be called with OverlappingLine Lines to assign nodes
void Geometry::AreasShareLine(unsigned int Area0, unsigned int Area1) {
    for (int i = 0; i < (*AreaPointers[Area0]).Lines.size(); i++) {
        for (int j = 0; j < (*AreaPointers[Area1]).Lines.size(); j++) {
            if ((*AreaPointers[Area0]).Lines(i) == (*AreaPointers[Area1]).Lines(j)) {
                // Void and write shared lines direct?
                (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).OverlappingLines.emplace_back((*AreaPointers[Area1]).Lines(j));
                (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).OverlappingLines.emplace_back((*AreaPointers[Area0]).Lines(i));
            }
        }
    }
}

// Don't need to do all permutations for 4 v 4 since North only ever borders South, East West (leaves 4 permutations NS, SN, EW, WE))
// [0, 2], [2, 0], [1, 3], [3, 1]

void Geometry::AreaLinesOverlap(unsigned int Area0, unsigned int Area1) {
    for (int i = 0; i < (*AreaPointers[Area0]).Lines.size(); i++) {
        for (int j = 0; j < (*AreaPointers[Area1]).Lines.size(); j++) {
            if (Geometry::LinesCollinear((*AreaPointers[Area0]).Lines(i), (*AreaPointers[Area1]).Lines(j))) {
                std::vector<double> KeypointXPositions;
                std::vector<double> KeypointYPositions;
                Eigen::VectorXi SharedKeypointIndices;
                Geometry::CheckIfLinesShareKeypoint(&SharedKeypointIndices, (*AreaPointers[Area0]).Lines(i), (*AreaPointers[Area1]).Lines(j));
                if (SharedKeypointIndices.size() == 2) { // Lines share a keypoint
                    int OtherIndex = Geometry::ReturnOtherEndpoint(SharedKeypointIndices(0));
                    unsigned int OtherKeypoint = (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).Endpoints(OtherIndex);
                    std::tie(KeypointXPositions, KeypointYPositions) = Geometry::PackageKeypointPositions((*AreaPointers[Area1]).Lines(j));
                    if (Geometry::SharedKeypointInRangeOfLine(KeypointXPositions, KeypointYPositions, OtherKeypoint)) {
                        (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).OverlappingLines.emplace_back((*AreaPointers[Area1]).Lines(j));
                        (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).OverlappingLines.emplace_back((*AreaPointers[Area0]).Lines(i));
                        //return true;
                    }
                    OtherIndex = Geometry::ReturnOtherEndpoint((*AreaPointers[Area1]).Lines(j)); // Switch to other point of Line1
                    OtherKeypoint = (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).Endpoints(OtherIndex);
                    std::tie(KeypointXPositions, KeypointYPositions) = Geometry::PackageKeypointPositions((*AreaPointers[Area0]).Lines(i));
                    if (Geometry::SharedKeypointInRangeOfLine(KeypointXPositions, KeypointYPositions, OtherKeypoint)) {
                        (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).OverlappingLines.emplace_back((*AreaPointers[Area1]).Lines(j));
                        (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).OverlappingLines.emplace_back((*AreaPointers[Area0]).Lines(i));
                        //return true;
                    }
                }
                else {
                    // Lines do not share a keypoint
                    std::tie(KeypointXPositions, KeypointYPositions) = Geometry::PackageKeypointPositions((*AreaPointers[Area1]).Lines(j));
                    if (Geometry::KeypointInRangeOfLine(KeypointXPositions, KeypointYPositions, (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).Endpoints(0))) {
                        (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).OverlappingLines.emplace_back((*AreaPointers[Area1]).Lines(j));
                        (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).OverlappingLines.emplace_back((*AreaPointers[Area0]).Lines(i));
                        //return true;
                    }
                    if (Geometry::KeypointInRangeOfLine(KeypointXPositions, KeypointYPositions, (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).Endpoints(1))) {
                        (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).OverlappingLines.emplace_back((*AreaPointers[Area1]).Lines(j));
                        (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).OverlappingLines.emplace_back((*AreaPointers[Area0]).Lines(i));
                        //return true;
                    }
                    // No common keypoint, first line not in range of other
                    std::tie(KeypointXPositions, KeypointYPositions) = Geometry::PackageKeypointPositions((*AreaPointers[Area0]).Lines(i));
                    if (Geometry::KeypointInRangeOfLine(KeypointXPositions, KeypointYPositions, (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).Endpoints(0))) {
                        (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).OverlappingLines.emplace_back((*AreaPointers[Area1]).Lines(j));
                        (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).OverlappingLines.emplace_back((*AreaPointers[Area0]).Lines(i));
                        //return true;
                    }
                    if (Geometry::KeypointInRangeOfLine(KeypointXPositions, KeypointYPositions, (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).Endpoints(1))) {
                        (*LinePointers[(*AreaPointers[Area0]).Lines(i)]).OverlappingLines.emplace_back((*AreaPointers[Area1]).Lines(j));
                        (*LinePointers[(*AreaPointers[Area1]).Lines(j)]).OverlappingLines.emplace_back((*AreaPointers[Area0]).Lines(i));
                        //return true;
                    }
                }
            }
        }
    }
    //return false;
}

bool Geometry::LinesCollinear(unsigned int Line0, unsigned int Line1) {
    std::vector<int> PotentialLines = { (int)Line0, (int)Line1 };
    std::vector<int> Keypoints = Geometry::PackageKeypointsOfPotentialLines(&PotentialLines);

    std::sort(Keypoints.begin(), Keypoints.end());
    Keypoints.erase(std::unique(Keypoints.begin(), Keypoints.end()), Keypoints.end());

    if (Keypoints.size() == 3) {
        Eigen::Matrix3d Collinear;
        Collinear << (*KeypointPointers[Keypoints[0]]).Position(0), (*KeypointPointers[Keypoints[0]]).Position(1), 1,
            (*KeypointPointers[Keypoints[1]]).Position(0), (*KeypointPointers[Keypoints[1]]).Position(1), 1,
            (*KeypointPointers[Keypoints[2]]).Position(0), (*KeypointPointers[Keypoints[2]]).Position(1), 1;
        if (abs(Collinear.determinant()) < 1e-17) {
            return true;
        }
        return false;
    }
}

void Geometry::CheckIfLinesShareKeypoint(Eigen::VectorXi* SharedKeypointIndices, unsigned int Line0, unsigned int Line1) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            if ((*LinePointers[Line0]).Endpoints(i) == (*LinePointers[Line1]).Endpoints(j)) {
                (*SharedKeypointIndices).resize(2);
                (*SharedKeypointIndices)(0) = i;
                (*SharedKeypointIndices)(1) = j;
                return;
            }
        }
    }
    (*SharedKeypointIndices).resize(3);
}

int Geometry::ReturnOtherEndpoint(unsigned int id) {
    return (id + 1) % 2;
}

std::tuple<std::vector<double>, std::vector<double>> Geometry::PackageKeypointPositions(unsigned int Line) {
    std::vector<double> KeypointXPositions;
    std::vector<double> KeypointYPositions;
    KeypointXPositions.reserve(2);
    KeypointYPositions.reserve(2);

    KeypointXPositions.emplace_back((*KeypointPointers[(*LinePointers[Line]).Endpoints(0)]).Position(0));
    KeypointXPositions.emplace_back((*KeypointPointers[(*LinePointers[Line]).Endpoints(1)]).Position(0));
    KeypointYPositions.emplace_back((*KeypointPointers[(*LinePointers[Line]).Endpoints(0)]).Position(1));
    KeypointYPositions.emplace_back((*KeypointPointers[(*LinePointers[Line]).Endpoints(1)]).Position(1));

    std::sort(KeypointXPositions.begin(), KeypointXPositions.end());
    std::sort(KeypointYPositions.begin(), KeypointYPositions.end());
    return std::make_tuple(KeypointXPositions, KeypointYPositions);
}

bool Geometry::SharedKeypointInRangeOfLine(std::vector<double> KeypointXPositions, std::vector<double> KeypointYPositions, unsigned int Keypoint) {
    if (Geometry::SharedXPosition(Keypoint, &KeypointXPositions)) {
        if ((*KeypointPointers[Keypoint]).Position(1) >= KeypointYPositions[0] && (*KeypointPointers[Keypoint]).Position(1) <= KeypointYPositions[1]) {
            return true;
        }
    }
    else if (Geometry::SharedYPosition(Keypoint, &KeypointYPositions)) {
        if ((*KeypointPointers[Keypoint]).Position(0) >= KeypointXPositions[0] && (*KeypointPointers[Keypoint]).Position(0) <= KeypointXPositions[1]) {
            return true;
        }
    }
    else {
        if (((*KeypointPointers[Keypoint]).Position(0) >= KeypointXPositions[0] && (*KeypointPointers[Keypoint]).Position(0) <= KeypointYPositions[0]) &&
            ((*KeypointPointers[Keypoint]).Position(1) >= KeypointYPositions[0] && (*KeypointPointers[Keypoint]).Position(1) <= KeypointYPositions[1])) {
            return true;
        }
    }
    return false;
}

bool Geometry::KeypointInRangeOfLine(std::vector<double> KeypointXPositions, std::vector<double> KeypointYPositions, unsigned int Keypoint) {
    if (Geometry::SharedXPosition(Keypoint, &KeypointXPositions)) {
        if ((*KeypointPointers[Keypoint]).Position(1) > KeypointYPositions[0] && (*KeypointPointers[Keypoint]).Position(1) < KeypointYPositions[1]) {
            return true; // Need different function for if sharedkeypoint since one value will always be equal, never passes &&
        }
    }
    else if (Geometry::SharedYPosition(Keypoint, &KeypointYPositions)) {
        if ((*KeypointPointers[Keypoint]).Position(0) > KeypointXPositions[0] && (*KeypointPointers[Keypoint]).Position(0) < KeypointXPositions[1]) {
            return true;
        }
    }
    else {
        if (((*KeypointPointers[Keypoint]).Position(0) > KeypointXPositions[0] && (*KeypointPointers[Keypoint]).Position(0) < KeypointYPositions[0]) &&
            ((*KeypointPointers[Keypoint]).Position(1) > KeypointYPositions[0] && (*KeypointPointers[Keypoint]).Position(1) < KeypointYPositions[1])) {
            return true;
        }
    }
    return false;
}

bool Geometry::SharedXPosition(unsigned int Keypoint, std::vector<double>* KeypointXPositions) {
    if ((*KeypointXPositions)[0] - (*KeypointXPositions)[1] == 0 && (*KeypointXPositions)[0] == (*KeypointPointers[Keypoint]).Position(0)) {
        return true;
    }
    return false;
}

bool Geometry::SharedYPosition(unsigned int Keypoint, std::vector<double>* KeypointYPositions) {
    if ((*KeypointYPositions)[0] - (*KeypointYPositions)[1] == 0 && (*KeypointYPositions)[0] == (*KeypointPointers[Keypoint]).Position(1)) {
        return true;
    }
    return false;
}

void Geometry::DeleteArea(unsigned int id) {
    Geometry::AddDeletedAreaID(id); // Could do with sorting when new added
    std::deque<Area>::iterator Location = std::find_if(AreasArray.begin(), AreasArray.end(), [id](Area& i) {return i.ID == id; });
    if (Location != AreasArray.end()) {
        int DistanceToStart = Location - AreasArray.begin();
        int DistanceToEnd = AreasArray.end() - Location - 1;
        AreasArray.erase(Location);  // Need to delete elements and nodes if present, nodes also need erasing from lines
        Geometry::AdjustAreaPointers(id, DistanceToStart, DistanceToEnd);
    }
}

void Geometry::AddDeletedAreaID(unsigned int id) {
    DeletedAreaIDs.emplace_back(id);
}

void Geometry::AdjustAreaPointers(unsigned int id, int DistanceToStart, int DistanceToEnd) {
    if (DistanceToStart != 0 && DistanceToEnd != 0) {
        int Size = AreaPointers.size() - 1;
        if (DistanceToStart < DistanceToEnd) {
            delete(AreaPointers[0]);
            for (int i = 0; i < DistanceToStart; i++) {
                AreaPointers[i] = AreaPointers[i + 1];
            }
            AreaPointers[id] = nullptr;
        }
        else {
            delete(AreaPointers.back());
            for (int i = 0; i < DistanceToEnd; i++) {
                AreaPointers[Size - i] = AreaPointers[Size - (i + 1)];
            }
            AreaPointers[id] = nullptr;
        }
    }
    else {
        AreaPointers[id] = nullptr;
    }
}

void Geometry::PrintAreaLines() {
    std::cout << "AREA LINES" << std::endl;
    for (unsigned int i = 0; i < AreasArray.size(); i++) {
        std::cout << "Area " << AreasArray[i].ID << ":" << std::endl;
        std::cout << AreasArray[i].Lines(0) << ", " << AreasArray[i].Lines(1) <<
            ", " << AreasArray[i].Lines(2) << ", " << AreasArray[i].Lines(3) << std::endl;
    }
    std::cout << " " << std::endl;
}

void Geometry::PrintAreaLines(unsigned int id) {
    std::cout << "AREA LINES" << std::endl;
    std::cout << "Area " << id << ":" << std::endl;
    std::cout << (*AreaPointers[id]).Lines(0) << ", " << (*AreaPointers[id]).Lines(1) <<
        ", " << (*AreaPointers[id]).Lines(2) << ", " << (*AreaPointers[id]).Lines(3) << std::endl;
    std::cout << "" << std::endl;
}

void Geometry::PrintAreaNodes() {
    std::cout << "AREA NODES" << std::endl;
    for (unsigned int i = 0; i < AreasArray.size(); i++) {
        std::cout << "Area " << AreasArray[i].ID << ":" << std::endl;
        std::string NodeList;
        for (unsigned int j = 0; j < AreasArray[i].Nodes.size(); j++) {
            int CurrentNode = AreasArray[i].Nodes[j];
            std::string CurrentNodeString = std::to_string(CurrentNode);
            NodeList.append(CurrentNodeString);
            NodeList.append(", ");
        }
        NodeList.pop_back();
        NodeList.pop_back();
        std::cout << NodeList << std::endl;
    }
    std::cout << " " << std::endl;
}

void Geometry::PrintAreaNodes(unsigned int id) {
    std::cout << "AREA NODES" << std::endl;
    std::cout << "Area " << id << ":" << std::endl;
    std::string NodeList;
    for (unsigned int j = 0; j < (*AreaPointers[id]).Nodes.size(); j++) {
        int CurrentNode = (*AreaPointers[id]).Nodes[j]; // Same as above
        std::string CurrentNodeString = std::to_string(CurrentNode);
        NodeList.append(CurrentNodeString);
        NodeList.append(", ");
    }
    NodeList.pop_back();
    NodeList.pop_back();
    std::cout << NodeList << std::endl;
    std::cout << " " << std::endl;
}

void Geometry::GenerateAreaMesh(unsigned int id) {
    if (Geometry::AreaLineDivisionsValid(id)) {
        Eigen::MatrixXi AreaLines = Geometry::PackageAreaLineEndpoints(id);
        std::vector<Eigen::MatrixXd> AreaLineNodes = Geometry::CreateAreaLineNodes(id, &AreaLines);
        std::vector<Eigen::MatrixXd> NodePositions = Geometry::InitialiseAreaNodes(id);
        Geometry::InputAreaBoundaryNodes(&AreaLineNodes, &NodePositions);
        Geometry::InitialiseAreaGeometry(&NodePositions);
        Geometry::CreateAreaNodes(id, &NodePositions); // Boundary nodes may not need doing if line has neighbour
        Geometry::CreateAreaElements(id);
        Geometry::CheckMeshQuality();
    }
}

bool Geometry::AreaLineDivisionsValid(unsigned int id) {
    if ((*LinePointers[(*AreaPointers[id]).Lines(0)]).Divisions == (*LinePointers[(*AreaPointers[id]).Lines(1)]).Divisions
        && (*LinePointers[(*AreaPointers[id]).Lines(2)]).Divisions == (*LinePointers[(*AreaPointers[id]).Lines(3)]).Divisions) {
        return true;
    }
    return false;
}

Eigen::MatrixXi Geometry::PackageAreaLineEndpoints(unsigned int id) {
    Eigen::MatrixXi Lines((*AreaPointers[id]).Lines.size(), 2); // Might throw up error for unsigned int "size"
    Lines << (*LinePointers[(*AreaPointers[id]).Lines(0)]).Endpoints(0), (*LinePointers[(*AreaPointers[id]).Lines(0)]).Endpoints(1),
        (*LinePointers[(*AreaPointers[id]).Lines(1)]).Endpoints(0), (*LinePointers[(*AreaPointers[id]).Lines(1)]).Endpoints(1),
        (*LinePointers[(*AreaPointers[id]).Lines(2)]).Endpoints(0), (*LinePointers[(*AreaPointers[id]).Lines(2)]).Endpoints(1),
        (*LinePointers[(*AreaPointers[id]).Lines(3)]).Endpoints(0), (*LinePointers[(*AreaPointers[id]).Lines(3)]).Endpoints(1);
    return Lines;
}

std::vector<Eigen::MatrixXd> Geometry::CreateAreaLineNodes(unsigned int id, Eigen::MatrixXi* AreaLines) {
    unsigned int divX = (*LinePointers[(*AreaPointers[id]).Lines(0)]).Divisions;
    unsigned int divY = (*LinePointers[(*AreaPointers[id]).Lines(2)]).Divisions;

    Eigen::MatrixXd BoundaryNS(divX + 1, 4);
    Eigen::MatrixXd BoundaryEW(divY + 1, 4);
    Eigen::MatrixXd AreaLineS(divX + 1, 2), AreaLineW(divY + 1, 2),
        AreaLineN(divX + 1, 2), AreaLineE(divY + 1, 2);
    std::vector<Eigen::MatrixXd> AreaLineNodes;

    // Constructs line for each pair of points, PUSH_BACK!
    for (int i = 0; i < (*AreaPointers[id]).Lines.size(); i++) {
        Eigen::RowVectorXd EndpointPositions = Geometry::PackageAreaLineEndpointPositions((*AreaLines).row(i));
        Eigen::RowVectorXd ConstantLineParameters = Geometry::CalculateConstantLineParameters(&EndpointPositions);
        Eigen::VectorXd SpacingRatioIntervals = Geometry::GenerateSpacingRatioIntervals((*LinePointers[(*AreaPointers[id]).Lines(i)]).Divisions, (*LinePointers[(*AreaPointers[id]).Lines(i)]).SpacingRatio);

        if (Geometry::AreaLineGradientFinite(&EndpointPositions)) {
            if (Geometry::ComputingNSLines(i)) {

                Geometry::CreateNSLineNodesFinite(&EndpointPositions, &ConstantLineParameters, &SpacingRatioIntervals, &BoundaryNS, i);
            }
            else {

                Geometry::CreateEWLineNodesFinite(&EndpointPositions, &ConstantLineParameters, &SpacingRatioIntervals, &BoundaryEW, i);
            }
        }
        else {
            if (Geometry::ComputingNSLines(i)) {

                Geometry::CreateNSLineNodesInfinite(&EndpointPositions, &SpacingRatioIntervals, &BoundaryNS, i);
            }
            else {

                Geometry::CreateEWLineNodesInfinite(&EndpointPositions, &SpacingRatioIntervals, &BoundaryEW, i);
            }
        }

        if (Geometry::ComputingNSLines(i)) {
            BoundaryNS(divX, 2 * i) = EndpointPositions(2);
            BoundaryNS(divX, 2 * i + 1) = EndpointPositions(3);
        }
        else {
            BoundaryEW(divY, 2 * (i - 2)) = EndpointPositions(2);
            BoundaryEW(divY, 2 * (i - 2) + 1) = EndpointPositions(3);
        }

    }

    AreaLineS.col(0) = BoundaryNS.col(0);
    AreaLineS.col(1) = BoundaryNS.col(1);
    AreaLineW.col(0) = BoundaryEW.col(0);
    AreaLineW.col(1) = BoundaryEW.col(1);
    AreaLineN.col(0) = BoundaryNS.col(2);
    AreaLineN.col(1) = BoundaryNS.col(3);
    AreaLineE.col(0) = BoundaryEW.col(2);
    AreaLineE.col(1) = BoundaryEW.col(3);

    AreaLineNodes.push_back({ AreaLineS });
    AreaLineNodes.push_back({ AreaLineN });
    AreaLineNodes.push_back({ AreaLineW });
    AreaLineNodes.push_back({ AreaLineE });
    return AreaLineNodes;
}

Eigen::RowVectorXd Geometry::PackageAreaLineEndpointPositions(Eigen::RowVectorXi AreaLine) {
    Eigen::RowVectorXd EndpointPositions(4);
    double x0 = (*KeypointPointers[(AreaLine)(0)]).Position(0);
    double y0 = (*KeypointPointers[(AreaLine)(0)]).Position(1);
    double x1 = (*KeypointPointers[(AreaLine)(1)]).Position(0);
    double y1 = (*KeypointPointers[(AreaLine)(1)]).Position(1);

    EndpointPositions << x0, y0, x1, y1;
    return EndpointPositions;
}

bool Geometry::AreaLineGradientFinite(Eigen::RowVectorXd* EndpointPositions) {
    if ((*EndpointPositions)(0) != (*EndpointPositions)(2)) {
        return true;
    }
    return false;
}

Eigen::RowVectorXd Geometry::CalculateConstantLineParameters(Eigen::RowVectorXd* EndpointPositions) {
    Eigen::RowVectorXd ConstantLineParameters(4);
    double m = ((*EndpointPositions)(3) - (*EndpointPositions)(1)) / ((*EndpointPositions)(2) - (*EndpointPositions)(0));
    double R = pow(pow((*EndpointPositions)(2) - (*EndpointPositions)(0), 2) + pow((*EndpointPositions)(3) - (*EndpointPositions)(1), 2), 0.5);
    double a = 1.0 + pow(m, 2);
    double b = 2.0 * (*EndpointPositions)(0) * (1.0 - m);
    ConstantLineParameters << m, R, a, b;
    return ConstantLineParameters;
}

Eigen::VectorXd Geometry::GenerateSpacingRatioIntervals(unsigned int div, float SpacingRatio) { // Could take out first term
    Eigen::VectorXd SpacingRatioIntervals = Eigen::VectorXd::Zero(div);
    if (SpacingRatio != 1.0) {
        Eigen::VectorXd SpacingRatioDistances = Eigen::VectorXd::Zero(div + 1);
        for (int i = 1; i < div + 1; i++) { // Works for SR < 1 also
            SpacingRatioDistances(i) = div - 1 + (SpacingRatio - 1) * (i - 1);
        }
        SpacingRatioDistances /= (double(div) - 1);

        double TotalFraction = SpacingRatioDistances.sum();
        SpacingRatioIntervals = SpacingRatioDistances.block(0, 0, div, 1);

        int i = 2;
        while (i < div) {
            SpacingRatioIntervals(i) += SpacingRatioDistances.block(1, 0, i - 1, 1).sum();
            i++;
        }
        SpacingRatioIntervals /= TotalFraction;
    }
    else {
        for (int i = 0; i < div; i++) {
            SpacingRatioIntervals(i) = i / div;
        }
    }
    return SpacingRatioIntervals;
}

bool Geometry::ComputingNSLines(int i) {
    if (i < 2) {
        return true;
    }
    return false;
}

void Geometry::CreateNSLineNodesFinite(Eigen::RowVectorXd* EndpointPositions, Eigen::RowVectorXd* ConstantLineParameters, Eigen::VectorXd* SpacingRatioIntervals, Eigen::MatrixXd* BoundaryNS, int i) {
    unsigned int divX = (*BoundaryNS).rows() - 1;
    double c;
    for (unsigned int j = 0; j < divX; j++) {
        c = -pow((*ConstantLineParameters)(1) * (*SpacingRatioIntervals)(j), 2);
        (*BoundaryNS)(j, 2 * i) = (-(*ConstantLineParameters)(3) + pow(pow((*ConstantLineParameters)(3), 2) - 4.0 * (*ConstantLineParameters)(2) * c, 0.5)) / (2.0 * (*ConstantLineParameters)(2));
        (*BoundaryNS)(j, 2 * i + 1) = (*ConstantLineParameters)(0) * ((*BoundaryNS)(j, 2 * i) - (*EndpointPositions)(0)) + (*EndpointPositions)(1);
    }
}

void Geometry::CreateEWLineNodesFinite(Eigen::RowVectorXd* EndpointPositions, Eigen::RowVectorXd* ConstantLineParameters, Eigen::VectorXd* SpacingRatioIntervals, Eigen::MatrixXd* BoundaryEW, int i) {
    unsigned int divY = (*BoundaryEW).rows() - 1;
    double c;
    for (unsigned int j = 0; j < divY; j++) {
        c = -pow((*ConstantLineParameters)(1) * (*SpacingRatioIntervals)(j), 2);
        (*BoundaryEW)(j, 2 * (i - 2)) = (-(*ConstantLineParameters)(3) + pow(pow((*ConstantLineParameters)(3), 2) - 4.0 * (*ConstantLineParameters)(2) * c, 0.5)) / (2.0 * (*ConstantLineParameters)(2));
        (*BoundaryEW)(j, 2 * (i - 2) + 1) = (*ConstantLineParameters)(0) * ((*BoundaryEW)(j, 2 * (i - 2)) - (*EndpointPositions)(0)) + (*EndpointPositions)(1);
    }
}

void Geometry::CreateNSLineNodesInfinite(Eigen::RowVectorXd* EndpointPositions, Eigen::VectorXd* SpacingRatioIntervals, Eigen::MatrixXd* BoundaryNS, int i) {
    unsigned int divX = (*BoundaryNS).rows() - 1;
    for (unsigned int j = 0; j < divX; j++) {
        (*BoundaryNS)(j, 2 * i) = (*EndpointPositions)(0);
        (*BoundaryNS)(j, 2 * i + 1) = (*EndpointPositions)(1) + ((*EndpointPositions)(3) - (*EndpointPositions)(1)) * (*SpacingRatioIntervals)(j);
    }
}

void Geometry::CreateEWLineNodesInfinite(Eigen::RowVectorXd* EndpointPositions, Eigen::VectorXd* SpacingRatioIntervals, Eigen::MatrixXd* BoundaryEW, int i) {
    unsigned int divY = (*BoundaryEW).rows() - 1;
    for (unsigned int j = 0; j < divY; j++) {
        (*BoundaryEW)(j, 2 * (i - 2)) = (*EndpointPositions)(0);
        (*BoundaryEW)(j, 2 * (i - 2) + 1) = (*EndpointPositions)(1) + ((*EndpointPositions)(3) - (*EndpointPositions)(1)) * (*SpacingRatioIntervals)(j);
    }
}

// Resize the nodal matrices of an area
std::vector<Eigen::MatrixXd> Geometry::InitialiseAreaNodes(unsigned int id) {
    std::vector<Eigen::MatrixXd> NodePositions;
    NodePositions.reserve(2);
    Eigen::MatrixXd NodesX((*LinePointers[(*AreaPointers[id]).Lines(2)]).Divisions + 1, (*LinePointers[(*AreaPointers[id]).Lines(0)]).Divisions + 1);
    Eigen::MatrixXd NodesY((*LinePointers[(*AreaPointers[id]).Lines(2)]).Divisions + 1, (*LinePointers[(*AreaPointers[id]).Lines(0)]).Divisions + 1);
    NodePositions.emplace_back(NodesX);
    NodePositions.emplace_back(NodesY);
    /*
    for (int i = 0; i < NodesX.size(); i++) {
        AreasArray[id].Nodes.push_back({ NodesArray.size() + i }); // NOT TO BE USED AFTER DELETING AREAS/NODES
    }

    int NumberOfNodesToLastLine = NodesX.cols() * (NodesX.rows() - 1);
    for (int i = 0; i < NodesX.cols(); i++) {
        LinesArray[AreasArray[id].Lines(0)].Nodes.emplace_back(NodesArray.size() + i );
        LinesArray[AreasArray[id].Lines(1)].Nodes.emplace_back(NodesArray.size() + NumberOfNodesToLastLine + i);
    }

    int NumberOfNodesAlongX = NodesX.cols();
    for (int i = 0; i < NodesX.rows(); i++) {
        LinesArray[AreasArray[id].Lines(2)].Nodes.emplace_back(NodesArray.size() + (NumberOfNodesAlongX) * i);
        LinesArray[AreasArray[id].Lines(3)].Nodes.emplace_back(NodesArray.size() + (NumberOfNodesAlongX) *  i + NumberOfNodesAlongX - 1);
    }
    */
    return NodePositions;
}

void Geometry::InputAreaBoundaryNodes(std::vector<Eigen::MatrixXd>* AreaLineNodes, std::vector<Eigen::MatrixXd>* NodePositions) {
    unsigned int divX = (*NodePositions)[0].cols() - 1;
    unsigned int divY = (*NodePositions)[0].rows() - 1;

    for (unsigned int i = 0; i < divX + 1; i++) {
        (*NodePositions)[0](0, i) = (*AreaLineNodes)[0](i, 0);
        (*NodePositions)[1](0, i) = (*AreaLineNodes)[0](i, 1);
        (*NodePositions)[0](divY, i) = (*AreaLineNodes)[1](i, 0);
        (*NodePositions)[1](divY, i) = (*AreaLineNodes)[1](i, 1);
    }
    for (unsigned int i = 0; i < divY + 1; i++) {
        (*NodePositions)[0](i, 0) = (*AreaLineNodes)[2](i, 0);
        (*NodePositions)[1](i, 0) = (*AreaLineNodes)[2](i, 1);
        (*NodePositions)[0](i, divX) = (*AreaLineNodes)[3](i, 0);
        (*NodePositions)[1](i, divX) = (*AreaLineNodes)[3](i, 1);
    }
}

// Constructs initial geometry matrix using Transfinite Interpolation, next would be smoothen or poisson/elliptical smoothen
void Geometry::InitialiseAreaGeometry(std::vector<Eigen::MatrixXd>* NodePositions) {
    unsigned int divX = (*NodePositions)[0].cols() - 1;
    unsigned int divY = (*NodePositions)[0].rows() - 1;
    double x, y;
    for (unsigned int i = 1; i < divX; i++) {
        x = i / (divX - 1.0);
        for (unsigned int j = 1; j < divY; j++) {
            y = j / (divY - 1.0);
            (*NodePositions)[0](j, i) = (1.0 - y) * (*NodePositions)[0](j, 0) + y * (*NodePositions)[0](j, divX) + (1.0 - x) * (*NodePositions)[0](0, i) + x * (*NodePositions)[0](divY, i) -
                ((1.0 - x) * (1.0 - y) * (*NodePositions)[0](0, 0) + x * y * (*NodePositions)[0](divY, divX) + x * (1.0 - y) * (*NodePositions)[0](divY, 0) + (1.0 - x) * y * (*NodePositions)[0](0, divX));
            (*NodePositions)[1](j, i) = (1.0 - y) * (*NodePositions)[1](j, 0) + y * (*NodePositions)[1](j, divX) + (1.0 - x) * (*NodePositions)[1](0, i) + x * (*NodePositions)[1](divY, i) -
                ((1.0 - x) * (1.0 - y) * (*NodePositions)[1](0, 0) + x * y * (*NodePositions)[1](divY, divX) + x * (1.0 - y) * (*NodePositions)[1](divY, 0) + (1.0 - x) * y * (*NodePositions)[1](0, divX));
        }
    }
}

void Geometry::CreateAreaNodes(unsigned int id, std::vector<Eigen::MatrixXd>* NodePositions) {
    int Xsize = (*NodePositions)[0].cols();
    int Ysize = (*NodePositions)[0].rows();
    for (int i = 0; i < Ysize; i++) {
        for (int j = 0; j < Xsize; j++) {
            Geometry::CreateNode((*NodePositions)[0](i, j), (*NodePositions)[1](i, j));
            (*AreaPointers[id]).Nodes.emplace_back(NodesArray.back().ID);
            if (i == 0) {
                (*LinePointers[(*AreaPointers[id]).Lines(0)]).Nodes.emplace_back(NodesArray.back().ID);
            }
            if (i == Ysize - 1) {
                (*LinePointers[(*AreaPointers[id]).Lines(1)]).Nodes.emplace_back(NodesArray.back().ID);
            }
            if (j == 0) {
                (*LinePointers[(*AreaPointers[id]).Lines(2)]).Nodes.emplace_back(NodesArray.back().ID);
            }
            if (j == Xsize - 1) {
                (*LinePointers[(*AreaPointers[id]).Lines(3)]).Nodes.emplace_back(NodesArray.back().ID);
            }
        }
    }
    Geometry::UpdateNodeCount();
}

void Geometry::CreateNode(double x, double y) {
    N.ID = Geometry::AssignNodeID();
    N.Position << x, y;
    NodesArray.emplace_back(N);
    Geometry::AssignNodePointer(N.ID);
}

unsigned int Geometry::AssignNodeID() {
    unsigned int NodeID;
    if (DeletedNodeIDs.size() > 0) {
        NodeID = DeletedNodeIDs[0];
        DeletedNodeIDs.erase(DeletedNodeIDs.begin());
        return NodeID;
    }
    else {
        NodeID = NextNewNodeID;
        Geometry::UpdateNextNewNodeID();
        return NodeID;
    }
}

void Geometry::UpdateNextNewNodeID() {
    if (NodesArray.size() != 0) {
        if (DeletedNodeIDs.size() == 0) {
            NextNewNodeID++;
        }
    }
    else {
        NextNewNodeID++;
    }
}

void Geometry::AssignNodePointer(unsigned int id) {
    if (id == NextNewNodeID - 1) {
        NodePointers.emplace_back(&(NodesArray.back()));
    }
    else {
        NodePointers[id] = &(NodesArray.back());
    }
}

void Geometry::UpdateNodeCount() {
    NodeCount = NodesArray.size();
}

void Geometry::PrintNodes() {
    std::cout << "NODE POSITIONS" << std::endl;
    for (unsigned int i = 0; i < NodesArray.size(); i++) {
        std::cout << "Node " << NodesArray[i].ID << ":" << std::endl;
        std::cout << NodesArray[i].Position(0) << ", " << NodesArray[i].Position(1) << std::endl;
    }
    std::cout << "" << std::endl;
}

void Geometry::PrintNodes(unsigned int id) {
    std::cout << "NODE POSITION" << std::endl;
    std::cout << "Node " << id << ":" << std::endl;
    std::cout << (*NodePointers[id]).Position(0) << ", " << (*NodePointers[id]).Position(1) << std::endl;
    std::cout << "" << std::endl;
}

void Geometry::InitialiseElements() {
    unsigned int StrainMatrixSize = ReturnStrainMatrixSize();
    E.Nodes.resize(ElementType);
    E.Strains.resize(ElementType, StrainMatrixSize);
    E.Stresses.resize(ElementType, StrainMatrixSize + 1);
}

unsigned int Geometry::ReturnStrainMatrixSize() {
    unsigned int result;
    if (Case == 0 || Case == 1) {
        result = 3;
        return result;
    }
    else if (Case == 2) {
        result = 6;
        return result;
    }
}

void Geometry::CreateAreaElements(unsigned int id) {
    if (Geometry::ElementsInitialised()) {
        int NodesX = (*LinePointers[(*AreaPointers[id]).Lines(0)]).Divisions + 1;
        int NodesY = (*LinePointers[(*AreaPointers[id]).Lines(2)]).Divisions + 1;
        Eigen::RowVectorXi ElementNodes(E.Nodes.size());
        for (int i = 0; i < NodesY - 1; i++) {
            for (int j = 0; j < NodesX - 1; j++) {
                ElementNodes(0) = (*AreaPointers[id]).Nodes[i * NodesX + j];
                ElementNodes(1) = ElementNodes(0) + 1;
                ElementNodes(2) = (*AreaPointers[id]).Nodes[(i + 1) * NodesX + j];
                ElementNodes(3) = ElementNodes(2) + 1;
                CreateElement(&ElementNodes);
            }
        }
    }
    Geometry::UpdateElementCount();
}

bool Geometry::ElementsInitialised() {
    if (E.Nodes.size() == ElementType) {
        return true;
    }
    return false;
}

void Geometry::CreateElement(Eigen::RowVectorXi* ElementNodes) {
    E.ID = Geometry::AssignElementID();
    E.Nodes = *ElementNodes;
    ElementsArray.emplace_back(E);
    Geometry::AssignElementPointer(E.ID);
}

unsigned int Geometry::AssignElementID() {
    unsigned int ElementID;
    if (DeletedElementIDs.size() > 0) {
        ElementID = DeletedElementIDs[0];
        DeletedElementIDs.erase(DeletedElementIDs.begin());
        return ElementID;
    }
    else {
        ElementID = NextNewElementID;
        Geometry::UpdateNextNewElementID();
        return ElementID;
    }
}

void Geometry::UpdateNextNewElementID() {
    if (ElementsArray.size() != 0) {
        if (DeletedElementIDs.size() == 0) {
            NextNewElementID++;
        }
    }
    else {
        NextNewElementID++;
    }
}

void Geometry::AssignElementPointer(unsigned int id) {
    if (id == NextNewElementID - 1) {
        ElementPointers.emplace_back(&(ElementsArray.back()));
    }
    else {
        ElementPointers[id] = &(ElementsArray.back());
    }
}

void Geometry::UpdateElementCount() {
    ElementCount = ElementsArray.size();
}

void Geometry::CheckMeshQuality() {
    Eigen::MatrixXd ElementLineVectors(ElementType, 2);
    Eigen::VectorXd ElementLineLengths(ElementType);
    Eigen::VectorXd LineAngles(2);
    Eigen::VectorXd ElementAspectRatios(ElementCount);
    Eigen::VectorXd ElementAreas(ElementCount);
    Eigen::VectorXd ElementSkewnesses(ElementCount);
    double AverageAspectRatio, MaxAspectRatio, AverageArea, MaxArea, AverageSkewness, MaxSkewness;

    for (int i = 0; i < ElementCount; i++) {
        Geometry::ReturnElementLineVectors(&ElementLineVectors, i);
        Geometry::ReturnElementLineLengths(&ElementLineVectors, &ElementLineLengths);
        Geometry::ReturnLineAngles(&ElementLineVectors, &ElementLineLengths, &LineAngles);
        Geometry::CheckAspectRatio(&ElementAspectRatios, &ElementLineLengths, i);
        Geometry::CheckSkewness(&ElementSkewnesses, &ElementLineLengths, &ElementAreas, i);
    }
    AverageAspectRatio = ElementAspectRatios.mean();
    MaxAspectRatio = ElementAspectRatios.maxCoeff();
    AverageArea = ElementAreas.mean();
    MaxArea = ElementAreas.maxCoeff();
    AverageSkewness = ElementSkewnesses.mean();
    MaxSkewness = ElementSkewnesses.maxCoeff();

    std::cout << "MESH QUALITY CHECK" << std::endl;
    std::cout << "Average element aspect ratio: " << AverageAspectRatio << std::endl;
    std::cout << "Maximum element aspect ratio: " << MaxAspectRatio << std::endl;
    std::cout << "Average element area: " << AverageArea << std::endl;
    std::cout << "Maximum element area: " << MaxArea << std::endl;
    std::cout << "Average element skewness: " << AverageSkewness << std::endl;
    std::cout << "Maximum element skewness: " << MaxSkewness << std::endl;
    std::cout << " " << std::endl;
}

void Geometry::ReturnElementLineVectors(Eigen::MatrixXd* LineVectors, unsigned int id) {
    (*LineVectors).row(0) = NodesArray[ElementsArray[id].Nodes(1)].Position - NodesArray[ElementsArray[id].Nodes(0)].Position;
    (*LineVectors).row(1) = NodesArray[ElementsArray[id].Nodes(2)].Position - NodesArray[ElementsArray[id].Nodes(3)].Position;
    (*LineVectors).row(2) = NodesArray[ElementsArray[id].Nodes(2)].Position - NodesArray[ElementsArray[id].Nodes(0)].Position;
    (*LineVectors).row(3) = NodesArray[ElementsArray[id].Nodes(1)].Position - NodesArray[ElementsArray[id].Nodes(3)].Position;
}

void Geometry::ReturnElementLineLengths(Eigen::MatrixXd* LineVectors, Eigen::VectorXd* LineLengths) { // Only good for 4 noded element
    (*LineLengths)(0) = (*LineVectors).row(0).norm();
    (*LineLengths)(1) = (*LineVectors).row(1).norm();
    (*LineLengths)(2) = (*LineVectors).row(2).norm();
    (*LineLengths)(3) = (*LineVectors).row(3).norm();
}

void Geometry::ReturnElementArea(Eigen::VectorXd* LineLengths, Eigen::VectorXd* LineAngles, Eigen::VectorXd* Areas, unsigned int id) {
    (*Areas)(id) = 0.5 * ((*LineLengths)(0) * (*LineLengths)(2) * sin((*LineAngles)(0)) + (*LineLengths)(1) * (*LineLengths)(3) * sin((*LineAngles)(1)));
}

void Geometry::CheckAspectRatio(Eigen::VectorXd* AspectRatios, Eigen::VectorXd* LineLengths, unsigned int id) {
    double AspectRatio = ((*LineLengths)(0) + (*LineLengths)(1)) / ((*LineLengths)(2) + (*LineLengths)(3));
    if (AspectRatio < 1) {
        AspectRatio = pow(AspectRatio, -1);
    }
    (*AspectRatios)(id) = AspectRatio;

} // Would be nice to know percentage of elements with high ARs

void Geometry::ReturnLineAngles(Eigen::MatrixXd* LineVectors, Eigen::VectorXd* LineLengths, Eigen::VectorXd* LineAngles) {
    (*LineAngles)(0) = acos(((*LineVectors)(0) * (*LineVectors)(2)) / ((*LineLengths)(0) * (*LineLengths)(2)));
    (*LineAngles)(1) = acos(((*LineVectors)(3) * (*LineVectors)(1)) / ((*LineLengths)(3) * (*LineLengths)(1)));
}

void Geometry::CheckSkewness(Eigen::VectorXd* Skewnesses, Eigen::VectorXd* LineLengths, Eigen::VectorXd* Areas, unsigned int id) {
    double SemiPerimeter = 0.5 * (*LineLengths).sum();
    (*Areas)(id) = pow(((SemiPerimeter - (*LineLengths)(0)) * (SemiPerimeter - (*LineLengths)(1)) * (SemiPerimeter - (*LineLengths)(2)) * (SemiPerimeter - (*LineLengths)(3))), 0.5);
    double R = 0.25 * pow(((*LineLengths)(0) * (*LineLengths)(3) + (*LineLengths)(1) * (*LineLengths)(2))
        * ((*LineLengths)(0) * (*LineLengths)(1) + (*LineLengths)(2) * (*LineLengths)(3)) * ((*LineLengths)(0)
            * (*LineLengths)(2) + (*LineLengths)(1) * (*LineLengths)(3)), 0.5) / (*Areas)(id);
    double PerfectArea = 2.0 * pow(R, 2);
    (*Skewnesses)(id) = (PerfectArea - (*Areas)(id)) / PerfectArea;
} // Acceptable under ~0.8

void Geometry::PrintElements() {
    std::cout << "ELEMENT NODES" << std::endl;
    for (unsigned int i = 0; i < ElementsArray.size(); i++) {
        std::cout << "Element " << ElementsArray[i].ID << ":" << std::endl;
        std::cout << ElementsArray[i].Nodes(0) << ", " << ElementsArray[i].Nodes(1) <<
            ", " << ElementsArray[i].Nodes(2) << ", " << ElementsArray[i].Nodes(3) << std::endl;
    }
    std::cout << "" << std::endl;
}

void Geometry::PrintElements(unsigned int id) {
    std::cout << "ELEMENT NODES" << std::endl;
    std::cout << "Element " << id << ":" << std::endl;
    std::cout << (*ElementPointers[id]).Nodes(0) << ", " << (*ElementPointers[id]).Nodes(1) << ", " << (*ElementPointers[id]).Nodes(2) << ", " << (*ElementPointers[id]).Nodes(3) << std::endl;
    std::cout << "" << std::endl;
}
