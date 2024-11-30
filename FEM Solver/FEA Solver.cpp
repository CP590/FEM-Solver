#include <iostream>
#include <iterator>
#include <chrono>
#include "Geometry.h"
#include "Analysis.h"


double E = 200e9;
double nu = 0.3;
double thickness = 36e-3;

unsigned int ElementType = 4;
unsigned int Case = 0;


int main() {
    auto start = std::chrono::steady_clock::now();

    Geometry::CreateKeypoint(0.0, 0.0);
    Geometry::CreateKeypoint(800e-3, 0.0);
    Geometry::CreateKeypoint(0.0, 42e-3);
    Geometry::CreateKeypoint(800e-3, 42e-3);
    //Geometry::CreateKeypoint(1.0, 0.0);
    //Geometry::CreateKeypoint(800e-3, 100e-3);
    //Geometry::CreateKeypoint(1.0, 100e-3);
    Geometry::CreateLine(0, 1);
    Geometry::CreateLine(2, 3);
    Geometry::CreateLine(0, 2);
    Geometry::CreateLine(1, 3);

    //Geometry::CreateLine(1, 4);
    //Geometry::CreateLine(5, 6);
    //Geometry::CreateLine(1, 5);
    //Geometry::CreateLine(4, 6);
    Geometry::AssignLineDivision(0, 5);
    Geometry::AssignLineDivision(1, 5);
    Geometry::AssignLineDivision(2, 1);
    Geometry::AssignLineDivision(3, 1);
    //Geometry::AssignLineSpacingRatio(0, 0.5);
    //Geometry::AssignLineSpacingRatio(1, 0.5);
    Geometry::Create4SidedArea(3, 1, 2, 0);
    //Geometry::Create4SidedArea(4, 5, 6, 7);
    //Geometry::AssignAreaNeighbours(1);
    Geometry::InitialiseElements();
    Geometry::GenerateAreaMesh(0);
    Geometry::PrintNodes();
    ApplyConstraintsToLine(2, true, true);
    ApplyYForcesToLine(3, -1050.0);
    SetMaterial(E, nu, thickness);
    SolveSystem();
    GenerateNodeStrainsAndStresses();
    PrintNodeDisplacements();
    //PrintElementStrains(0);
    PrintElementStresses();

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_time = end - start;
    std::cout << elapsed_time.count() << std::endl;
    return 0;
    /*
        end = std::chrono::steady_clock::now();
        double solve_time = elapsed_time_1.count() + elapsed_time_3.count() + elapsed_time_4.count() + elapsed_time_5.count() + elapsed_time_6.count() + elapsed_time_7.count();

        std::cout << "Mesh generation time: " << elapsed_time_1.count() << "s (" << double(100.0 * elapsed_time_1.count() / solve_time) << "%)" << std::endl;
        std::cout << "Local k matrix computation and global K matrix assembly time: " << elapsed_time_3.count() << "s (" << double(100.0 * elapsed_time_3.count() / solve_time) << "%)" << std::endl;
        std::cout << "Global K matrix reduction time: " << elapsed_time_4.count() << "s (" << double(100.0 * elapsed_time_4.count() / solve_time) << "%)" << std::endl;
        std::cout << "System solve time: " << elapsed_time_5.count() << "s (" << double(100.0 * elapsed_time_5.count() / solve_time) << "%)" << std::endl;
        std::cout << "Original u matrix population time: " << elapsed_time_6.count() << "s (" << double(100.0 * elapsed_time_6.count() / solve_time) << "%)" << std::endl;
        std::cout << "Strain and stress calculation time: " << elapsed_time_7.count() << "s (" << double(100.0 * elapsed_time_7.count() / solve_time) << "%)" << std::endl;
        std::cout << "Overall run time: " << solve_time << "s" << std::endl;
        */
}
