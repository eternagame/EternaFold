/////////////////////////////////////////////////////////////////
// MakeCoords.hpp
/////////////////////////////////////////////////////////////////
#ifndef MAKECOORDS_HPP
#define MAKECOORDS_HPP

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include "LBFGS.hpp"
#include "Utilities.hpp"
#include "SStruct.hpp"

/////////////////////////////////////////////////////////////////
// Constants
/////////////////////////////////////////////////////////////////

const double LOOP_STRENGTH = 30;
const double REPULSIVE_STRENGTH = 0.2;
const double BACKBONE_STRENGTH = 50;
const double PI = 3.141592653589793238462643383279502884197169399;
const int STEM_WIDTH = 3;

/////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////

struct Point
{
    double x, y;
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
};

enum ConstraintType
{
    ConstraintType_LENGTH,
    ConstraintType_REPULSIVE
};

struct Constraint
{
    ConstraintType type;
    int i;
    int j;
    double dist;
    double strength;
    
    Constraint(ConstraintType type, int i, int j, double dist, double strength) : 
        type(type), i(i), j(j), dist(dist), strength(strength) {}
};

/////////////////////////////////////////////////////////////////
// class MakeCoords
/////////////////////////////////////////////////////////////////

class MakeCoords : public LBFGS
{
    const SStruct &sstruct;
    std::vector<Constraint> constraints;
    
public:    
    std::vector<Point> coords;
    std::vector<Point> gradients;
    
    std::vector<double> GetParams(const std::vector<Point> &points) const;
    void SetParams(const std::vector<double> &values);
    void PrintParams(const std::string &filename) const;
    
    std::vector<int> ComputeLoop(const std::vector<int> &mapping, int left) const;
    Point ComputeLoopCenter(Point p1, Point p2, int k, int n) const;
    std::vector<Point> ComputeLoopPositions(Point p1, Point center, int n) const;
    void InitialPlacement();
    
    double Distance(Point p, Point q) const;
    void AddConstraints();
    
    
    MakeCoords(const SStruct &sstruct);
    void ComputeGradient(std::vector<double> &gradient, const std::vector<double> &values);
    double ComputeFunction(const std::vector<double> &values);
    void Report(int iteration, const std::vector<double> &theta, double objective, double step_length);
    void Report(const std::string &s);
};

#endif
