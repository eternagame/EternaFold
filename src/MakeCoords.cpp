/////////////////////////////////////////////////////////////////
// MakeCoords.cpp
/////////////////////////////////////////////////////////////////

#include "MakeCoords.hpp"

/////////////////////////////////////////////////////////////////
// MakeCoords::GetParams()
//
// Retrieve unrolled parameter vector.
/////////////////////////////////////////////////////////////////

std::vector<double> MakeCoords::GetParams(const std::vector<Point> &points) const
{
    std::vector<double> res;
    for (size_t i = 1; i < points.size(); i++)
    {
        res.push_back(points[i].x);
        res.push_back(points[i].y);
    }
    return res;
}

/////////////////////////////////////////////////////////////////
// MakeCoords::SetParams()
//
// Load unrolled parameter vector.
/////////////////////////////////////////////////////////////////

void MakeCoords::SetParams(const std::vector<double> &values)
{
    coords = std::vector<Point>(sstruct.GetLength()+1);
    gradients = std::vector<Point>(sstruct.GetLength()+1);
    
    int k = 0;
    for (size_t i = 1; i < coords.size(); i++)
    {
        coords[i].x = values[k++];
        coords[i].y = values[k++];
    }
}

/////////////////////////////////////////////////////////////////
// MakeCoords::PrintParams()
//
// Print parameters.
/////////////////////////////////////////////////////////////////

void MakeCoords::PrintParams(const std::string &filename) const
{
    std::ofstream outfile(filename.c_str());
    if (!outfile) Error("Could not open coordinates file for writing.");
    for (int i = 1; i <= sstruct.GetLength(); i++)
    {
        outfile << coords[i].x << " " << coords[i].y << std::endl;
    }
    outfile.close();
}

/////////////////////////////////////////////////////////////////
// ComputeLoop()
//
// Compute the residues involved in a loop.
/////////////////////////////////////////////////////////////////

std::vector<int> MakeCoords::ComputeLoop(const std::vector<int> &mapping, int left) const
{
    std::vector<int> ret(1, left);
    int level = 0;
    
    // look for other letters in the loop
    
    for (int k = left+1; k < mapping[left]; k++)
    {
        if (mapping[k] == 0)
        {
            if (level == 0) ret.push_back(k);
        }
        else
        {
            if (mapping[k] < k) level--;
            if (level == 0) ret.push_back(k);
            if (mapping[k] > k)
            {
                if (level == 0) 
                    for (int i = 0; i < STEM_WIDTH-1; i++) 
                        ret.push_back(0);
                level++;
            }
        }
    }  
    
    ret.push_back(mapping[left]);
    for (int i = 0; i < STEM_WIDTH-1; i++) ret.push_back(0);
  
    return ret;
}

/////////////////////////////////////////////////////////////////
// MakeCoords::ComputeLoopCenter()
//
// Given two points (x1,y1) and (x2,y2) that span a central
// angle of 2*PI*k/n (going counterclockwise), compute the
// center (xc,yc).
/////////////////////////////////////////////////////////////////

Point MakeCoords::ComputeLoopCenter(Point p1, Point p2, int k, int n) const
{
    double theta = PI*(0.5 - double(k)/n);
    return Point(p1.x + 0.5 * ((p2.x-p1.x) - std::tan(theta) * (p2.y-p1.y)),
                 p1.y + 0.5 * (std::tan(theta) * (p2.x-p1.x) + (p2.y-p1.y)));
}

/////////////////////////////////////////////////////////////////
// MakeCoords::ComputeLoopPositions()
//
// Given the initial and center point of a loop, and the number
// of residues in the loop, compute the positions of each residue
// in the loop (clockwise).
/////////////////////////////////////////////////////////////////

std::vector<Point> MakeCoords::ComputeLoopPositions(Point p1, Point center, int n) const
{
    std::vector<Point> ret;
    double alpha = 2*PI/n;
    double dx = p1.x - center.x;
    double dy = p1.y - center.y;
    for (int i = 0; i < n; i++)
        ret.push_back(Point(center.x + std::cos(-alpha*i)*dx - std::sin(-alpha*i)*dy,
                            center.y + std::sin(-alpha*i)*dx + std::cos(-alpha*i)*dy));
    return ret;
}

/////////////////////////////////////////////////////////////////
// MakeCoords::InitialPlacement()
//
// Initial placement of all coordinates using deterministic
// procedure.
/////////////////////////////////////////////////////////////////

void MakeCoords::InitialPlacement()
{
    coords = std::vector<Point>(sstruct.GetLength()+1);
    gradients = std::vector<Point>(sstruct.GetLength()+1);
    
    std::vector<int> mapping = sstruct.GetMapping();
    const int L = int(mapping.size()) - 1;
    
    // set any unknown pairings to be unpaired
    
    for (int i = 1; i <= L; i++) mapping[i] = std::max(0, mapping[i]);
    
    // temporarily force a pairing of the first and last residues
    
    if (mapping[1] != 0) mapping[mapping[1]] = 0;
    if (mapping[L] != 0) mapping[mapping[L]] = 0;
    mapping[1] = L;
    mapping[L] = 1;
    
    // set coordinates of first and last residues
    
    coords[1].x = 0;
    coords[1].y = 0;
    coords[L].x = STEM_WIDTH;
    coords[L].y = 0;
    
    // iteratively find internal loops and grow a structure from them
    
    for (int i = 1; i <= L; i++)
    {
        if (mapping[i] > i)
        {
            std::vector<int> loop = ComputeLoop(mapping, i);
            int j = mapping[i];
            Point center = ComputeLoopCenter(coords[i], coords[j], STEM_WIDTH, loop.size());
            std::vector<Point> positions = ComputeLoopPositions(coords[i], center, loop.size());
            for (int k = 0; k < int(loop.size()); k++)
            {
                if (loop[k] > 0) coords[loop[k]] = positions[k];
            }
        }
    }
}

/////////////////////////////////////////////////////////////////
// MakeCoords::Distance()
//
// Compute squared distance between two points.
/////////////////////////////////////////////////////////////////

double MakeCoords::Distance(Point p, Point q) const
{
    p.x -= q.x;
    p.y -= q.y;
    return std::sqrt(p.x*p.x + p.y*p.y);
}

/////////////////////////////////////////////////////////////////
// MakeCoords::AddConstraints()
//
// Add all constraints associated with loops in RNA structure.
/////////////////////////////////////////////////////////////////

void MakeCoords::AddConstraints()
{
    std::vector<int> mapping = sstruct.GetMapping();
    const int L = int(mapping.size()) - 1;
    
    // set any unknown pairings to be unpaired
    
    for (int i = 1; i <= L; i++) mapping[i] = std::max(0, mapping[i]);
    
    // add constraints for loops
    
    for (int i = 1; i <= L; i++)
    {
        if (mapping[i] > i)
        {
            std::vector<int> loop = ComputeLoop (mapping, i);
            for (size_t k1 = 0; k1 < loop.size(); k1++)
            {
                if (loop[k1] <= 0) continue;
                for (size_t k2 = k1+1; k2 < loop.size(); k2++)
                {
                    if (loop[k2] <= 0) continue;
                    double dist = Distance(coords[loop[k1]], coords[loop[k2]]);
                    constraints.push_back(Constraint(ConstraintType_LENGTH, loop[k1], loop[k2], dist, LOOP_STRENGTH / std::max(1, int(loop.size()) - 7)));
                }
            }
        }
    }
        
    // add backbone constraints
    
    for (int i = 1; i < L; i++)
    {
        double dist = Distance(coords[i], coords[i+1]);
        constraints.push_back(Constraint(ConstraintType_LENGTH, i, i+1, dist, BACKBONE_STRENGTH));
    }
    
    // add repulsive constraints
    
    for (int i = 1; i <= L; i++)
        for (int j = i+1; j <= L; j++)
            constraints.push_back(Constraint(ConstraintType_REPULSIVE, i, j, 0, REPULSIVE_STRENGTH));
}

/////////////////////////////////////////////////////////////////
// MakeCoords::MakeCoords()
//
// Constructor.
/////////////////////////////////////////////////////////////////

MakeCoords::MakeCoords(const SStruct &sstruct) : 
    LBFGS<double>(), sstruct(sstruct)
{
    InitialPlacement();
    AddConstraints();
}

/////////////////////////////////////////////////////////////////
// MakeCoords::ComputeGradient()
//
// Compute objective gradient.
/////////////////////////////////////////////////////////////////

void MakeCoords::ComputeGradient(std::vector<double> &gradient, const std::vector<double> &values)
{
    SetParams(values);
  
    for (size_t i = 0; i < constraints.size(); i++)
    {
        
        // precompute points involved in constraint
        
        const Point &p1 = coords[constraints[i].i];
        const Point &p2 = coords[constraints[i].j];
        const double r = Distance(p1, p2);
        
        // check type of constraint
        
        switch (constraints[i].type)
        {
            
            case ConstraintType_LENGTH: 
            {
                double a = constraints[i].strength * constraints[i].dist;
                double b = constraints[i].strength / (2 * constraints[i].dist * constraints[i].dist);
                double dVdr = -a/(r*r) + 2*b*r;
                
                gradients[constraints[i].i].x += dVdr * (p1.x - p2.x) / r;
                gradients[constraints[i].i].y += dVdr * (p1.y - p2.y) / r;
                gradients[constraints[i].j].x += dVdr * (p2.x - p1.x) / r;
                gradients[constraints[i].j].y += dVdr * (p2.y - p1.y) / r;
            }
            break;
            
            case ConstraintType_REPULSIVE: 
            {
                double a = constraints[i].strength;
                double dVdr = -a/(r*r);
                
                gradients[constraints[i].i].x += dVdr * (p1.x - p2.x) / r;
                gradients[constraints[i].i].y += dVdr * (p1.y - p2.y) / r;
                gradients[constraints[i].j].x += dVdr * (p2.x - p1.x) / r;
                gradients[constraints[i].j].y += dVdr * (p2.y - p1.y) / r;
            }
            break;
            
            default:
                Error("Invalid constraint type.");
        }
    }
    
    gradient = GetParams(gradients);
}

/////////////////////////////////////////////////////////////////
// MakeCoords::ComputeFunction()
//
// Compute objective function.
/////////////////////////////////////////////////////////////////

double MakeCoords::ComputeFunction(const std::vector<double> &values)
{
    SetParams(values);
    
    double ret = 0;
    for (size_t i = 0; i < constraints.size(); i++)
    {
        // precompute points involved in constraint
        
        const Point &p1 = coords[constraints[i].i];
        const Point &p2 = coords[constraints[i].j];
        const double r = Distance(p1, p2);
        const double d = constraints[i].dist;
        
        // check type of constraint
        
        switch (constraints[i].type)
        {
            
            case ConstraintType_LENGTH: 
            {
                double a = constraints[i].strength * constraints[i].dist;
                double b = constraints[i].strength / (2 * constraints[i].dist * constraints[i].dist);
                double V = a/r + b*r*r - (a/d + b*d*d);
                
                ret += V;
            }
            break;
            
            case ConstraintType_REPULSIVE: 
            {
                double a = constraints[i].strength;
                double V = a/r;
                
                ret += V;
            }
            break;
            
            default:
                Error("Invalid constraint type.");
        }
    }
    
    return ret;
}

/////////////////////////////////////////////////////////////////
// MakeCoords::Report()
//
// Print feedback from minimizer.
/////////////////////////////////////////////////////////////////

void MakeCoords::Report(int iteration, const std::vector<double> &theta, double objective, double step_length)
{
    std::cerr << "Minimizer: iteration " << iteration 
              << ", f = " << objective 
              << ", step length = " << step_length << std::endl;
}

void MakeCoords::Report(const std::string &s)
{
    std::cerr << "Minimizer: " << s << std::endl;
}  

/////////////////////////////////////////////////////////////////
// main()
//
// Main program.
/////////////////////////////////////////////////////////////////
 
int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::cerr << std::endl
                  << "Usage: " << argv[0] << " BPSEQFILE OUTFILE" << std::endl
                  << std::endl
                  << "       where BPSEQFILE    is the name of the input BPSEQ file" << std::endl
                  << "             OUTFILE      is the name of the output coordinates file" << std::endl
                  << std::endl;
        exit (1);
    }
    
    SStruct sstruct(argv[1]);
    MakeCoords folder(sstruct);
    std::vector<double> theta = folder.GetParams(folder.coords);
    folder.Minimize(theta);
    folder.SetParams(theta);
    folder.PrintParams(argv[2]);
}
