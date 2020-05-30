//////////////////////////////////////////////////////////////////////
// EncapsulatedPostScript.hpp
//////////////////////////////////////////////////////////////////////

#ifndef ENCAPSULATEDPOSTSCRIPT_HPP
#define ENCAPSULATEDPOSTSCRIPT_HPP

#include <fstream>
#include <string>
#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// class EncapsulatedPostScript
//////////////////////////////////////////////////////////////////////

class EncapsulatedPostScript
{
    std::ostream &out;
    double image_width;
    double image_height;
    int font_size;
    bool done;
    
public:
    EncapsulatedPostScript(std::ostream &out, double image_width, double image_height, int font_size);
    void SetRGBColor(double r, double g, double b);
    void SetColorBlack();
    void DrawString(double x, double y, const std::string &s, int pos_x, int pos_y);
    void DrawLine(double sx, double sy, double ex, double ey, double width);
    void DrawCircle(double x, double y, double r);
    void Close();
};

#endif
