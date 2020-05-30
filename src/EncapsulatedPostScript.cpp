//////////////////////////////////////////////////////////////////////
// EncapsulatedPostScript.cpp
//
// The routines shown here for creating encapsulated postscript
// figures were adapted from a modification of PlotRNA by
// Marc Parisien.
//////////////////////////////////////////////////////////////////////

#include "EncapsulatedPostScript.hpp"

//////////////////////////////////////////////////////////////////////
// EncapsulatedPostScript::EncapsulatedPostScript()
//
// Constructor.  Write PostScript prolog.
//////////////////////////////////////////////////////////////////////

EncapsulatedPostScript::EncapsulatedPostScript(std::ostream &out, double image_width, double image_height, int font_size) :
    out(out), image_width(image_width), image_height(image_height), font_size(font_size), done(false)
{
    out << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl
        << "%%BoundingBox: 0 0 " << int(image_width) << " " << int(image_height) << std::endl
        << "1.0 1.0 scale" << std::endl
        << "0 0 translate" << std::endl
        << "/HelveticaBold findfont" << std::endl
        << font_size << " scalefont" << std::endl
        << "setfont" << std::endl;        
}

//////////////////////////////////////////////////////////////////////
// EncapsulatedPostScript::SetRGBColor()
//
// Set current color explicitly.
//////////////////////////////////////////////////////////////////////

void EncapsulatedPostScript::SetRGBColor(double r, double g, double b)
{
    Assert(0.0 <= r && r <= 1.0, "Out-of-range.");
    Assert(0.0 <= g && g <= 1.0, "Out-of-range.");
    Assert(0.0 <= b && b <= 1.0, "Out-of-range.");
    
    out << std::setprecision(3) << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << r << " " << g << " " << b << " setrgbcolor" << std::endl;
}

//////////////////////////////////////////////////////////////////////
// EncapsulatedPostScript::SetColorBlack()
//
// Set current color back to black.
//////////////////////////////////////////////////////////////////////

void EncapsulatedPostScript::SetColorBlack()
{
    out << "0 0 0 setrgbcolor" << std::endl;
}

//////////////////////////////////////////////////////////////////////
// EncapsulatedPostScript::DrawString()
//
// Write a text string.  Adapted from:
// 
// http://www.nipr.ac.jp/~uap-mon/uapm/src.bak/pltSyowaMag_save.c
//
// Text alignment:
//   pos_x : x-align  0:left   1:center 2:right
//   pos_y : y-align  0:bottom 1:center 2:top
//////////////////////////////////////////////////////////////////////

void EncapsulatedPostScript::DrawString(double x, double y, const std::string &s, int pos_x, int pos_y)
{
    if (done) Error("EPS file already closed.");

    int kx = 0, ky = 0;

    out << std::setprecision(3) << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << x << " " << image_height - y << " moveto" << std::endl
        << "(" << s << ")" << std::endl;

    if (pos_x == 1) kx = 2;
    if (pos_x == 2) kx = 1;
    if (pos_y == 1) ky = 2;
    if (pos_y == 2) ky = 1;

    if (pos_x == 1 || pos_x == 2)
    {
        out << "dup stringwidth pop " << kx << " div neg 0 rmoveto" << std::endl;
    }
    
    if( pos_y == 1 || pos_y == 2 )
    {
        out << "gsave" << std::endl
            << "newpath" << std::endl
            << "0 0 moveto" << std::endl
            << "(" << s << ") true charpath flattenpath" << std::endl
            << "pathbbox /charheight exch def pop pop pop" << std::endl
            << "closepath" << std::endl
            << "grestore" << std::endl
            << "0 charheight " << ky << " div neg rmoveto" << std::endl;
    }

    out << "show" << std::endl;
}

//////////////////////////////////////////////////////////////////////
// EncapsulatedPostScript::DrawLine()
//
// Draw a line from (sx,sy) to (ex,ey) with given width.
//////////////////////////////////////////////////////////////////////

void EncapsulatedPostScript::DrawLine(double sx, double sy, double ex, double ey, double width)
{
    if (done) Error("EPS file already closed.");
    out << std::setprecision(3) << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << width << " setlinewidth" << std::endl
        << sx << " " << image_height - sy << " moveto " << ex << " " << image_height - ey
        << " lineto stroke" << std::endl;
}

//////////////////////////////////////////////////////////////////////
// EncapsulatedPostScript::DrawCircle()
//
// Draw a circle at (x,y) with given radius.
//////////////////////////////////////////////////////////////////////

void EncapsulatedPostScript::DrawCircle(double x, double y, double r)
{
    if (done) Error("EPS file already closed.");
    out << std::setprecision(3) << std::setiosflags(std::ios::showpoint | std::ios::fixed)
        << x << " " << image_height - y << " moveto" << std::endl
        << x << " " << image_height - y << " " << r << " 0 360 arc closepath fill" << std::endl;
}

//////////////////////////////////////////////////////////////////////
// EncapsulatedPostScript::Close()
//
// Finish EPS file.
//////////////////////////////////////////////////////////////////////

void EncapsulatedPostScript::Close()
{
    if (done) Error("EPS file already closed.");

    out << "showpage" << std::endl
        << "%EOF" << std::endl;
    done = true;
}
