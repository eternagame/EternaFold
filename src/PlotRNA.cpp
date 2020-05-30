/////////////////////////////////////////////////////////////////
// PlotRNA.cpp
//
// Plot an RNA secondary structure.
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
// Header files
/////////////////////////////////////////////////////////////////

#include "gd.h"
#include "gdfontmb.h"
#include "gdfonts.h"
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include "Utilities.hpp"
#include "SStruct.hpp"
#include "EncapsulatedPostScript.hpp"

/////////////////////////////////////////////////////////////////
// Constants
/////////////////////////////////////////////////////////////////

const int HIGH_CONFIDENCE_COLOR[] = { 0, 0, 0 };
const int LOW_CONFIDENCE_COLOR[] = { 216, 216, 255 };

const double SCALE = 12;
const int BORDER = 50;
const int FONT_SIZE = 10;
const double PI = 3.141592653589793238462643383279502884197169399;

/////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////

struct Point
{
    double x, y;
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
};

/////////////////////////////////////////////////////////////////
// Globals
/////////////////////////////////////////////////////////////////

int HEADER = 0;
std::string png_output_filename = "";
std::string eps_output_filename = "";

/////////////////////////////////////////////////////////////////
// ReadCoords()
//
// Read coordinates file.
/////////////////////////////////////////////////////////////////

std::vector<Point> ReadCoords(const std::string &filename)
{
    std::vector<Point> res(1);
    std::ifstream infile(filename.c_str());
    if (infile.fail()) Error("Could not read coordinates file.");
    
    double first, second;
    while (infile >> first >> second)
        res.push_back(Point(first, second));
    
    infile.close();
    return res;  
}

/////////////////////////////////////////////////////////////////
// ReadPosteriors()
//
// Read posteriors.
/////////////////////////////////////////////////////////////////

std::vector<double> ReadPosteriors(const std::string &filename, const SStruct &sstruct)
{
    const std::vector<int> &mapping = sstruct.GetMapping();
    std::vector<double> posteriors(sstruct.GetLength()+1, 1.0);
    if (filename == "") return posteriors;
    
    std::ifstream infile(filename.c_str());
    if (infile.fail()) Error("Could not read posteriors file.");
    
    std::string s;
    int length = 0;
    while (std::getline(infile, s))
    {
        std::istringstream iss(s);
        
        char let;
        int from;
        if (iss >> from >> let)
        {
            length++;
            if (from != length) Error("Bad nucleotide numbering in posteriors file.");
            if (!isalpha(let)) Error("Bad letter in posteriors file.");
            while (iss >> s)
            {
                std::string::size_type idx = s.find(':');
                if (idx == std::string::npos) Error("Badly formatted line in posteriors file.");
                int to = atoi(s.substr(0, idx).c_str());
                if (to < 0) Error("Negative mapping indices not allowed in posteriors file.");
                double value = atof(s.substr(idx + 1).c_str());
                if (value < -0.10 || value > 1.10) Error("Invalid value in posteriors file.");
                value = Clip(value, 0.0, 1.0);
                
                if (from < 1 || from > sstruct.GetLength()) Error("Index in posteriors file does not match BPSEQ length.");
                if (to < 1 || to > sstruct.GetLength()) Error("Index in posteriors file does not match BPSEQ length.");
                if (mapping[from] == to) posteriors[from] = value;
                if (mapping[to] == from) posteriors[to] = value;
                if (mapping[from] == 0) posteriors[from] -= value;
                if (mapping[to] == 0) posteriors[to] -= value;
            }
        }
        else
        {
            Error("Badly formatted line in posteriors file.");
        }
    }
    
    infile.close();
    
    return posteriors;
}

/////////////////////////////////////////////////////////////////
// DrawRNA()
//
// Draw RNA secondary structure.
/////////////////////////////////////////////////////////////////

void DrawRNA(const SStruct &sstruct,
             std::vector<Point> coords,
             const std::vector<double> &posteriors,
             std::string title)
{
    bool toggle_png = (png_output_filename != std::string(""));
    bool toggle_eps = (eps_output_filename != std::string(""));
        
    const std::string &sequence = sstruct.GetSequences()[0];
    const std::vector<int> &mapping = sstruct.GetMapping();
    
    // transform and recenter coordinates
    
    for (size_t i = 1; i < coords.size(); i++)
        coords[i].y *= -1;
    
    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    
    for (size_t i = 1; i < coords.size(); i++)
    {
        min_x = std::min(min_x, coords[i].x);
        min_y = std::min(min_y, coords[i].y);
        max_x = std::max(max_x, coords[i].x);
        max_y = std::max(max_y, coords[i].y);
    }
    
    int width = 0;
    int height = 0;
    for (size_t i = 1; i < coords.size(); i++)
    {
        coords[i].x = std::floor((coords[i].x - min_x) * SCALE);
        coords[i].y = std::floor((coords[i].y - min_y) * SCALE);
        width = std::max(width, int(coords[i].x)+1);
        height = std::max(height, int(coords[i].y)+1);
    }
    
    // create image
    
    int image_width = width + 2*BORDER;
    int image_height = height + 2*BORDER + HEADER;
    
    gdImagePtr image = NULL;
    EncapsulatedPostScript *eps = NULL;
    std::ofstream out;    
    
    if (toggle_eps)
    {
        out.open(eps_output_filename.c_str());
        eps = new EncapsulatedPostScript(out, image_width, image_height, FONT_SIZE);
    }
    if (toggle_png)
    {
        image = gdImageCreateTrueColor(image_width, image_height);
    }
    
    // allocate palette

    int white = 0;
    int black = 0;
    std::vector<int> confidence;

    if (toggle_png)
    {
        white = gdImageColorAllocate(image, 255, 255, 255);
        black = gdImageColorAllocate(image, 0, 0, 0);
        
        confidence.resize(101);
        for (int i = 0; i <= 100; i++)
        {
            int r = (LOW_CONFIDENCE_COLOR[0] * (100 - i) + HIGH_CONFIDENCE_COLOR[0] * i) / 100;
            int g = (LOW_CONFIDENCE_COLOR[1] * (100 - i) + HIGH_CONFIDENCE_COLOR[1] * i) / 100;
            int b = (LOW_CONFIDENCE_COLOR[2] * (100 - i) + HIGH_CONFIDENCE_COLOR[2] * i) / 100;
            confidence[i] = gdImageColorAllocate(image, r, g, b);
        }
    }
    
    // set white background

    if (toggle_png)
    {
        gdImageFilledRectangle(image, 0, 0, image_width-1, image_height-1, white);
    }
    
    // obtain font

    gdFontPtr text_font = NULL;
    gdFontPtr number_font = NULL;
    
    if (toggle_png)
    {
        text_font = gdFontGetMediumBold();
        number_font = gdFontGetSmall();
    }
        
    // draw RNA
    
    for (size_t i = 1; i < coords.size(); i++)
    {
        
        // choose color
        
        double p = posteriors[i];
        int index = Clip(int(p * 100), 0, 100);
        
        if (toggle_eps)
        {
            eps->SetRGBColor(double((LOW_CONFIDENCE_COLOR[0] * (100 - index) + HIGH_CONFIDENCE_COLOR[0] * index) / 100) / 256.0,
                             double((LOW_CONFIDENCE_COLOR[1] * (100 - index) + HIGH_CONFIDENCE_COLOR[1] * index) / 100) / 256.0,
                             double((LOW_CONFIDENCE_COLOR[2] * (100 - index) + HIGH_CONFIDENCE_COLOR[2] * index) / 100) / 256.0);
        }

        // draw letter
        
        int x = BORDER + int(coords[i].x);
        int y = HEADER + BORDER + int(coords[i].y);

        if (toggle_eps)
        {
            eps->DrawString(x, y, SPrintF("%c", sequence[i]), 1, 1);
        }

        if (toggle_png)
        {
            gdImageSetAntiAliased(image, confidence[index]);
            gdImageChar(image, text_font, x - text_font->w/2, y - text_font->h/2, sequence[i], confidence[index]);
        }
        
        // draw bonds
        
        if (mapping[i] > int(i))
        {
            int ox = BORDER + int(coords[mapping[i]].x);
            int oy = HEADER + BORDER + int(coords[mapping[i]].y);
            
            char left = toupper(sequence[i]);
            char right = toupper(sequence[mapping[i]]);
            
            // draw line for a Watson-Crick pair
            
            if (left == 'A' && (right == 'T' || right == 'U') ||
                left == 'C' && right == 'G' ||
                left == 'G' && right == 'C' ||
                (left == 'T' || left == 'U') && right == 'A')
            {

                if (toggle_eps)
                {
                    eps->DrawLine(0.75*x + 0.25*ox, 0.75*y + 0.25*oy,
                                  0.25*x + 0.75*ox, 0.25*y + 0.75*oy, 2.0);
                }

                if (toggle_png)
                {
                    gdImageLine(image, 
                                int(0.75*x + 0.25*ox), int(0.75*y + 0.25*oy),
                                int(0.25*x + 0.75*ox), int(0.25*y + 0.75*oy),
                                gdAntiAliased);
                }                
            } 
            
            // draw dot for GU pair
            
            else if (left == 'G' && (right == 'T' || right == 'U') ||
                     (left == 'T' || left == 'U') && right == 'G')
            {

                if (toggle_eps)
                {
                    eps->DrawCircle(0.5*x + 0.5*ox, 0.5*y + 0.5*oy, 3.0);
                }

                if (toggle_png)
                {
                    gdImageFilledArc(image,
                                     (int)(0.5*x + 0.5*ox), (int)(0.5*y + 0.5*oy),
                                     5, 5, 0, 360, gdAntiAliased, gdArc);
                }
            }

            // draw non-canonical base-pairings
            
            else
            {
                if (toggle_eps)
                {
                    eps->DrawLine(0.75*x + 0.25*ox, 0.75*y + 0.25*oy,
                                  0.25*x + 0.75*ox, 0.25*y + 0.75*oy, 2.0);
                }
            }
        }

        if (toggle_eps)
        {
            eps->SetColorBlack();
        }
        
        // draw number if needed
        
        if (i % 10 == 0)
        {
            std::string number = SPrintF ("%d", i);
            int lx = BORDER + int(coords[i-1].x) - x;
            int ly = HEADER + BORDER + int(coords[i-1].y) - y;
            int rx = (i + 1 < coords.size()) ? BORDER + int(coords[i+1].x) - x : -lx;
            int ry = (i + 1 < coords.size()) ? HEADER + BORDER + int(coords[i+1].y) - y : -ly;
            double al = std::atan2(double(ly), double(lx));
            double ar = std::atan2(double(ry), double(rx));
            if (al > ar) al += 2 * PI;
            double an = (al + ar) / 2.0;
            int nx = int(1.5 * SCALE * std::cos(an));
            int ny = int(1.5 * SCALE * std::sin(an));

            if (toggle_eps)
            {
                eps->DrawString(x + nx, y + ny, number, 1, 1);
            }

            if (toggle_png)
            {
                gdImageString(image, number_font, 
                              x + nx - int(number.length()) * number_font->w/2, 
                              y + ny - number_font->h/2, 
                              reinterpret_cast<unsigned char *>(const_cast<char *>(number.c_str())), black);
            }
        }
    }
    
    if (title != "")
    {
        
        // write title
        
        if (toggle_eps)
        {
            eps->DrawString(image_width / 2, HEADER/2, title, 1, 1);
        }

        if (toggle_png)
        {
            int max_chars = int(0.75 * image_width / number_font->w);
            if (int(title.length()) > max_chars)
            {
                title = title.substr(0, std::max(0, max_chars - 3)) + "...";
            }            

            gdImageString(image, number_font, image_width / 2 - int(title.length()) * number_font->w/2,
                          HEADER/2, reinterpret_cast<unsigned char *>(const_cast<char *>(title.c_str())), black);
        }
            
    }

    if (toggle_eps)
    {
        eps->Close();
        out.close();
    }

    if (toggle_png)
    {
        // write to output file
        
        FILE *pngout = fopen(png_output_filename.c_str(), "wb");
        gdImagePng(image, pngout);
        fclose(pngout);
        gdImageDestroy(image);
    }
}

/////////////////////////////////////////////////////////////////
// main()
//
// Main program.
/////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if (argc < 4)
    {
        std::cerr << std::endl
                  << "Usage: " << argv[0] << " [OPTION]... BPSEQFILE COORDSFILE" << std::endl
                  << std::endl
                  << "       where BPSEQFILE           is the name of the input BPSEQ file" << std::endl
                  << "             COORDFILE           is the name of the input coordinates file" << std::endl
                  << std::endl
                  << "Miscellaneous arguments:" << std::endl
                  << "  --posteriors POSTERIORSFILE    is an optional posteriors file" << std::endl
                  << "  --title TITLE                  is an optional title" << std::endl
                  << "  --eps FILENAME                 specifies EPS format output" << std::endl
                  << "  --png FILENAME                 specifies PNG format output" << std::endl
                  << std::endl;
        exit(1);
    }
    
    // parse arguments
    
    std::string posteriors_filename;
    std::string title;

    std::vector<std::string> default_args;
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (std::string(argv[i]) == "--posteriors")
            {
                if (i == argc-1) Error("Argument required after --posteriors.");
                posteriors_filename = argv[++i];
            }
            else if (std::string(argv[i]) == "--title")
            {
                if (i == argc-1) Error("Argument required after --title.");
                title = argv[++i];
                HEADER = 50;
            }
            else if (std::string(argv[i]) == "--eps")
            {
                if (i == argc-1) Error("Argument required after --eps.");
                eps_output_filename = argv[++i];
            }
            else if (std::string(argv[i]) == "--png")
            {
                if (i == argc-1) Error("Argument required after --png.");
                png_output_filename = argv[++i];
            }
            else
            {
                std::cerr << "Unknown argument: " << argv[i] << std::endl;
                exit (1);
            }
        }
        else
        {
            default_args.push_back(argv[i]);
        }
    }
    
    if (default_args.size() != 2) Error("Incorrect number of arguments.");
    if (eps_output_filename == std::string("") &&
        png_output_filename == std::string(""))
        Error("At least one output format (eps, png) must be specified.");
    
    SStruct sstruct(default_args[0]);
    std::vector<Point> coords = ReadCoords(default_args[1]);
    std::vector<double> posteriors = ReadPosteriors(posteriors_filename, sstruct);
    
    DrawRNA(sstruct, coords, posteriors, title);
}
