/*
  Copyright (C) 2014 Benjamin Batistic

  This file is part of Billiards Numerical Library.

  Billiards Numerical Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Billiards Numerical Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Billiards Numerical Library.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __HISTOGRAM_H
#define __HISTOGRAM_H

#include <vector>
#include <array>
#include <numeric>
#include <functional>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>


class Histogram {

    public: 

    struct HistBeam;

    Histogram() {}

    template<typename T>
    Histogram (const T& data, int n) {computeHistogram(data, n);}

    std::vector<HistBeam> beam;

    int ndata;
    double mean;
    double var;
    double min;
    double max;

    template<typename T>
    void print(T&);

    void print() {print (std::cout);};

    template<typename T>
    void computeHistogram(const T&, int);

};

struct Histogram::HistBeam {
    HistBeam() : count {0} {}
    HistBeam(int a) : count {a} {}
    void increaseCount() {++count;}
    double mid;
    double width;
    double probability;
    double density;
    int count;
};

template<typename T>
void Histogram::computeHistogram (const T& data, int nbeams)
{
    beam.clear();
    beam.assign(nbeams, 0);
    ndata = data.size();
    mean = accumulate(data.begin(), data.end(), 0.0) / ndata;
    double q = accumulate(data.begin(), data.end(), 0.0, 
            [](double acc, double x) -> double {return acc + x * x;});
    var = (q / ndata) - (mean * mean);
    min = *std::min_element(data.begin(), data.end());
    max = *std::max_element(data.begin(), data.end());
    double range = 10.0 * sqrt (var);
    double hmin = mean - 0.5 * range;
    hmin = (hmin < min) ? min : hmin;
    double hmax = mean + 0.5 * range;
    hmax = (hmax > max) ? max : hmax;
    range = hmax - hmin;
    double bwidth = range / nbeams;

    for (int i = 0; i < ndata; ++i) {
        double x = data[i];
        if (x > hmin && x < hmax) {
            size_t j = floor ((x - hmin) / bwidth);
            beam[j].increaseCount();
        }
    }
    for (int i = 0; i < nbeams; ++i) {
        beam[i].mid = hmin + (i + 0.5) * bwidth;
        beam[i].width = bwidth;
        beam[i].probability = ((double) beam[i].count) / ((double) ndata);
        beam[i].density = beam[i].probability / bwidth;
    }
}

template<typename T>
void Histogram::print(T& file)
{
    file << std::setw(25) << "mids";
    file << std::setw(25) << "density";
    file << std::endl;
    for (auto b : beam) {
        file << std::setw(25) << std::setprecision(8) << b.mid;
        file << std::setw(25) << std::setprecision(8) << b.density;
        file << std::endl;
    }
}

#endif
