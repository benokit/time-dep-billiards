#ifndef __ROUTINES_H
#define __ROUTINES_H

#include <cstdio>
#include <vector>
#include <array>
#include <numeric>
#include <functional>
#include <cmath>
#include <fstream>
#include <iostream>

struct HistBeam {
    HistBeam() : count {0} {}
    HistBeam(int a) : count {a} {}
    void increaseCount() {++count;}
    double middle;
    double width;
    double probability;
    double density;
    int count;
};

struct Histogram {

    Histogram() {}

    template<typename T>
    Histogram (const T& data, int n) {computeHistogram(data, n);}

    std::vector<HistBeam> beam;

    int ndata;
    double mean;
    double var;
    double min;
    double max;

    void print(std::ofstream*);
    void print();

    template<typename T>
    void computeHistogram(const T&, int);
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
    double bwidth = range / nbeams;
    double hmin = mean - 0.5 * range;
    hmin = (hmin < min) ? min : hmin;
    double hmax = mean + 0.5 * range;
    hmax = (hmax > max) ? max : hmax;

    for (int i = 0; i < ndata; ++i) {
        double x = data[i];
        if (x > hmin && x < hmax) {
            size_t j = floor ((x - hmin) / bwidth);
            beam[j].increaseCount();
        }
    }
    for (int i = 0; i < nbeams; ++i) {
        beam[i].middle = hmin + (i + 0.5) * bwidth;
        beam[i].width = bwidth;
        beam[i].probability = ((double) beam[i].count) / ((double) ndata);
        beam[i].density = beam[i].probability / bwidth;
    }
}

void Histogram::print(std::ofstream* file)
{
    char buffer[200];
    sprintf(buffer, "%25s%25s\n", "middle", "density");
    *file << buffer;
    for (auto b : beam) {
        sprintf(buffer, "%25.8f%25.8f\n", b.middle, b.density);
        *file << buffer;
    }
}

void Histogram::print()
{
    printf("%25s%25s\n", "middle", "density");
    for (auto b : beam) {printf("%25.8f%25.8f\n", b.middle, b.density);}
}

#endif
