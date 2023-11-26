#ifndef __STATISTICS_H
#define __STATISTICS_H

#include<vector>

class Statistics
{
    public:

    Statistics() {};

    template<typename M>
    Statistics (const M& data) {compute (data);};

    template<typename M, typename F>
    Statistics (const M& data, const F& f) {compute (data, f);};
        
    struct Info {
        double mean;
        double var;
    };

    std::vector<Info> statvec;

    template<typename M>
    void compute (const M&);

    template<typename M, typename F>
    void compute (const M&, const F&);

    template<typename S, typename T>
    void print (const S&, T&);

    template<typename S>
    void print (const S& s) {print (s, std::cout);};
};

template <typename M>
void Statistics::compute (const M& m)
{
    int nensemble = m.size();
    int nsteps = m[0].size();
    statvec.clear();
    statvec.assign(nsteps, (Info) {0,0});
    int i, j;
    for (i = 0; i < nensemble; i++)
        for (j = 0; j < nsteps; j++) {
            statvec[j].mean += m[i][j];
            statvec[j].var += m[i][j] * m[i][j];
        }
    for (j = 0; j < nsteps; j++) {
        statvec[j].mean /= (double) nensemble;
        statvec[j].var -= nensemble * statvec[j].mean * statvec[j].mean;
        statvec[j].var /= (double) nensemble - 1;
    }
}

template <typename M, typename F>
void Statistics::compute (const M& m, const F& f)
{
    int nensemble = m.size();
    int nsteps = m[0].size();
    statvec.clear();
    statvec.assign(nsteps, (Info) {0,0});
    int i, j;
    for (i = 0; i < nensemble; i++)
        for (j = 0; j < nsteps; j++) {
            double x = f (m[i][j]);
            statvec[j].mean += x;
            statvec[j].var += x * x;
        }
    for (j = 0; j < nsteps; j++) {
        statvec[j].mean /= (double) nensemble;
        statvec[j].var -= nensemble * statvec[j].mean * statvec[j].mean;
        statvec[j].var /= (double) nensemble - 1;
    }
}

template <typename S, typename T>
void Statistics::print (const S& steps, T& file)
{
    file << std::setw(25) << "step";
    file << std::setw(25) << "mean";
    file << std::setw(25) << "var";
    file << std::endl;
    for (unsigned i = 0; i < steps.nsteps; i++)
    {
        file << std::setw(25) << std::setprecision(8) << steps.cum_step(i);
        file << std::setw(25) << std::setprecision(8) << statvec[i].mean;
        file << std::setw(25) << std::setprecision(8) << statvec[i].var;
        file << std::endl;
    } 
}

#endif