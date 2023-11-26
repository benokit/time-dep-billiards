#ifndef __STATISTICS_H
#define __STATISTICS_H

#include <vector>

class Statistics {
    public:
        Statistics(){};

        template <typename M>
        Statistics(const M &data) { compute(data); };

        template <typename M, typename F>
        Statistics(const M &data, const F &f) { compute(data, f); };

        struct Info
        {
            double mean;
            double var;
        };

        std::vector<Info> statistics_info;

        template <typename M>
        void compute(const M &);

        template <typename M, typename F>
        void compute(const M &, const F &);

        template <typename S, typename T>
        void print(const S &, T &);

        template <typename S>
        void print(const S &s) { print(s, std::cout); };
};

template <typename M>
void Statistics::compute(const M &m)
{
    int n_ensemble = m.size();
    int n_steps = m[0].size();
    statistics_info.clear();
    statistics_info.assign(n_steps, (Info){0, 0});
    int i, j;
    for (i = 0; i < n_ensemble; i++)
        for (j = 0; j < n_steps; j++)
        {
            statistics_info[j].mean += m[i][j];
            statistics_info[j].var += m[i][j] * m[i][j];
        }
    for (j = 0; j < n_steps; j++)
    {
        statistics_info[j].mean /= (double)n_ensemble;
        statistics_info[j].var -= n_ensemble * statistics_info[j].mean * statistics_info[j].mean;
        statistics_info[j].var /= (double)n_ensemble - 1;
    }
}

template <typename M, typename F>
void Statistics::compute(const M &m, const F &f)
{
    int n_ensemble = m.size();
    int n_steps = m[0].size();
    statistics_info.clear();
    statistics_info.assign(n_steps, (Info){0, 0});
    int i, j;
    for (i = 0; i < n_ensemble; i++)
        for (j = 0; j < n_steps; j++)
        {
            double x = f(m[i][j]);
            statistics_info[j].mean += x;
            statistics_info[j].var += x * x;
        }
    for (j = 0; j < n_steps; j++)
    {
        statistics_info[j].mean /= (double)n_ensemble;
        statistics_info[j].var -= n_ensemble * statistics_info[j].mean * statistics_info[j].mean;
        statistics_info[j].var /= (double)n_ensemble - 1;
    }
}

template <typename S, typename T>
void Statistics::print(const S &steps, T &file)
{
    file << std::setw(25) << "step";
    file << std::setw(25) << "mean";
    file << std::setw(25) << "var";
    file << std::endl;
    for (unsigned i = 0; i < steps.n_steps; i++)
    {
        file << std::setw(25) << std::setprecision(8) << steps.cum_step(i);
        file << std::setw(25) << std::setprecision(8) << statistics_info[i].mean;
        file << std::setw(25) << std::setprecision(8) << statistics_info[i].var;
        file << std::endl;
    }
}

#endif