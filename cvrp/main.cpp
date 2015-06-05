#pragma once

#include <stdio.h>
#include <assert.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <iterator>

using namespace std;
const int cMax = 999;

class Solver
{
public:
    Solver(vector <vector <int> >& matr) : m(matr)
    {
        const int n = m.size();
        vRow.assign(n, true);
        vCol.assign(n, true);
        d1.assign(n, 0);
        d2.assign(n, 0);
    }
    int solve(vector <int>& res, int hStop);
private:
    bool isSkipRow(int i) const  { return !vRow[i]; }
    bool isSkipCol(int j) const  { return !vCol[j]; }
    void setSkipRow(int i, bool bSkip) { vRow[i] = !bSkip; }
    void setSkipCol(int j, bool bSkip) { vCol[j] = !bSkip; }

    void findEdge(int& iT, int& jT);
    void reduce();
    void calcdij(bool bSkipZero);
private:
    vector <vector <int> > m;
    vector <bool> vRow;
    vector <bool> vCol;
    std::vector <int> d1;
    std::vector <int> d2;
};

void Solver::findEdge(int& iT, int& jT)
{
    int ijMax = -1;
    const int n = m.size();
    for(int i = 0; i < n; ++i)
    {
        if(isSkipRow(i))
            continue;
        for(int j = 0; j < n; ++j)
        {
            if(isSkipCol(j))
                continue;
            if(m[i][j] == 0 && ijMax < d1[i] +d2[j])
            {
                iT = i;
                jT = j;
                ijMax = d1[i] +d2[j];
            }
        }
    }
}

void Solver::reduce()
{
    const int n = m.size();
    for(int i = 0; i < n; ++i)
    {
        if(isSkipRow(i))
            continue;
        for(int j = 0; j < n; ++j)
        {
            if(isSkipCol(j))
                continue;
            m[i][j] -= d1[i] + d2[j];
            if(m[i][j] < 0)
                m[i][j] = 0;
        }
    }
}

void Solver::calcdij(bool bSkipZero)
{
    const int n = m.size();
    // update d1 d2
    for(int i = 0; i < n; ++i)
    {
        int min = cMax;
        int nZero = 0;
        if(isSkipRow(i))
            continue;
        for(int j = 0; j < n; ++j)
        {
            if(isSkipCol(j))
                continue;
            const int val = m[i][j];
            if(bSkipZero && val == 0)
            {
                ++nZero;
                if(nZero == 1)
                    continue;
            }
            if(val < min)
                min = val;
        }
        d1[i] = min;
    }

    for(int j = 0; j < n; ++j)
    {
        int min = cMax;
        int nZero = 0;
        if(isSkipCol(j))
            continue;
        for(int i = 0; i < n; ++i)
        {
            if(isSkipRow(i))
                continue;
            int val = m[i][j];
            if(!bSkipZero)
                val -= d1[i];
            if(bSkipZero && val == 0)
            {
                ++nZero;
                if(nZero == 1)
                    continue;
            }
            if(val < min)
                min = val;
        }
        d2[j] = min;
    }
}

int Solver::solve(vector <int>& res, int hStop)
{
    res.clear();
    const int n = m.size();
    calcdij(false);
    int h = 0;
    for(int i = 0; i < n; ++i)
        h += d1[i]+d2[i];
    std::vector <std::pair<int, int> > edges;
    // reduction
    while(true)
    {
        if(h >= hStop)
            return cMax;
        int nSkip = 0;
        for(int i = 0; i < n; ++i)
            nSkip += isSkipRow(i);
        if(nSkip == n)
            break;
        reduce();
        calcdij(true);
        int iMax, jMax;
        findEdge(iMax, jMax);
        // exclude edge
        int hExclude = h + d1[iMax] + d2[jMax];
        // include edge
        int hInclude = h;
        {
            int tmpIj = m[jMax][iMax];
            m[jMax][iMax] = cMax;
            setSkipRow(iMax, true);
            setSkipCol(jMax, true);
            calcdij(false);
            for(int i = 0; i < n; ++i)
            {
                if(!isSkipRow(i))
                    hInclude += d1[i];
                if(!isSkipCol(i))
                    hInclude += d2[i];
            }
            // restore
            setSkipRow(iMax, false);
            setSkipCol(jMax, false);
            m[jMax][iMax] = tmpIj;
        }
       
        // compare branches
        if(hInclude <= hExclude)
        {
            edges.push_back(std::pair<int, int>(iMax, jMax));
            m[jMax][iMax] = cMax; // skip edge j i
            setSkipRow(iMax, true);
            setSkipCol(jMax, true);
            h = hInclude;
        }
        else
        {
            m[iMax][jMax] = cMax; // skip edge i j 
            h = hExclude;
        }
        calcdij(false);
    }
    std::vector<int> vtx;
    // extract the lists of vertices from the list of edges
    res.push_back(0);
    while(true)
    {
        for(int j = 0; j < edges.size(); ++j)
        {
            if(res.back() == edges[j].first)
            {
                res.push_back(edges[j].second);
                // mark edge as processed
                edges[j].first  = n+1;
                edges[j].second = n+1;
                break;
            }
        }
        if(res.front() == res.back())
            break;
    }
    return h;
}

int calcVal(const vector <vector <int> >& m, vector <int>& vtx)
{
    int h = 0;
    for(int i = 1; i < vtx.size(); ++i)
        h += m[vtx[i-1]][vtx[i]]; // weight of the edge
    return h;
}

void extractCars(const vector <vector <int> >& m, const vector <int>& resCommon, int nCar, vector <vector <int> >& res, vector <int>& weights)
{
    int n = m.size() - nCar + 1;
    vector<int> v;
    v.push_back(resCommon[0]);
    assert(resCommon.front() == 0 && resCommon.back() == 0);
    for(int j = 1; j < resCommon.size(); ++j)
    {
        if(resCommon[j] >= n || resCommon[j] == 0)
        {
            v[0] = 0;
            v.push_back(0);
            res.push_back(v);
            weights.push_back(calcVal(m, v));
            v.clear();
            v.push_back(resCommon[j]);
        }
        else
            v.push_back(resCommon[j]);
    }
}

void printResult(int h, const vector <vector <int> >& res, const vector <int>& weights)
{
    cout << "Value total: " << h << endl;
    for(int j = 0; j < res.size(); ++j)
    {
        cout << "Car: " << j+1 << " value: "<< weights[j] << endl;
        cout << "Path: ";
        copy(res[j].begin(), res[j].end(), ostream_iterator<int>(cout, " "));
        cout << endl;
    }
    cout << endl;
}

int main(int argc, char* argv[])
{
    if(argc != 2)
        return 1;
    int nTown = 0, nCarMax = 0, limit = 0;
    // 1. read data
    const char* fName = argv[1];
    ifstream infile(fName);
    if(!infile.is_open())
        return 1;
    infile >> nTown;
    infile >> nCarMax;
    infile >> limit;
    if(nCarMax > nTown)
        nCarMax = nTown;
    vector <vector <int> > m(nTown, vector<int>(nTown, 0));
    for(int i = 0; i < nTown; ++i)
    {
        for(int j = 0; j < nTown; ++j)
            infile >> m[i][j];
    }
    infile.close();
    // 2. solve
    vector <vector <int> > resOptimal;
    vector <int > weightsOptimal;
    int hMin = cMax;
    for(int nCar = 1; nCar <= nCarMax; ++nCar)
    {
        if(nCar > 1)
        {
            // add new car
            int n = m.size();
            // add column
            for(int i = 0; i < n; ++i)
                m[i].push_back(m[i][0]);
            // add row
            m.push_back(m[0]);
        }
        vector <int> resCommon;
        Solver task(m);
        int h = task.solve(resCommon, hMin);
        if(resCommon.empty())
            continue;
        vector <vector <int> > res;
        vector <int > weights;
        extractCars(m, resCommon, nCar, res, weights);
        printResult(h, res, weights); // use this code for debug
        int hCarMax = 0;
        for(int j = 0; j < weights.size(); ++j)
        {
            int hCar = weights[j];
            if(hCar > hCarMax)
                hCarMax = hCar;
        }
        if(hCarMax > limit)
            continue;
        if(h < hMin)
        {
            hMin = h;
            resOptimal = res;
            weightsOptimal = weights;
        }
    }

    // 3. print results
    if(hMin >= cMax)
        cout << "can not find the solution" << endl;
    else
    {
        cout << "Optimal solution is" << endl;
        printResult(hMin, resOptimal, weightsOptimal);
    }
    return 0;
}