//
// Created by Colin Small on 8/2/21.
//

#include <Rcpp.h>
using namespace Rcpp;

#include "_jenks.h"
#include <vector>
#include <algorithm>

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> jenksInitMatrices(std::vector<double> data, int nClasses){
    int k = nClasses + 1;
    std::vector<std::vector<double>> lowerClassLimits(data.size(), std::vector<double>(k, 0.0));
    std::vector<std::vector<double>> varianceCombinations(data.size(), std::vector<double>(k, std::numeric_limits<double>::infinity()));
    varianceCombinations[0] = std::vector<double>(k, 0.0);

    for(int i = 1; i < nClasses+1; i++){
        lowerClassLimits[1][i] = 1;
    }

    return std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>{lowerClassLimits, varianceCombinations};
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> jenksMatrices(std::vector<double> data, int nClasses){
    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> results = jenksInitMatrices(data, nClasses);
    std::vector<std::vector<double>> lowerClassLimits = std::get<0>(results);
    std::vector<std::vector<double>> varianceCombinations = std::get<1>(results);

    double variance = 0.0;

    for(int l = 1; l < data.size(); l++){
        double sum = 0.0;
        double sumSquares = 0.0;
        double w = 0.0;

        for(int m = 0; m < l; m++){
            int lowerClassLimit = l-m+1;
            double val = data[lowerClassLimit-1];
            w += 1;
            sum += val;
            sumSquares += val * val;
            variance = sumSquares - (sum*sum) / w;
            int i4 = lowerClassLimit-1;
            if(i4 != 0)
            {
                for(int j = 2; j < nClasses+1; j++)
                {
                    if(varianceCombinations[l][j] >= (variance + varianceCombinations[i4][j-1])){
                        lowerClassLimits[l][j] = lowerClassLimit;
                        varianceCombinations[l][j] = variance + varianceCombinations[i4][j-1];
                    }
                }
            }
        }
        lowerClassLimits[l][1] = 1.0;
        varianceCombinations[l][1] = variance;
    }
    return std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>{lowerClassLimits, varianceCombinations};
}

std::vector<double> getJenksBreaks(std::vector<double> data, std::vector<std::vector<double>> lowerClassLimits, int nClasses){
    int k = data.size() - 1;
    std::vector<double> kClass(nClasses+1, 0.0);
    int countNum = nClasses;
    kClass[0] = data[0];
    kClass[nClasses] = data[k];
    while(countNum > 1){
        int elt = lowerClassLimits[k][countNum]-2;
        kClass[countNum-1] = data[elt];
        k = lowerClassLimits[k][countNum]-1;
        countNum -= 1;
    }
    return kClass;
}

// TODO: preallocate memory for clusters and access with [] instead of pushing back
std::vector<std::vector<double>> separateData(std::vector<double> data, std::vector<double> boundaries){
    std::vector<std::vector<double>> clusters = {};
    std::vector<double> cluster = {};

    boundaries.erase(boundaries.begin());
    int i = 0;

    for(double& d : data){
        if(d <= boundaries[i])
            cluster.push_back(d);
        else{
            clusters.push_back(cluster);
            cluster = {d};
            i += 1;
        }
    }

    clusters.push_back(cluster);
    return clusters;

}

std::tuple<std::vector<std::vector<double>>, std::vector<double>> jenks(std::vector<double> data, int nClasses){
//    assert(nClasses > data.size());
    if (nClasses > data.size())
        std::cout << nClasses << " > " << data.size() << std::endl;
    std::sort(data.begin(), data.end());
    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> results = jenksMatrices(data, nClasses);
    std::vector<std::vector<double>> lowerClassLimits = std::get<0>(results);
    std::vector<std::vector<double>> x = std::get<1>(results);
    std::vector<double> boundaries = getJenksBreaks(data, lowerClassLimits, nClasses);
    return std::tuple<std::vector<std::vector<double>>, std::vector<double>>{separateData(data, boundaries), boundaries};
}
