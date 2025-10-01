/*
 * Decoder.h
 *
 * Any decoder must have the format below, i.e., implement the method decode(vector< double >&)
 * returning a double corresponding to the fitness of that vector. If parallel decoding is to be
 * used in the BRKGA framework, then the decode() method _must_ be thread-safe; the best way to
 * guarantee this is by adding 'const' to the end of decode() so that the property will be checked
 * at compile time.
 *
 * The chromosome inside the BRKGA framework can be changed if desired. To do so, just use the
 * first signature of decode() which allows for modification. Please use double values in the
 * interval [0,1) when updating, thus obeying the BRKGA guidelines.
 *
 *  Created on: Jan 14, 2011
 *      Author: rtoso
 */

#ifndef DECODER_H
#define DECODER_H

#include "data.h"
#include "util.h"

double evaluate(vector<int> &permutation, string solutionFileOut = "") {
    
    vector<unsigned> tour;
    vector<unsigned> items;
    
    int prev = 1;
    double currentW = 0.0;
    double currentTime = 0.0;
    double sumProfit = 0.0;
    vector<bool> inserted(data.numCities+1, false);
    
    int minWeight = 987654321;
    for(unsigned i = 0; i < data.items.size(); ++i) {
        minWeight = min(minWeight, data.items[i].weight);
    }
    
    for(unsigned i = 0; i < permutation.size(); ++i) {
        
        if(data.capacityOfKnapsack - currentW < minWeight) break;
        
        int IDItem = permutation[i];
        int IDCity = data.items[IDItem].idCity;
        
        if(inserted[IDCity] == true) {
        
            int _prev = 1;
            double _currentW = 0.0;
            double _currentTime = 0.0;
            double _sumProfit = 0.0;
            for(unsigned i = 0; i < tour.size(); ++i) {
                if(tour[i] == IDCity) continue;
                _currentTime += data.distance[_prev][tour[i]]/(data.maxSpeed - data.v*_currentW);
                for(unsigned j = 0; j < items.size(); ++j) {
                    if(data.items[items[j]].idCity == tour[i]) {
                        _currentW += data.items[items[j]].weight;
                        _sumProfit += data.items[items[j]].profit;
                    }
                }
                _prev = tour[i];
            }
            
            _currentTime += data.distance[_prev][IDCity] / (data.maxSpeed - data.v * _currentW);                    
            
            for(unsigned j = 0; j < items.size(); ++j) {
                if(data.items[items[j]].idCity == IDCity) {
                    _currentW += data.items[items[j]].weight;
                    _sumProfit += data.items[items[j]].profit;
                }
            }           
            _currentW += data.items[IDItem].weight;
            _sumProfit += data.items[IDItem].profit;
            
            if  (
                    (
                    _currentTime + data.distance[IDCity][data.numCities] / (data.maxSpeed - data.v * _currentW)<= data.maxTime
                    )
                    &&
                    (
                    _currentW  <= data.capacityOfKnapsack                    
                    )
                )
                {
                
                _prev = IDCity;
                
                tour.erase(remove(tour.begin(), tour.end(), IDCity), tour.end());
                tour.push_back(IDCity);
                items.push_back(IDItem);
                currentTime = _currentTime;
                currentW = _currentW;
                sumProfit = _sumProfit;
                prev = _prev;
            }                    
        }
        else {
            
            if  (
                    (
                    currentTime + data.distance[prev][IDCity]/(data.maxSpeed - data.v*currentW) + data.distance[IDCity][data.numCities]/(data.maxSpeed - data.v * (currentW + data.items[IDItem].weight)) <= data.maxTime
                    )
                    &&
                    (
                    currentW + data.items[IDItem].weight <= data.capacityOfKnapsack                    
                    )
                ) {
        
                currentTime += data.distance[prev][IDCity]/(data.maxSpeed - data.v*currentW);
                currentW += data.items[IDItem].weight;    
                sumProfit += data.items[IDItem].profit;
                
                prev = IDCity;
                
                items.push_back(IDItem);
                tour.push_back(IDCity);
                inserted[IDCity] = true;
            }
        }
    }
    
    currentTime += data.distance[prev][data.numCities]/(data.maxSpeed - data.v*currentW);

    if(!solutionFileOut.empty()) {
        
        ofstream fout(solutionFileOut.c_str());
        
        fout<<"[";
        //cout<<"[";
        if(tour.size() > 0) {
            fout<<tour[0];
            //cout<<tour[0];
        }
        for(unsigned i = 1; i < tour.size(); ++i) {
            fout<<','<<tour[i];
            //cout<<','<<tour[i];
        }
        fout<<"]"<<endl;
        //cout<<"]"<<endl;
        
        sort(items.begin(), items.end());
        
        fout<<"[";
        //cout<<"[";
        if(items.size() > 0) {
            fout<<items[0];
            //cout<<items[0];
        }
        for(unsigned i = 1; i < items.size(); ++i) {
            fout<<','<<items[i];
            //cout<<','<<items[i];
        }
        fout<<"]"<<endl;
        //cout<<"]"<<endl;
    }

    return sumProfit;
}    

double localSearch(vector<int> &permutation, double maxTime = 987654321, bool ls = false) {
        
    vector<int> bestPermutation = permutation;
    double bestValue = evaluate(permutation);
    
    if(ls == false) return bestValue;
    
    tEnd = high_resolution_clock::now();
    duration<double> exec_time = duration_cast<duration<double> >(tEnd - tBegin);

    if(exec_time.count() >= maxTime) return bestValue;

    while(1) {
            
        bool improvement = false;
    
        double bestNeighborValue = bestValue;
        vector<int> bestNeighborPermutation;
                        
        tEnd = high_resolution_clock::now();
        duration<double> exec_time = duration_cast<duration<double> >(tEnd - tBegin);

        if(exec_time.count() >= maxTime) break;
        
        for(int i=0;i<(int)permutation.size() && !improvement;i++) {
                
            tEnd = high_resolution_clock::now();
            exec_time = duration_cast<duration<double> >(tEnd - tBegin);

            if(exec_time.count() >= maxTime) break;
            
            for(int j=i+1;j<(int)permutation.size() && !improvement;j++) {
            
                tEnd = high_resolution_clock::now();
                exec_time = duration_cast<duration<double> >(tEnd - tBegin);

                if(exec_time.count() >= maxTime) break;
                
                swap(permutation[i], permutation[j]);
                    
                double aaa = evaluate(permutation);
                
                if(aaa > bestNeighborValue) {
                    bestNeighborValue = aaa;
                    bestNeighborPermutation = permutation;
                    improvement = true;
                }
                swap(permutation[i], permutation[j]);
            }
        }
        
        if(!improvement) break;
        
        permutation = bestNeighborPermutation;
        bestValue = bestNeighborValue;
    }
    
    return bestValue;
}

class Decoder {

public:
    
    double maxTime;
    int lsfreq;
        
	Decoder() {}
        
	~Decoder() {}

        double decode(vector< double >& chromosome, string solutionFileOut = "") const {

            vector< pair<double, int > > subset;
            
            for(int i = 0; i < (int) chromosome.size(); ++i) {                
                if(chromosome[i] < -0.5) continue;
                int IDItem = i + 1;                
                subset.push_back(make_pair(chromosome[i], IDItem));
            }
            
            sort(subset.begin(), subset.end());
            
            vector<int> subpermutation;
            for(int i=0;i<(int)subset.size();i++) {
                subpermutation.push_back(subset[i].second);
            }
            
            double result = localSearch(subpermutation, maxTime, (numIter % lsfreq == 0) && solutionFileOut.empty());
            
            if((numIter % lsfreq == 0)) {
                for(int i=0;i<(int)subpermutation.size();i++) {
                    chromosome[subpermutation[i]-1] = 0.00001 + ((1.0/1000.0) * i);
                }
            }
            
            if(!solutionFileOut.empty()) {
                evaluate(subpermutation, solutionFileOut);
            }
            
            return -result;
        }

private:
};

#endif
