#ifndef DATA_H
#define DATA_H

#include "util.h"


struct Item {
    int id;
    int profit;
    int weight;
    int idCity;
        
    Item(int _id, int _profit, int _weight, int _idCity) : id(_id), profit(_profit), weight(_weight), idCity(_idCity) {}
};

class Data {

    public:
            string problemName;
            string knapsackDataType;            
            int numCities;
            int numItems;
            int capacityOfKnapsack;
            double v, minSpeed, maxSpeed;
            double maxTime;
            double **distance;
            
            int totalProfit;
            
            vector<Item> items;

            void readData(string file) {

                string line;
                stringstream ss;

                ifstream fin(file.c_str());

                if(!fin) {
                    clog<<"ERROR!"<<endl;
                    exit(0);
                }
              
                getline(fin, line); // PROBLEM NAME:    a280-ThOP

                getline(fin, line); // KNAPSACK DATA TYPE: bounded strongly corr
    
                getline(fin, line); // DIMENSION:  280
                for(int j=0;j<(int)line.length();j++) {
                    if(line[j] < '0' || line[j] > '9') line[j] = ' ';
                }
                ss.clear();
                ss<<line;
                ss>>numCities;

                getline(fin, line); // NUMBER OF ITEMS:    279
                for(int j=0;j<(int)line.length();j++) {
                    if(line[j] < '0' || line[j] > '9') line[j] = ' ';
                }
                ss.clear();
                ss<<line;
                ss>>numItems;
   
                getline(fin, line); // CAPACITY OF KNAPSACK:   25936
                for(int j=0;j<(int)line.length();j++) {
                    if(line[j] < '0' || line[j] > '9') line[j] = ' ';
                }
                ss.clear();
                ss<<line;
                ss>>capacityOfKnapsack;
                
                getline(fin, line); // MAX TIME:   xxx
                for(int j=0;j<(int)line.length();j++) {
                    if(line[j] < '0' || line[j] > '9') line[j] = ' ';
                }
                ss.clear();
                ss<<line;
                ss>>maxTime;
                
                getline(fin, line); // MIN SPEED:  0.1
                for(int j=0;j<(int)line.length();j++) {
                    if((line[j] < '0' || line[j] > '9') && line[j] != '.') line[j] = ' ';
                }
                ss.clear();
                ss<<line;
                ss>>minSpeed;
                
                getline(fin, line); // MAX SPEED:  1
                for(int j=0;j<(int)line.length();j++) {
                    if((line[j] < '0' || line[j] > '9') && line[j] != '.') line[j] = ' ';
                }
                ss.clear();
                ss<<line;
                ss>>maxSpeed;
                
                getline(fin, line); // EDGE_WEIGHT_TYPE:   CEIL_2D
                
                getline(fin, line); // NODE_COORD_SECTION  (INDEX, X, Y): 
                
                v = (maxSpeed - minSpeed)/capacityOfKnapsack;
                
		vector<pair<double, double> > vertex;
                vertex.resize(numCities+1);
                
                distance = new double*[numCities+1];
                for(int i=0;i<numCities+1;i++) {
                    distance[i] = new double[numCities+1];
                }

		int id;
                for(int i=1;i<=numCities;i++) {
                    fin>>id;
                    fin>>vertex[i].first;
                    fin>>vertex[i].second;
                }
                
                for(int i=1;i<=numCities;i++) {
                    for(int j=i;j<=numCities;j++) {
                        distance[i][j] = distance[j][i] = (int) ceil(   
                                                                        sqrt(
                                                                            ((vertex[i].first-vertex[j].first)*(vertex[i].first-vertex[j].first)) +
                                                                            ((vertex[i].second-vertex[j].second)*(vertex[i].second-vertex[j].second))
                                                                        )
                                                                    );
                    }
                }

                getline(fin, line); // '\n'
                getline(fin, line); // ITEMS SECTION    (INDEX, PROFIT, WEIGHT, ASSIGNED NODE NUMBER):

                int id_city;
                double p, w;

                items.push_back(Item(0, -INF, INF, -1));
                totalProfit = 0;
                for(int i=1;i<=numItems;i++) {
                    fin>>id>>p>>w>>id_city;
                    items.push_back(Item(id, p, w, id_city));
                    totalProfit += p;
                }

                fin.close();
            }
};

Data data;

#endif
