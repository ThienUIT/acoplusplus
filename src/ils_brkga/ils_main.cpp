#include "data.h"
#include "util.h"
#include "MTRand.h"

MTRand rng;				// initialize the random number generator          
 
std::pair<bool, int> checkSolution(string solutionFile) {
    
    string line_tour, line_items;
    stringstream ss;

    ifstream fin(solutionFile.c_str());

    if(!fin) {
        clog<<"ERROR!"<<endl;
        exit(0);
    }

    getline(fin, line_tour);
    getline(fin, line_items);
    
    fin.close();
    
    for(int j=0;j<(int)line_tour.length();j++) {
        if(line_tour[j] < '0' || line_tour[j] > '9') line_tour[j] = ' ';
    }
    ss.clear();
    ss<<line_tour;
    
    set<int> st_tour;
    vector<int> tour;
    int v;
    
    tour.push_back(1);
    st_tour.insert(1);
    while(ss>>v) {
        if(v <= 1 || v >= data.numCities) {
            return make_pair(false, -1);
        }
        tour.push_back(v);
        st_tour.insert(v);
    }    
    tour.push_back(data.numCities);
    st_tour.insert(data.numCities);
    
     if((int)st_tour.size() != (int)tour.size()) {
        clog<<"Erro! as cidades repetem "<< (int)st_tour.size() << " != "<<(int)tour.size()<<endl;
        return make_pair(false, -1);
    }
    
    for(int j=0;j<(int)line_items.length();j++) {
        if(line_items[j] < '0' || line_items[j] > '9') line_items[j] = ' ';
    }
    ss.clear();
    ss<<line_items;
    
    set<int> st_items;
    vector<int> items;
    while(ss>>v) {
        items.push_back(v);
        st_items.insert(v);
    }    
    
    if((int)st_items.size() != (int)items.size()) {
        clog<<"Erro! os items repetem "<< (int)st_items.size() << " != "<<(int)items.size()<<endl;
        return make_pair(false, -1);   
    }
    
    vector<pair<int, int> > weight_profit;
    weight_profit.resize((int)tour.size()+1);
    
    for(int i=0;i<(int)tour.size();i++) {
        weight_profit[i].first = weight_profit[i].second = 0;
        for(int j=0;j<(int)items.size();j++) {
            if(data.items[items[j]].idCity == tour[i]) {
                weight_profit[i].first += data.items[items[j]].weight;
                weight_profit[i].second += data.items[items[j]].profit;
            }
        }
    }
    
    for(int i=0;i<(int)tour.size();i++) {
        weight_profit[i].first = weight_profit[i].second = 0;
        for(int j=0;j<(int)items.size();j++) {
            if(data.items[items[j]].idCity == tour[i]) {
                weight_profit[i].first += data.items[items[j]].weight;
                weight_profit[i].second += data.items[items[j]].profit;
            }
        }
    }
    for(int j=0;j<(int)items.size();j++) {
        bool ok = false;
        for(int i=0;i<(int)tour.size();i++) {
            if(data.items[items[j]].idCity == tour[i]) {
                ok = true;
            }
        }
        if(!ok) {
            clog<<"Erro! nenhuma cidade passa onde o item "<<items[j]<<" estah"<<endl;
            return make_pair(false, -1);    
        }
    }

    double currentCapacityOfKnapsack = 0;
    double currentTime = 0.0;
    double currentProfit = 0.0;
    
    int prev = 1;
    
    for(int i=1;i<(int)tour.size();i++) {
        currentTime += data.distance[prev][tour[i]]/(data.maxSpeed - data.v*currentCapacityOfKnapsack);
        currentCapacityOfKnapsack += weight_profit[i].first;
        currentProfit += weight_profit[i].second;
        prev = tour[i];
    }
    
    if(currentCapacityOfKnapsack - EPS > data.capacityOfKnapsack) {
        clog<<"Error! capacidade da mochila extrapolada  "<<currentCapacityOfKnapsack<<' '<<data.capacityOfKnapsack<<endl;
        return make_pair(false, -1);
    }
           
    if(currentTime - EPS > data.maxTime) {
        clog<<"Error!  tempo extrapolado "<<currentTime<<' '<<data.maxTime<<endl;
        return make_pair(false, -1);
    }
    
    return make_pair(true, currentProfit);
}

struct Solution {
    
    vector<int> permutation;
    double fitness;
    
    int lastChange;
    
    void randomSolution() {
        for(int i=1;i<data.items.size();i++) {
            permutation.push_back(i);
        }
        
        int rnd1, rnd2;
        
        for(int i=1;i<data.items.size()/2;i++) {        
            do {
                rnd1 = rng.randInt(permutation.size()-1);
                rnd2 = rng.randInt(permutation.size()-1);
            } while (rnd1 == rnd2);
            
            swap(permutation[rnd1], permutation[rnd2]);
        }
        
        evaluate();
    }
    
    Solution shaking(int k) {
    
        std::uniform_int_distribution<int> distribution(0, permutation.size()-1);
       
        int rnd1, rnd2;        
        
        Solution nS = *this;
        
        while(k--) {
        
            do {
                rnd1 = rng.randInt(permutation.size()-1);
                rnd2 = rng.randInt(permutation.size()-1);
            } while (rnd1 == rnd2);
            
            swap(nS.permutation[rnd1], nS.permutation[rnd2]);
        }
        
        nS.evaluate();
        
        return nS;               
    }
        
    double localSearch(double runtime) {
        
        vector<int> bestPermutation = permutation;
        double bestValue = evaluate();

        tEnd = high_resolution_clock::now();
        duration<double> exec_time = duration_cast<duration<double> >(tEnd - tBegin);

        if(exec_time.count() >= runtime) return bestValue;
        
        while(1) {
            
            bool improvement = false;
        
            double bestNeighborValue = bestValue;
            vector<int> bestNeighborPermutation;
                            
            tEnd = high_resolution_clock::now();
            duration<double> exec_time = duration_cast<duration<double> >(tEnd - tBegin);
    
            if(exec_time.count() >= runtime) break;
            
            for(int i=0;i<(int)permutation.size() && !improvement;i++) {
                    
                tEnd = high_resolution_clock::now();
                exec_time = duration_cast<duration<double> >(tEnd - tBegin);

                if(exec_time.count() >= runtime) break;
                
                for(int j=i+1;j<(int)permutation.size() && !improvement;j++) {
                
                    tEnd = high_resolution_clock::now();
                    exec_time = duration_cast<duration<double> >(tEnd - tBegin);

                    if(exec_time.count() >= runtime) break;
                    
                    swap(permutation[i], permutation[j]);
                        
                    double aaa = evaluate();
                    
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
        
        fitness = bestValue;
        return bestValue;
    }
    
    double evaluate(string solutionFileOut = "") {
        
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
                    
                    lastChange = i;
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
                    
                    lastChange = i;
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

        fitness = sumProfit;
        return sumProfit;
    }    
};

Solution iteratedLocalSearch(double runtime, string logFileName) {

    using namespace std::chrono;
    tBegin = high_resolution_clock::now();
    
    Solution s;
    
    s.randomSolution();
    s.localSearch(runtime);
    
    Solution bestSolution = s;
    
    ofstream fout(logFileName.c_str());
    fout << "Profit Time(s)" << endl;
    tEnd = high_resolution_clock::now();
    exec_time = duration_cast<duration<double> >(tEnd - tBegin);
    fout << fixed << setprecision(0) << bestSolution.fitness << ' ' << setprecision(5) << exec_time.count() << endl;
    
    Solution nS;
    
    int k = 2;
    
    while(1) {
        
        tEnd = high_resolution_clock::now();
        exec_time = duration_cast<duration<double> >(tEnd - tBegin);
        
        if(exec_time.count() >= runtime) break;
     
        nS = s.shaking(k);
        nS.localSearch(runtime);        
        
        if(nS.fitness > s.fitness) {
            s = nS;
            bestSolution = s;            
            tEnd = high_resolution_clock::now();
            exec_time = duration_cast<duration<double> >(tEnd - tBegin);
            fout << fixed << setprecision(0) << bestSolution.fitness << ' ' << setprecision(5) << exec_time.count() << endl;
            k = 2;
        }
        else {
            k++;
        }           
           
        //clog << "here" << endl;
    }
    
    fout.close();
    
    return bestSolution;
}

void usage() {
    clog << "Usage ./ilsthop --inputfile <input_file_name> --seed <seed_rnd_numbers> --time <runtime_secs> --outputfile <output_file_name>" << endl;
    exit(0);
}

int main(int argc, char* argv[]) {
    
    if(argc < 9) usage();
    
    string instanceFile;
    string solutionFile;
    int rngSeed, runtime;
    
    char parameterStr[1000];
    
    int check_parameters = 0;    
    for(int i = 1; i < argc; i += 2) { 
        if(strcmp(argv[i], "--inputfile") == 0) { sscanf(argv[i+1],"%s", parameterStr); instanceFile = string(parameterStr); check_parameters += 1; }
        else if(strcmp(argv[i], "--seed") == 0) { sscanf(argv[i+1],"%d", &rngSeed); check_parameters += 1; }
        else if(strcmp(argv[i], "--time") == 0) { sscanf(argv[i+1],"%d", &runtime); check_parameters += 1; }
        else if(strcmp(argv[i], "--outputfile") == 0) { sscanf(argv[i+1],"%s", parameterStr); solutionFile = string(parameterStr); check_parameters += 1; } 
        else check_parameters = -INF;        
    }
    
    if(check_parameters != 4) usage();

    data.readData(instanceFile);
    
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    rng = MTRand(rngSeed);				// initialize the random number generator       
    Solution solution = iteratedLocalSearch(runtime, solutionFile + ".log");

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
        
    solution.evaluate(solutionFile);
    
    pair<bool, int> result = checkSolution(solutionFile);

    if(result.first && result.second == solution.fitness) {}       
    else { clog << "Error! " << result.second << ' ' << solution.fitness << endl; exit(0); }
    
    return 0;
}
