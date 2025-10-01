#include "data.h"
#include "util.h"
#include <iostream>
#include "Decoder.h"
#include "MTRand.h"
#include "BRKGA.h"

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

void usage() {
    clog << "Usage ./brkgathop --inputfile <input_file_name> --p <pop_size> --pe <elite_pop_size> --pm <mutant_pop_size> --rhoe <inheritance_prob> --seed <seed_rnd_numbers> --time <runtime_secs> --outputfile <output_file_name>" << endl;
    exit(0);
}

int main(int argc, char* argv[]) {
    
    if(argc < 19) usage();
    
    string instanceFile;
    string solutionFile;
    int rngSeed, runtime, lsfreq;       
    unsigned p;           // size of population
    double pe;            // fraction of population to be the elite-set
    double pm;            // fraction of population to be replaced by mutants
    double rhoe;          // probability that offspring inherit an allele from elite parent
    
    const unsigned K = 1;              // number of independent populations
    const unsigned MAXT = 1;           // number of threads for parallel decoding
    
    bool calibrationModeOn = false;
    
    char parameterStr[1000];
    
    int check_parameters = 0;    
    for(int i = 1; i < argc; i += 2) { 
        if(strcmp(argv[i], "--inputfile") == 0) { sscanf(argv[i+1],"%s", parameterStr); instanceFile = string(parameterStr); check_parameters += 1; }
        else if(strcmp(argv[i], "--p") == 0) { sscanf(argv[i+1],"%u", &p); check_parameters += 1; }
        else if(strcmp(argv[i], "--pe") == 0) { sscanf(argv[i+1],"%lf", &pe); check_parameters += 1; }
        else if(strcmp(argv[i], "--pm") == 0) { sscanf(argv[i+1],"%lf", &pm); check_parameters += 1; }
        else if(strcmp(argv[i], "--rhoe") == 0) { sscanf(argv[i+1],"%lf", &rhoe); check_parameters += 1; }
        else if(strcmp(argv[i], "--lsfreq") == 0) { sscanf(argv[i+1],"%d", &lsfreq); check_parameters += 1; }
        else if(strcmp(argv[i], "--seed") == 0) { sscanf(argv[i+1],"%d", &rngSeed); check_parameters += 1; }
        else if(strcmp(argv[i], "--time") == 0) { sscanf(argv[i+1],"%d", &runtime); check_parameters += 1; }
        else if(strcmp(argv[i], "--outputfile") == 0) { sscanf(argv[i+1],"%s", parameterStr); solutionFile = string(parameterStr); check_parameters += 1; } 
        else if(strcmp(argv[i], "--calibration") == 0) { calibrationModeOn = true; i -= 1; } 
        else check_parameters = -INF;        
    }
    
    if(check_parameters != 9) usage();
    
    data.readData(instanceFile);
        
    const unsigned n = data.numItems;  // size of chromosomes    
    
    if(calibrationModeOn) runtime = ceil(data.numItems/10.0);
    
    Decoder decoder;            // initialize the decoder
    decoder.maxTime = runtime;
    decoder.lsfreq = lsfreq;
    
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
        
    MTRand rng(rngSeed);                // initialize the random number generator
    
    // initialize the BRKGA-based heuristic
    BRKGA< Decoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
    
    unsigned generation = 0;        // current generation
    unsigned generationNoImprovement = 0; // number of predecessor generations without improvement
    const unsigned X_INTVL = 50;    // exchange best individuals at every X_INTVL generations
    const unsigned X_NUMBER = 1;    // exchange top X_NUMBER best
    const unsigned NUM_GENS_NO_IMPROVEMENT = 100;    // run while while not reaching NUM_GENS_NO_IMPROVEMENT generations without improvement

    using namespace std::chrono;
    tBegin = high_resolution_clock::now();
    
    double prevBestFitness = -1;
    
    ofstream fout;
    if(!calibrationModeOn) {
        std::string logFileName = solutionFile + ".log";
        fout.open(logFileName.c_str());
        fout << "Profit Time(s)" << endl;
        tEnd = high_resolution_clock::now();
        exec_time = duration_cast<duration<double> >(tEnd - tBegin);
        fout << fixed << setprecision(0) << -algorithm.getBestFitness() << ' ' << setprecision(5) << exec_time.count() << endl;
        prevBestFitness = -algorithm.getBestFitness();
    }
        
    while(1) {
        
        tEnd = high_resolution_clock::now();
        exec_time = duration_cast<duration<double> >(tEnd - tBegin);

        if(exec_time.count() >= runtime) break;
        
        numIter = generation + 1; 
        algorithm.evolve();    // evolve the population for one generation
        
        if(!calibrationModeOn && prevBestFitness != -algorithm.getBestFitness()) {
            prevBestFitness = -algorithm.getBestFitness();
            tEnd = high_resolution_clock::now();
            exec_time = duration_cast<duration<double> >(tEnd - tBegin);
            fout << fixed << setprecision(0) << -algorithm.getBestFitness() << ' ' << setprecision(5) << exec_time.count() << endl;
        }
        
        ++generation;
    }                 
    
    if(!calibrationModeOn) fout.close();
    
    std::pair<double, std::vector< double > > solution = std::make_pair(-algorithm.getBestFitness(), algorithm.getBestChromosome());
    
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
    
    if(calibrationModeOn) cout << -solution.first << endl;
    else {
        decoder.decode(solution.second, solutionFile);
        pair<bool, int> result = checkSolution(solutionFile);
        if(result.first && result.second == solution.first) {}       
        else { clog << "Error! " << result.second << ' ' << solution.first << endl; exit(0); }
    }

    return 0;
}
