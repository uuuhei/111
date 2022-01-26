#include <iostream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <numeric>
#include <random>
#include <stdexcept>
#include <cstdio>
#include <cmath>

#include <boost/program_options.hpp>

// Zack's library definining the class 'Individual'
#include "individual.h"


using namespace std;
namespace po = boost::program_options;

// Convenience functions
double randomdouble(double a, double b);

// Evolving a population and the traits of its constituent individuals
std::vector<double> getMutList(int mutCount, double traitStDev); //记得修改trait standard deviation
void pickParents(vector<Individual> &population, double totalFitness, vector<Individual> &parents, double target, vector<double> &fitnessArr);
Individual mateParents(vector<Individual> parents, std::vector<double> mutList);
bool checkBaby(Individual baby);
void evolvePop(vector<Individual> &population, double target, int popSize, const std::vector<double> mutList);

// Extracting data
double getTotalFitness(vector<Individual> &population, double target, vector<double> &fitnessArr);
double getAverageLB(vector<Individual> &population);
double getVmor(vector<Individual> &population);
string getData(vector<Individual> &population);


int main(int argc, char** argv)
{
    
    // Standard deviation of log body size
    double traitStDev = 0.587405;
    
    po::options_description desc("simulate jaw evolution");
    desc.add_options()
        ("help", "Show help message")
        ("mutCount", po::value<int>()->default_value(1000), "Number of mutations we will have")
        ("popSize", po::value<int>()->default_value(100), "Size of the populations")
        ("endpointsensitivity", po::value<double>()->default_value(0.01), "How close average fitness gets to 1")
        ("reps", po::value<int>()->default_value(1), "Number of simulation repetitions")
        ("gen_limit", po::value<int>()->default_value(30000), "Maximum number of generations before cutting short simulation run")
        ("min_lb", po::value<double>()->default_value(1), "log of minimum body size")
        ("max_lb", po::value<double>()->default_value(3), "log of maximum body size")
        ("burnLength", po::value<double>()->default_value(0), "Length of burnin and burnout")
	 // ("output_prefix", po::value<int>()->default_value(""), "Prefix for output files")
        ("output_suffix", po::value<string>()->default_value(""), "Prefix for output files");



    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);
    const int mutCount = vm["mutCount"].as<int>();
    const int popSize = vm["popSize"].as<int>();
    const double endpointsensitivity = vm["endpointsensitivity"].as<double>(); // used 0.001 in Zack's sims -- using 0.05 you get average population fitness of <0.05
    const int reps = vm["reps"].as<int>();
    const int gen_limit = vm["gen_limit"].as<int>();
    const double min_lb = vm["min_lb"].as<double>();
    const double max_lb = vm["max_lb"].as<double>();
    const double burnLength = vm["burnLength"].as<double>(); // originally set to popSize/2
    // const string output_prefix = vm["output_prefix"].as<string>();
    const string output_suffix = vm["output_suffix"].as<string>();
    
    // targets
    const double lowerBound = 2.176091; // log_10(150)
    const double upperBound = 2.176091;

    // outputing the start parameter data
    ofstream pFile("parameters_" + output_suffix + ".txt");
    pFile << "popSize\thigh_bd\tlow_bd\tmutCount\tendpointsensitivity" << endl; // 
    pFile << popSize << "\t" << pow(10,upperBound) << "\t" << pow(10,lowerBound) << "\t" << mutCount << "\t" << endpointsensitivity << endl;
    pFile.close();

    // output stream
    ofstream outputFile("endpoints_" + output_suffix + ".txt");
    outputFile << "type\tvar_bd\tavg_bd\tavg_fit\tgens_passed" << endl;

    // looping whole sim
    for(int i = 0; i < reps; i++) {
        srand (time(0));
        cout << " " << endl;
        cout << "simulation rep: " << output_suffix << " ------------------------------------------------------------------------" << endl;
        cout << " " << endl;
        std::vector<double> mutEffectsHigh(mutCount, 0);
        std::vector<double> mutEffectsLow(mutCount, 0);
        
        // generation files opening
        
        std::ostringstream h;
        std::ostringstream l;
        h << "high_" << output_suffix << ".txt";
        l << "low_" << output_suffix << ".txt";

        string hName = h.str();
        string lName = l.str();

        ofstream highFile(hName.c_str());
        ofstream lowFile(lName.c_str());
        highFile << "gen\tavg_bd\tvar_bd\n";
        lowFile << "gen\tavg_bd\tvar_bd\n";

        // end of generation output file data


        // get a mutation list
        std::vector<double> mutList;
        mutList = getMutList(mutCount, traitStDev);

        // get phenotypes for starting population
        double start_lb = 2;
        
        // create populations of popSize
        Individual squirrel(start_lb);
        std::vector<Individual> highPop(popSize, squirrel);
        std::vector<Individual> lowPop(popSize, squirrel);

        // burn-in the populations

        
        if (burnLength != 0) {
            for(int i = 0; i < burnLength; i++) {
                evolvePop(highPop, start_lb, popSize, mutList);
                evolvePop(lowPop, start_lb, popSize, mutList);

                if(i == 0) {
                    cout << "first burnin generation: " << endl;
                    cout << "   low body_size burnin: " << pow(10, (getAverageLB(lowPop)));
                    std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness low pop: " << getTotalFitness(lowPop, lowerBound, fitnessArr) / popSize << endl;
                    cout << "   high body_size burnin: " << pow(10, (getAverageLB(highPop)));
                    // std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness high pop: " << getTotalFitness(highPop, upperBound, fitnessArr) / popSize << endl;
                    cout << " " << endl;
                }

                if(i == burnLength - 1) {
                    cout << "final burnin generation " << "(" << i + 1 << ")" << ": " << endl;
                    cout << "   low body size burnin: " << pow(10, (getAverageLB(lowPop)));
                    std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness low pop: " << getTotalFitness(lowPop, lowerBound, fitnessArr) / popSize << endl;
                    cout << "   high body size burnin: " << pow(10, (getAverageLB(highPop)));
                    // std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness high pop: " << getTotalFitness(highPop, upperBound, fitnessArr) / popSize << endl;
                    cout << " " << endl;
                    cout << "starting simulation" << endl;
                }
            }
        }

        // evolve the populations

        int count = 0;
        int lowCount = 0;
        int highCount = 0;
        bool lowFinished = false;
        bool highFinished = false;
        bool done = false;
        int burnoutLength = burnLength;
        int outputFrequency = 50; // how often (in terms of generations) to print out data
        time_t start = time(0);
        
        // while the average body size is not at the target
        while(!done && count < gen_limit) {
            
            // outputing all the generation level data:
            if(count % outputFrequency == 0) {
                if(!lowFinished) {
                    string lowData = getData(lowPop);
                    lowFile << count << "\t" << lowData;
                }
                if(!highFinished) {
                    string highData = getData(highPop);
                    highFile << count << "\t" << highData;
                }
            }

            // checking if low pop has finished evolving
            if(fabs(pow(10, getAverageLB(lowPop)) - pow(10, lowerBound)) > endpointsensitivity) {
                evolvePop(lowPop, lowerBound, popSize, mutList);
                if(count % outputFrequency == 0) {
                    cout << "generation: " << count;
                    cout << "   body size low pop: " << pow(10, (getAverageLB(lowPop))); //编一个新的function，得到正常bodysize的数据 - 搞定
                    std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness low pop: " << getTotalFitness(lowPop, lowerBound, fitnessArr) / popSize << endl;
                    
                    // warnings
                    double lowTraitVal = getAverageLB(lowPop);
                    if(lowTraitVal < (log(10.1))) {
                        cout << "   Low pop. approaching lower bound on body size   ";
                    } else if(lowTraitVal > (log(990.0))) {
                        cout << "   Low pop. approaching upper bound on body size   ";
                    }
                }
            } else if(!lowFinished) {
                lowFinished = true;
                lowCount = count;
                if(lowFinished && highFinished) {
                    done = true;
                }
            }

            // checking if high pop has finished evolving
            if(fabs(pow(10, upperBound) - pow(10, getAverageLB(highPop))) > endpointsensitivity) {
                evolvePop(highPop, upperBound, popSize, mutList);
                if(count % outputFrequency == 0) {
                    cout << "generation: " << count;
                    cout << "   Body Size high pop: " << pow(10, (getAverageLB(highPop)));
                    std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness high pop: " << getTotalFitness(highPop, upperBound, fitnessArr) / popSize << endl;
                    
                    // warnings
                    double highTraitVal = getAverageLB(highPop);
                    if(highTraitVal < log(10.1)) {
                        cout << "   High pop. approaching lower bound on body size   ";
                    } else if(highTraitVal > log(990.0)) {
                        cout << "   High pop. approaching upper bound on body size   ";
                    }
                }
            } else if(!highFinished) {
                highCount = count;
                highFinished = true;
                if(lowFinished && highFinished) {
                    done = true;
                }
            }
            
            count++;
        }

        // print final stats 记得再写一个算式把原bodysize补上
        cout << "ending simulation" << endl;
        cout << " " << endl;
        cout << "final stats: " << endl;
        cout << "high pop body size: " << pow(10, (getAverageLB(highPop))) << " " << "   high pop generations: " << highCount << endl;
        cout << "low pop body size: " << pow(10, (getAverageLB(lowPop))) << " " << "   low pop generations: " << lowCount << endl;
        
        // print total time
        double diff = difftime(time(0), start);
        cout << "total time: " << diff << " seconds" << endl;

        // print endpoint data
        std::vector<double> fitnessArr(popSize, 0.0);
        outputFile << "low\t" << getVmor(lowPop) << "\t" << pow(10, (getAverageLB(lowPop)))<< "\t" << getTotalFitness(lowPop, lowerBound, fitnessArr) / popSize << "\t" << lowCount << endl;
        outputFile << "high\t" << getVmor(highPop) << "\t" << pow(10, (getAverageLB(highPop))) << "\t" << getTotalFitness(highPop, upperBound, fitnessArr) / popSize << "\t" << highCount << endl;

        // closing output generation files
        highFile.close();
        lowFile.close();
        
    }
    return 0;
}


// Get a random number uniformly distributed between a, b

double randomdouble(double a, double b) 
{
    double random = ((double) rand()) / (double) RAND_MAX;
    double diff = b - a;
    double r = random * diff;
    return a + r;
}


/* Get a mutation list. Note that all the mutations are drawn in advance from zero-mean normal
 * distributions; when applying a mutation to a newly generated baby, we just assign it a value
 * from this pre-assembled list instead of actually mutating its parents' values. This function
 * takes a standard deviation as its second argument and assumes a mean of 0.
 */

std::vector<double> getMutList(int mutCount, const double traitStDev) {
    
    // source of randomness for initializing random seed (http://stackoverflow.com/a/38245134)
    std::random_device randomnessSource;
    // define univariate normal distributions
    std::normal_distribution<double> lbMutMaker(0.0, traitStDev);
    // store draws from the normal distributions specified above using a for loop
    std::vector<double> mutList(mutCount, 0.0);
    for (int i = 0; i < mutCount; ++i) {
        // Seed obtained using the Mersenne twister pseudorandom number generator above
        std::mt19937 gen(randomnessSource());
        mutList[i] = lbMutMaker(gen);
    }
    
    return mutList;
}
    

/* Overloaded version of the previous function, which takes a variance-covariance matrix as
 * its second argument, meaning that mutations are drawn from a multivariate normal.
 */

//直接把多元正态分布删除了 因为只有一个变量


// Randomly pick two parents from the current generation

void pickParents(vector<Individual> &population, double totalFitness, vector<Individual> &parents, double target, vector<double> &fitnessArr) {
    int parentsPicked = 0;

    while(parentsPicked < 2) {
        double tempFitness = 0;
        double rand = randomdouble(0, totalFitness);
        for(int i = 0; i < population.size(); i++) {
            tempFitness = tempFitness + fitnessArr[i];
            // cout << tempFitness << " " << rand << endl;
            if(tempFitness > rand) {
                parents[parentsPicked] = population[i];
                parentsPicked++;
                break;
            }
        }
    }
}


// Make a baby from two parents

Individual mateParents(vector<Individual> parents, std::vector<double> mutList) {
    Individual baby(parents[0].phenotypeVal);
    int mutRate = 100000;
    int val;
    for(int i = 0; i < baby.lbMutations.size(); i++) {
        val = rand();
        if(val % mutRate == 0) {
            baby.lbMutations[i] = mutList[i];
        } else if(val %2 == 0){
            baby.lbMutations[i] = parents[0].lbMutations[i];
        } else{
            baby.lbMutations[i] = parents[1].lbMutations[i];
        }
    }

    return baby;
}


// Make sure that a baby satisfies the empirical constraints on trait values

bool checkBaby(Individual baby) {

    double lb = baby.getLB();

    
    if((lb < 1 || lb > 3)){
        int stopper = 0;
        return false;
    }
    return true;
}


// Step the population forward by 1 generation

void evolvePop(vector<Individual> &population, double target, int popSize, const std::vector<double> mutList) {

    std::vector<double> fitnessArr(popSize, 0.0);
    // get total fitness
    double totalFitness = getTotalFitness(population, target, fitnessArr);
    
    // loop for length of population, make babies
    int counter = 0;
    int stallCount = 0;
    Individual ind = population[0];
    std::vector<Individual> parents(2, ind);
    std::vector<Individual> newPop(popSize, ind);
    while(counter < popSize) {
        // pick parents
        pickParents(population, totalFitness, parents, target, fitnessArr);
        // make a baby
        Individual baby = mateParents(parents, mutList);

        if(!checkBaby(baby)) {
            stallCount++;
            if(stallCount > 500) {
                cout << "problem" << endl;
                stallCount = 0;
            }
        } else {
            newPop[counter] = baby;
            // add baby to newPop
            counter++;
        }
    }

    population = newPop;
    newPop.clear();
}


// Get the fitness of an entire population

double getTotalFitness(vector<Individual> &population, double target, vector<double> &fitnessArr) {
    double totalFitness = 0;
    for(int i = 0; i < population.size(); i++) {
        double currentFitness = population[i].getFitness(target);
        fitnessArr[i] = currentFitness;
        totalFitness += currentFitness;
    }
    return totalFitness;
}


// Return the average log(body size)

double getAverageLB(vector<Individual> &population) {
    double totalLB = 0;
    for(int i = 0; i < population.size(); i++){
        totalLB += population[i].getLB();
    }
    return totalLB/population.size();
}


// Get the total morphological variation (sum of trait variances) in the population

double getVmor(vector<Individual> &population) {

    int popSize = population.size();
    // creating a vector of vectors to hold all the phenotype vals
    std::vector<double> popVects(popSize, 0.0);
  
    // puting the phenotype vals in the vectors
    for(int i = 0; i < popSize; i++) {
        popVects[i] = population[i].getLB();

    }

    // getting body size mean
    double lbMean = (std::accumulate(popVects.begin(), popVects.end(), 0.0)) / popSize;


    // getting body size variance
    double lbVar = 0.0;
    for(int j = 0; j < popSize; j++) {
        lbVar += ((popVects[j] - lbMean) * (popVects[j] - lbMean)) / popSize;
    }

    return lbVar;
}


// Get the population means and variances of individual traits - 需要修改，已经不再是individual traits了

string getData(vector<Individual> &population) {

    std::stringstream s;
    int popSize = population.size();
    // creating a vector of vectors to hold all the phenotype vals
    std::vector<double> popVects(popSize, 0.0);
  
    // puting the phenotype vals in the vectors
    for(int i = 0; i < popSize; i++) {
        popVects[i] = population[i].getLB();
    }

    // getting body size mean
    double lbMean = (std::accumulate(popVects.begin(), popVects.end(), 0.0)) / popSize;

    s << lbMean << "\t";

    // getting body size variance 也没什么用，只有一个trait
    double lbVar = 0.0;
    for(int j = 0; j < popSize; j++) {
        lbVar += ((popVects[j] - lbMean) * (popVects[j] - lbMean)) / popSize;
    }

    s << lbVar << "\n";

    string data = s.str();
    return data;
}

