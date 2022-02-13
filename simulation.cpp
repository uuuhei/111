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
int randominteger(int a, int b);

// Evolving a population and the traits of its constituent individuals
std::vector<double> getMutList(int mutCount, double traitStDev); // trait standard deviation
Individual pickAndMateParents(std::vector<Individual> &population, double totalFitness, double target, std::vector<double> &fitnessArr, std::vector<double> mutList);
bool checkBaby(Individual baby, double lowerLim, double upperLim);
void evolvePop(vector<Individual> &population, double target, double lowerLim, double upperLim, int popSize, const std::vector<double> mutList, int numShelters, double calamFreq, double calamStrength, bool showShelterStats);

// Extracting data
double getTotalFitness(vector<Individual> &population, double target, vector<double> &fitnessArr);
double getAverageLB(vector<Individual> &population);
double getVmor(vector<Individual> &population);
string getData(vector<Individual> &population, double target);


int main(int argc, char** argv)
{
    
    // Standard deviation of log body size
    double traitStDev = 0.05;
    
    po::options_description desc("simulate jaw evolution");
    desc.add_options()
        ("help", "Show help message")
        ("mutCount", po::value<int>()->default_value(1000), "Number of mutations we will have")
        ("popSize", po::value<int>()->default_value(100), "Size of the populations")
        ("endpointsensitivity", po::value<double>()->default_value(0.01), "How close average fitness gets to 1")
        ("reps", po::value<int>()->default_value(1), "Number of simulation repetitions")
        ("gen_limit", po::value<int>()->default_value(30000), "Maximum number of generations before cutting short simulation run")
        ("min_lb", po::value<double>()->default_value(1), "Log of minimum body size")    // log_10(10)
        ("max_lb", po::value<double>()->default_value(3), "Log of maximum body size")    // log_10(1000)
        ("start_lb", po::value<double>()->default_value(2), "Log of starting body size") // log_10(100)
        ("nShelters", po::value<int>()->default_value(10), "Number of shelters available to each generation") // Set to 10% of popSize by default
        ("calamityFrequency", po::value<double>()->default_value(0.1), "Frequency of calamities affecting unsheltered individuals")
        ("calamityStrength", po::value<double>()->default_value(0.1), "Fraction of the population eliminated by a calamity")
        ("burnLength", po::value<double>()->default_value(0), "Length of burnin and burnout")
	 // ("output_prefix", po::value<int>()->default_value(""), "Prefix for output files")
        ("output_suffix", po::value<string>()->default_value(""), "Prefix for output files");



    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);
    const int mutCount = vm["mutCount"].as<int>();
    const int popSize = vm["popSize"].as<int>();
    const double endpointsensitivity = vm["endpointsensitivity"].as<double>();
    const int reps = vm["reps"].as<int>();
    const int gen_limit = vm["gen_limit"].as<int>();
    const double min_lb = vm["min_lb"].as<double>();
    const double max_lb = vm["max_lb"].as<double>();
    const double start_lb = vm["start_lb"].as<double>();
    const int nShelters = vm["nShelters"].as<int>();
    const double calamityFrequency = vm["calamityFrequency"].as<double>();
    const double calamityStrength = vm["calamityStrength"].as<double>();
    const double burnLength = vm["burnLength"].as<double>(); // this is just the *minimum* burnin length
    // const string output_prefix = vm["output_prefix"].as<string>();
    const string output_suffix = vm["output_suffix"].as<string>();
    
    // body size optimum (log scale)
    const double targetLbs = 2.176091; // log_10(150)

    // outputing the start parameter data
    ofstream pFile("parameters_" + output_suffix + ".txt");
    pFile << "popSize\tbodySize\tmutCount\tendpointsensitivity" << endl; //
    pFile << popSize << "\t" << pow(10, targetLbs) << "\t" << mutCount << "\t" << endpointsensitivity << endl;
    pFile.close();

    // output stream
    ofstream outputFile("endpoints_" + output_suffix + ".txt");
    outputFile << "var_bd\tavg_bd\tavg_fit\tgens_passed" << endl;

    // looping whole sim
    for(int i = 0; i < reps; i++) {
        srand (time(0));
        cout << " " << endl;
        cout << "simulation rep: " << output_suffix << " ------------------------------------------------------------------------" << endl;
        cout << " " << endl;
        std::vector<double> mutEffects(mutCount, 0);
        
        // generation files opening
        
        std::ostringstream params;
        params << "params_" << output_suffix << ".txt";
        string pName = params.str();
        ofstream genFile(pName.c_str());
        genFile << "gen\tavg_bd\tvar_bd\tavg_fit\n";

        // generation files closing

        // get a mutation list
        std::vector<double> mutList;
        mutList = getMutList(mutCount, traitStDev);

        /* Get phenotypes for starting populations. If a nonzero burnin period is specified, we
         * will start with a population of individuals sharing the same randomly generated value,
         * and let it evolve (using the same evolvePop() function used for the simulation proper)
         * toward a specified starting body size. The idea is to generate a population where all
         * individuals are phenotypically identical and equally distant from the fitness optimum,
         * but which already contains some genotypic variation for selection to act upon. If burnin
         * is skipped, the simulation starts with a bunch of individuals that are identical both
         * phenotypically and genotypically.
         */
         
        double init;
        
        if (burnLength != 0) {
            init = randomdouble(min_lb, max_lb);
        } else {
            init = start_lb;
        }
        
        // Create a population of specified size, with all individuals initialized as shelter-less
        Individual squirrel(init);
        std::vector<Individual> Pop(popSize, squirrel);
        
        int outputFrequency = 50; // how often (in terms of generations) to print out data

        /* Burn-in the populations. The burnin phase will end (1) when the average log body size
         * gets acceptably close to the specified starting value OR (2) when the target number of
         * generations is reached -- whichever comes last.
         */

        if (burnLength != 0) {
            int burnCount = 0;
            while(fabs(pow(10, getAverageLB(Pop)) - pow(10, start_lb)) > endpointsensitivity || burnCount < burnLength) {
                
                evolvePop(Pop, start_lb, min_lb, max_lb, popSize, mutList, 0, 0, 0, false);
            
                if(burnCount % outputFrequency == 0) {
                    cout << "burnin generation: " << burnCount << endl;
                    cout << "Avg body size: " << pow(10, (getAverageLB(Pop)));
                    std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness: " << getTotalFitness(Pop, targetLbs, fitnessArr) / popSize << endl;
                    cout << " " << endl;
                }
                
                burnCount++;
            }
        }
        
        // Randomly select individuals to receive shelter
        std::vector<int> shelterIdx;
        for(int i = 0; i < nShelters; i++) {
            int idx = randominteger(0, popSize - 1);
            shelterIdx.push_back(idx);
        }
        
        // Give shelter to the individuals selected
        for(int j = 0; j < popSize; j++) {
            if (std::find(shelterIdx.begin(), shelterIdx.end(), j) != shelterIdx.end()) {
                Pop[j].setShelter(true);
            }
        }

        // Evolve the population
        int count = 0;
        bool finished = false;
        time_t start = time(0);
        
        // While the average body size is not at the target
        while(!finished && count < gen_limit) {
            
            // Outputing all the generation level data:
            if(!finished && count % outputFrequency == 0) {
                string data = getData(Pop, targetLbs);
                genFile << count << "\t" << data;
            }

            // Checking if the population has finished evolving
            if(fabs(pow(10, getAverageLB(Pop)) - pow(10, targetLbs)) > endpointsensitivity) {
                
                if (count % outputFrequency == 0) {
                    evolvePop(Pop, targetLbs, min_lb, max_lb, popSize, mutList, nShelters, calamityFrequency, calamityStrength, true);
                    cout << "generation: " << count;
                    cout << "   Avg body size: " << pow(10, (getAverageLB(Pop))); //编一个新的function，得到正常bodysize的数据 - 搞定
                    std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness: " << getTotalFitness(Pop, targetLbs, fitnessArr) / popSize << endl;
                    
                    // Warnings
                    if(pow(10, (getAverageLB(Pop))) < 11.0) {
                        cout << "   Approaching lower bound on body size   ";
                    } else if(pow(10, (getAverageLB(Pop))) > 990.0) {
                        cout << "   Approaching upper bound on body size   ";
                    }
                } else {
                    // Evolve quietly
                    evolvePop(Pop, targetLbs, min_lb, max_lb, popSize, mutList, nShelters, calamityFrequency, calamityStrength, false);
                }
            } else {
                finished = true;
            }
            
            count++;
        }

        // print final stats 记得再写一个算式把原bodysize补上
        cout << "ending simulation" << endl;
        cout << " " << endl;
        cout << "final stats: " << endl;
        cout << "avg body size: " << pow(10, (getAverageLB(Pop))) << " " << "   generations: " << count << endl;
        
        // print total time
        double diff = difftime(time(0), start);
        cout << "total time: " << diff << " seconds" << endl;

        // print endpoint data
        std::vector<double> fitnessArr(popSize, 0.0);
        outputFile << getVmor(Pop) << "\t" << pow(10, (getAverageLB(Pop)))<< "\t" << getTotalFitness(Pop, targetLbs, fitnessArr) / popSize << "\t" << count << endl;

        // closing output generation files
        genFile.close();
        
    }
    return 0;
}


// Get a random real number uniformly distributed between a, b

double randomdouble(double a, double b) {
    // Source of randomness for initializing random seed (http://stackoverflow.com/a/38245134)
    std::random_device rd;
    // Define a uniform distribution over reals
    std::uniform_real_distribution<double> unifDouble(a, b);
    // Mersenne twister pseudorandom number generator, initialized with a seed generated above
    std::mt19937 prng(rd());
    // Make a random draw from the distribution defined above
    double random = unifDouble(prng);
    
    return random;
}


// Get a random integer uniformly distributed between a, b

int randominteger(int a, int b) {
    std::random_device rd;
    // Define a uniform distribution over integers
    std::uniform_int_distribution<int> unifInt(a, b);
    std::mt19937 prng(rd());
    int randInt = unifInt(prng);
    
    return randInt;
}


/* Get a mutation list. Note that all the mutations are drawn in advance from zero-mean normal
 * distributions; when applying a mutation to a newly generated baby, we just assign it a value
 * from this pre-assembled list instead of actually mutating its parents' values. This function
 * takes a standard deviation as its second argument and assumes a mean of 0.
 */

std::vector<double> getMutList(int mutCount, const double traitStDev) {
    
    std::random_device randomnessSource;
    // define univariate normal distributions
    std::normal_distribution<double> lbMutMaker(0.0, traitStDev);
    // store draws from the normal distributions specified above using a for loop
    std::vector<double> mutList(mutCount, 0.0);
    for (int i = 0; i < mutCount; ++i) {
        std::mt19937 gen(randomnessSource());
        mutList[i] = lbMutMaker(gen);
    }
    
    return mutList;
}


// Randomly pick two parents from the current generation and make a baby from them

Individual pickAndMateParents(std::vector<Individual> &population, double totalFitness, double target, std::vector<double> &fitnessArr, std::vector<double> mutList) {
    // Pick parents
    int parentsPicked = 0;
    std::vector<int> parentIdx;

    while(parentsPicked < 2) {
        double tempFitness = 0;
        double rand = randomdouble(0, totalFitness);
        for(int i = 0; i < population.size(); i++) {
            tempFitness = tempFitness + fitnessArr[i];
            // cout << tempFitness << " " << rand << endl;
            if (tempFitness > rand) {
                parentIdx.push_back(i);
                parentsPicked++;
                break;
            }
        }
    }
    
    int parentAidx = parentIdx[0];
    int parentBidx = parentIdx[1];
    
    // The baby is initialized with no shelter
    Individual baby(population[parentAidx].phenotypeVal);
    int mutRate = 20000;
    int val;
    
    // Genotype inheritance
    for(int i = 0; i < baby.lbMutations.size(); i++) {
        val = randominteger(0, RAND_MAX);
        if(val % mutRate == 0) {
            baby.lbMutations[i] = mutList[i];
        } else if(val %2 == 0){
            baby.lbMutations[i] = population[parentAidx].lbMutations[i];
        } else{
            baby.lbMutations[i] = population[parentBidx].lbMutations[i];
        }
    }
    
    /* Shelter inheritance: if one of the two parents has a shelter, it will pass it on
     * to its offspring (and cease to have it, meaning that only the first offspring
     * the parent produces gets to inherit that parent's shelter). If both parents
     * have it, the one that will pass it on to the offspring is selected at random.
     */
    
    if (population[parentAidx].hasShelter() == true && population[parentBidx].hasShelter() == false) {
        baby.setShelter(true);
        population[parentAidx].setShelter(false);
    } else if (population[parentAidx].hasShelter() == false && population[parentBidx].hasShelter() == true) {
        baby.setShelter(true);
        population[parentBidx].setShelter(false);
    } else if (population[parentAidx].hasShelter() == true && population[parentBidx].hasShelter() == true) {
        baby.setShelter(true);
        double coinToss = randomdouble(0, 1);
        if (coinToss <= 0.5) {
            population[parentAidx].setShelter(false);
        } else {
            population[parentBidx].setShelter(false);
        }
    }

    return baby;
}


// Make sure that a baby satisfies the empirical constraints on trait values

bool checkBaby(Individual baby, double lowerLim, double upperLim) {

    double lb = baby.getLB();

    if (lb < lowerLim || lb > upperLim) {
        int stopper = 0;
        return false;
    }
    return true;
}


// Step the population forward by 1 generation

void evolvePop(vector<Individual> &population, double target, double lowerLim, double upperLim, int popSize, const std::vector<double> mutList, int numShelters, double calamFreq, double calamStrength, bool showShelterStats) {

    std::vector<Individual> reproductivePop;
    
    // Calamity check
    double randomDraw = randomdouble(0, 1);
    if (randomDraw < calamFreq) {
        std::vector<Individual> shelteredAndSurvivors;
        std::vector<Individual> unshelteredPop;
        int deathToll = ceil(calamStrength * population.size());
        
        // Divide the population into sheltered and unsheltered individuals
        for(int i = 0; i < population.size(); i++) {
            if (population[i].hasShelter() == true) {
                shelteredAndSurvivors.push_back(population[i]);
            } else {
                unshelteredPop.push_back(population[i]);
            }
        }
        
        // Get the indices of n = deathToll individuals to be eliminated from unshelteredPop
        std::vector<int> toEliminate;
        for(int j = 0; j < deathToll; j++) {
            int killed = randominteger(0, (int) (unshelteredPop.size() - 1));
            toEliminate.push_back(killed);
        }
        
        // Uncomment the line below for more verbose output:
        // cout << "   A calamity eliminated " << toEliminate.size() << " unsheltered individuals." << endl;
        
        // Add the survivors (= unshelteredPop members not in toEliminate) to the sheltered individuals
        for(int k = 0; k < unshelteredPop.size(); k++) {
            if (std::find(toEliminate.begin(), toEliminate.end(), k) == toEliminate.end()) {
                shelteredAndSurvivors.push_back(unshelteredPop[k]);
            }
        }
        
        // Reproductive population set to the union of sheltered individuals and unsheltered survivors
        reproductivePop = shelteredAndSurvivors;
    } else {
        // If no calamity has occurred, reproductive population = total population
        reproductivePop = population;
    }
    
    // Get total fitness
    std::vector<double> fitnessArr(reproductivePop.size(), 0.0);
    double totalFitness = getTotalFitness(reproductivePop, target, fitnessArr);
    
    // Loop for length of population, make babies
    int counter = 0;
    int stallCount = 0;
    Individual ind = reproductivePop[0];
    std::vector<Individual> newPop(popSize, ind);
    while(counter < popSize) {
        // Pick parents and make a baby
        Individual baby = pickAndMateParents(reproductivePop, totalFitness, target, fitnessArr, mutList);

        if (!checkBaby(baby, lowerLim, upperLim)) {
            stallCount++;
            if(stallCount > 500) {
                cout << "problem" << endl;
                stallCount = 0;
            }
        } else {
            newPop[counter] = baby;
            // Add baby to newPop
            counter++;
        }
    }
    
    // Check for unassigned shelters
    std::vector<int> shelteredIndividuals;
    for(int l = 0; l < popSize; l++) {
        if (newPop[l].hasShelter() == true) {
            shelteredIndividuals.push_back(l);
        }
    }
    int availableShelters = (int) (numShelters - shelteredIndividuals.size());
    
    // Print out the number of assigned and unassigned shelters
    if (showShelterStats) {
        cout << "   Assigned shelters: " << shelteredIndividuals.size() << endl;
        cout << "   Unassigned shelters: " << availableShelters << endl;
    }
    
    // Assign remaining shelters at random
    while(availableShelters > 0) {
        // Select an individual at random
        int ind = randominteger(0, popSize - 1);
        // Only assign a shelter to this individual if it does not have one already
        if (std::find(shelteredIndividuals.begin(), shelteredIndividuals.end(), ind) == shelteredIndividuals.end()) {
            newPop[ind].setShelter(true);
            // Decrement the number of available shelters
            availableShelters--;
        }
    }

    population = newPop;
    reproductivePop.clear();
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

string getData(vector<Individual> &population, double target) {

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

    s << lbVar << "\t";
    
    // getting average fitness
    std::vector<double> fitnessArr(popSize, 0.0);
    double avgFitness = getTotalFitness(population, target, fitnessArr) / popSize;
    
    s << avgFitness << "\n";

    string data = s.str();
    return data;
}
