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

// Header file definining the class 'Individual'
#include "individual.h"


using namespace std;
namespace po = boost::program_options;

// Convenience functions
double randomdouble(double a, double b);
int randominteger(int a, int b);
double variance(std::vector<double> vals, double mean);

// Evolving a population and the traits of its constituent individuals
std::vector<std::vector<double> > getMutList(int mutCount, const double lbStDev, const double bpStDev);
Individual pickAndMateParents(std::vector<Individual> &population, double totalFitness, double target, std::vector<double> &fitnessArr, std::vector<std::vector<double> > mutList, bool evolveLB, bool evolveBP, int *parentA, int *parentB);
bool checkBaby(Individual baby, double lowerLim, double upperLim);
void evolvePop(vector<Individual> &population, double target, double selStrength, double lowerLim, double upperLim, int popSize, const std::vector<std::vector<double> > mutList, int numShelters, double calamFreq, double calamStrength, bool evolveLB, bool evolveBP, bool showShelterStats, int *uniqueParents);

// Extracting data
double getTotalFitness(std::vector<Individual> &population, double target, double selStrength, std::vector<double> &fitnessArr);
double getTraitMean(std::vector<Individual> &population, string trait);
double getTraitVar(std::vector<Individual> &population, string trait);
string getData(std::vector<Individual> &population, double target, double selStrength);


int main(int argc, char** argv)
{
    
    po::options_description desc("Simulate trait evolution under wealth inheritance");
    desc.add_options()
        ("help", "Show help message")
        ("mutCount", po::value<int>()->default_value(1000), "Number of mutations we will have")
        ("lbMutStrength", po::value<double>()->default_value(0.0005), "Standard deviation of the mutation model for log body size")
        ("bpMutStrength", po::value<double>()->default_value(0.005), "Standard deviation of the mutation model for bequeathal probability")
        ("popSize", po::value<int>()->default_value(100), "Number of individuals in the population")
        ("endpointsensitivity", po::value<double>()->default_value(0.02), "How close body size has to get to target value")
        ("reps", po::value<int>()->default_value(1), "Number of simulation repetitions")
        ("gen_limit", po::value<int>()->default_value(30000), "Maximum number of generations before cutting the simulation short")
        ("burnin", po::value<bool>()->default_value(true), "Should the simulation generate some initial variation first?")
        ("min_lb", po::value<double>()->default_value(2.302585), "Log of minimum body size")    // ln(10)
        ("max_lb", po::value<double>()->default_value(6.907755), "Log of maximum body size")    // ln(1000)
        ("start_lb", po::value<double>()->default_value(4.605171), "Log of starting body size") // ln(100)
        ("target_lb", po::value<double>()->default_value(5.010635), "Log of target body size")  // ln(150)
        ("selection_strength", po::value<double>()->default_value(65), "Penalty incurred by increasing distance from target body size")
        /* Default selection strength value: derived here by plugging an empirical value reported in a meta-analysis by Hoekstra et al.
         * (2001; doi:10.1073/pnas.161281098) -- 0.088, median standardized linear selection gradient for viability selection -- into
         * Eq. 1 of Alfaro et al. (2004; doi:10.1111/j.0014-3820.2004.tb01673.x). (1/2)*1/(0.088^2) = 64.566, which is rounded here to
         * the value above.
         */
        ("start_bp", po::value<double>()->default_value(0.05), "Starting bequeathal probability")
        ("shelterFraction", po::value<double>()->default_value(0.1), "Fraction of the population possessing shelter") // Set to 10% by default
        ("calamityFrequency", po::value<double>()->default_value(0.1), "Frequency of calamities affecting unsheltered individuals")
        ("calamityStrength", po::value<double>()->default_value(0.1), "Fraction of the population eliminated by a calamity")
	 // ("output_prefix", po::value<int>()->default_value(""), "Prefix for output files")
        ("output_suffix", po::value<string>()->default_value(""), "Prefix for output files");



    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);
    const int mutCount = vm["mutCount"].as<int>();
    const double bodySizeStDev = vm["lbMutStrength"].as<double>();
    const double beqProbStDev = vm["bpMutStrength"].as<double>();
    const int popSize = vm["popSize"].as<int>();
    const double endpointsensitivity = vm["endpointsensitivity"].as<double>();
    const int reps = vm["reps"].as<int>();
    const int gen_limit = vm["gen_limit"].as<int>();
    const bool burnin = vm["burnin"].as<bool>();
    const double min_lb = vm["min_lb"].as<double>();
    const double max_lb = vm["max_lb"].as<double>();
    const double start_lb = vm["start_lb"].as<double>();
    const double target_lb = vm["target_lb"].as<double>();
    const double selection_strength = vm["selection_strength"].as<double>();
    const double start_bp = vm["start_bp"].as<double>();
    const double shelterFraction = vm["shelterFraction"].as<double>();
    const double calamityFrequency = vm["calamityFrequency"].as<double>();
    const double calamityStrength = vm["calamityStrength"].as<double>();
 // const string output_prefix = vm["output_prefix"].as<string>();
    const string output_suffix = vm["output_suffix"].as<string>();

    // Output user-specified parameter values
    ofstream pFile("simulation_settings_" + output_suffix + ".txt");
    pFile << "popSize\tbodySizeTarget\tmutCount\tendpointSensitivity" << endl; //
    pFile << popSize << "\t" << exp(target_lb) << "\t" << mutCount << "\t" << endpointsensitivity << endl;
    pFile.close();

    // Output stream
    ofstream outputFile("endpoints_" + output_suffix + ".txt");
    outputFile << "averageBodySize\tbodySizeVariance\taverageBequeathalProb\tbequeathalProbVariance\taverageFitness\tgenerationsPassed" << endl;

    // Looping whole sim
    for(int i = 0; i < reps; i++) {
        srand (time(0));
        cout << " " << endl;
        cout << "simulation rep: " << output_suffix << " ------------------------------------------------------------------------" << endl;
        cout << " " << endl;
        
        // Generation files opening
        std::ostringstream params;
        params << "parameters_" << output_suffix << ".txt";
        string pName = params.str();
        ofstream genFile( pName.c_str() );
        genFile << "generation\taverageBodySize\tbodySizeVariance\taverageBequeathalProb\tbequeathalProbVariance\taverageFitness\tuniqueParents\n";

        // Get a mutation list
        std::vector<std::vector<double> > mutList;
        mutList = getMutList(mutCount, bodySizeStDev, beqProbStDev);

        /* Get phenotypes for starting populations. If burnin == true, we will start with a
         * population of individuals sharing the same randomly generated values, and let it evolve
         * (using the same evolvePop() function used for the simulation proper) toward a specified
         * starting body size. The idea is to generate a population where all individuals are
         * phenotypically identical and equally distant from the fitness optimum, but which already
         * contains some genotypic variation for selection to act upon. If burnin == false, the
         * simulation starts with a bunch of individuals that are identical both  phenotypically
         * and genotypically.
         */
         
        double startingPhenotype [2];
        
        if (burnin) {
            startingPhenotype[0] = randomdouble(min_lb, max_lb); // random initial log body size
        } else {
            startingPhenotype[0] = start_lb;
        }
        
        startingPhenotype[1] = start_bp;
        
        // Create a population of specified size, with all individuals initialized as shelter-less
        Individual squirrel(startingPhenotype);
        std::vector<Individual> Pop(popSize, squirrel);
        
        // How often (in terms of generations) to print out data
        int outputFrequency = 50;
        
        // Get the absolute number of shelters available to the population
        int nShelters = ceil(shelterFraction * popSize);

        /* Burn in the populations (optional). The burnin phase will end when the average log body
         * size gets acceptably close to the specified starting value.
         */
        int uniqueParents;

        if (burnin) {
            std::ostringstream tmp;
            tmp << "tmp_" << output_suffix << ".txt"; // generate names for temporary files to store burnin data in
            string tmpName = tmp.str();
            ofstream tmpFile( tmpName.c_str() );      // open the temporary files

            int burnCount = 0;
            while( fabs( exp(getTraitMean(Pop, "bodySize")) - exp(start_lb) ) > endpointsensitivity ) {
                
                string data = getData(Pop, target_lb, selection_strength);
                tmpFile << burnCount << "\t" << data; // the temporary files are structured the same way as the generation files
                evolvePop(Pop, start_lb, selection_strength, min_lb, max_lb, popSize, mutList, 0, 0, 0, true, false, false, &uniqueParents);
                tmpFile << uniqueParents << endl;
            
                if(burnCount % outputFrequency == 0) {
                    cout << "burnin generation: " << burnCount << endl;
                    cout << "Avg body size: " << exp(getTraitMean(Pop, "bodySize"));
                    cout << "   Avg bequeathal prob: " << getTraitMean(Pop, "beqProb");
                    std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness: " << getTotalFitness(Pop, start_lb, selection_strength, fitnessArr) / popSize << endl;
                    cout << " " << endl;
                }
                
                burnCount++;
            }
            
            tmpFile.close();                        // close the temporary files for writing
            
            // Credit: Joost Huizinga, https://stackoverflow.com/a/26372759
            
            std::string lastline = "";              // declare a string variable to store the last line of the temporary file in
            std::ifstream tmpf( tmpName.c_str() );  // open the temporary file for reading
            tmpf.seekg(0, std::ios_base::end);      // go to the last character of the temporary file
            char ch = ' ';                          // initial character not equal to newline
            while(ch != '\n') {
                tmpf.seekg(-2, std::ios_base::cur); // go two characters back
                if( (int)tmpf.tellg() <= 0 ) {
                    tmpf.seekg(0);                  // start of the last line
                    break;
                }
                tmpf.get(ch);                       // check the next character
            }
            
            std::getline(tmpf, lastline);           // store the result in the previously declared variable
            tmpf.close();                           // close the temporary file for reading
            remove( tmpName.c_str() );              // delete the temporary file
            
            /* We want to send the last line of the temporary file to the first post-header file of the generation file,
             * but before doing so, we will replace the generation number by 0.
             */
            std::string lastlineParsed = lastline.erase(0, lastline.find("\t") );
            genFile << "0\t" << lastlineParsed << endl;
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
            
            /* Outputing all the generation level data. If we ran burnin, we already got our generation 0
             * in the generation file, so we reset the counter here and start at 1. If not, we start at 0
             * and will include it in the file.
             */
            int printcount = (burnin ? count + 1 : count);
            
            if(!finished && printcount % outputFrequency == 0) {
                string data = getData(Pop, target_lb, selection_strength);
                genFile << printcount << "\t" << data;
            }

            // Checking if the population has finished evolving
            if( fabs( exp(getTraitMean(Pop, "bodySize")) - exp(target_lb) ) > endpointsensitivity ) {
                
                if (printcount % outputFrequency == 0) {
		   
                    evolvePop(Pop, target_lb, selection_strength, min_lb, max_lb, popSize, mutList, nShelters, calamityFrequency, calamityStrength, true, true, true, &uniqueParents);
                    genFile << uniqueParents << endl;
                    cout << "generation: " << printcount;
                    cout << "   Avg body size: " << exp(getTraitMean(Pop, "bodySize"));
                    cout << "   Avg bequeathal prob: " << getTraitMean(Pop, "beqProb");
                    std::vector<double> fitnessArr(popSize, 0.0);
                    cout << "   Avg fitness: " << getTotalFitness(Pop, target_lb, selection_strength, fitnessArr) / popSize << endl;
                    
                    // Warnings
                    if (exp(getTraitMean(Pop, "bodySize")) < 11.0) {
                        cout << "   Approaching lower bound on body size   ";
                    } else if (exp(getTraitMean(Pop, "bodySize")) > 990.0) {
                        cout << "   Approaching upper bound on body size   ";
                    }
                } else {
                    // Evolve quietly
                    evolvePop(Pop, target_lb, selection_strength, min_lb, max_lb, popSize, mutList, nShelters, calamityFrequency, calamityStrength, true, true, false, &uniqueParents);
                }
            } else {
                finished = true;
            }
            
            count++;
        }

        // Print final stats
        cout << "ending simulation" << endl;
        cout << " " << endl;
        cout << "final stats: " << endl;
        cout << "avg body size: " << exp(getTraitMean(Pop, "bodySize")) << "   avg bequeathal prob: " << getTraitMean(Pop, "beqProb") << "   generations: " << count << endl;
        
        // Print total time
        double diff = difftime(time(0), start);
        cout << "total time: " << diff << " seconds" << endl;

        // Print endpoint data
        std::vector<double> fitnessArr(popSize, 0.0);
        outputFile << exp(getTraitMean(Pop, "bodySize")) << "\t" << getTraitVar(Pop, "bodySize") << "\t" << getTraitMean(Pop, "beqProb") << "\t" << getTraitVar(Pop, "beqProb") << "\t" << getTotalFitness(Pop, target_lb, selection_strength, fitnessArr) / popSize << "\t" << count << endl;

        // Generation files closing
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
    double randDouble = unifDouble(prng);
    
    return randDouble;
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


// Get the variance of a vector of real numbers given a pre-calculated mean

double variance(std::vector<double> vals, double mean) {
    int length = vals.size();
    double var = 0.0;
    for(int i = 0; i < length; i++) {
        var += ((vals[i] - mean) * (vals[i] - mean)) / length;
    }

    return var;
}


/* Get a mutation list. Note that all the mutations are drawn in advance from zero-mean normal
 * distributions; when applying a mutation to a newly generated baby, we just assign it a value
 * from this pre-assembled list instead of actually mutating its parents' values. This function
 * assumes a mean of 0 and allows the two traits (log body size and bequeathal probability) to
 * have different standard deviations, passed in as the 2nd and 3rd argument, respectively.
 */

std::vector<std::vector<double> > getMutList(int mutCount, const double lbStDev, const double bpStDev) {
    
    std::random_device lbRandomnessSource;
    std::random_device bpRandomnessSource;
    // Define univariate normal distributions
    std::normal_distribution<double> lbMutMaker(0.0, lbStDev);
    std::normal_distribution<double> bpMutMaker(0.0, bpStDev);
    // Blank vector for initialization
    std::vector<double> zeros(2, 0.0);
    // Store draws from the normal distributions specified above using a for loop
    std::vector<std::vector<double> > mutList(mutCount, zeros);
    for (int i = 0; i < mutCount; ++i) {
        std::mt19937 lbGen(lbRandomnessSource());
        mutList[i][0] = lbMutMaker(lbGen);
        std::mt19937 bpGen(bpRandomnessSource());
        mutList[i][1] = bpMutMaker(bpGen);
    }
    
    return mutList;
}


// Randomly pick two parents from the current generation and make a baby from them

Individual pickAndMateParents(std::vector<Individual> &population, double totalFitness, double target, std::vector<double> &fitnessArr, std::vector<std::vector<double> > mutList, bool evolveLB, bool evolveBP, int *parentA, int *parentB) {
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
    double initPhenotypeVals [2] = {population[parentAidx].logBodySize, population[parentAidx].bequeathalProb};
    Individual baby(initPhenotypeVals);
    int mutRate = 20000;
    int val;
    
    // Genotype inheritance
    for(int i = 0; i < baby.lbMutations.size(); i++) {
        val = rand(); // using val = randominteger(0, RAND_MAX) causes a huge performance hit
        // double val = randomdouble(0, 1);
        // if (val < 1 / (double) mutRate) {
        if (val % mutRate == 0) {
            if (evolveLB) {
                baby.lbMutations[i] = mutList[i][0];
            }
            if (evolveBP) {
                baby.bpMutations[i] = mutList[i][1];
            }
        // } else if (val < 0.5) {
        } else if (val %2 == 0) {
            baby.lbMutations[i] = population[parentAidx].lbMutations[i];
            baby.bpMutations[i] = population[parentAidx].bpMutations[i];
        } else {
            baby.lbMutations[i] = population[parentBidx].lbMutations[i];
            baby.bpMutations[i] = population[parentBidx].bpMutations[i];
        }
    }
    
    /* Shelter inheritance: if a parent has a shelter, it will pass it on to its
     * offspring with a probability equal to that parent's bequeathalProb attribute,
     * and cease to have it (meaning that only the first offspring the parent produces
     * gets to inherit that parent's shelter). If both parents want to pass their
     * shelter on to their offspring, parent A will do so with a probability of
     * A.getBP() / (A.getBP() + B.getBP()) and parent B will do so with a probability
     * of B.getBP() / (A.getBP() + B.getBP()).
     */
    
    double randA = randomdouble(0, 1);
    double randB = randomdouble(0, 1);
    double probA = population[parentAidx].getBP();
    double probB = population[parentBidx].getBP();
    
    if (population[parentAidx].hasShelter() == true && population[parentBidx].hasShelter() == false) {
        if (randA <= probA) {
            baby.setShelter(true);
            population[parentAidx].setShelter(false);
        }
    } else if (population[parentAidx].hasShelter() == false && population[parentBidx].hasShelter() == true) {
        if (randB <= probB) {
            baby.setShelter(true);
            population[parentBidx].setShelter(false);
        }
    } else if (population[parentAidx].hasShelter() == true && population[parentBidx].hasShelter() == true) {
        if (randA <= probA && randB > probB) {
            baby.setShelter(true);
            population[parentAidx].setShelter(false);
        } else if (randA > probA && randB <= probB) {
            baby.setShelter(true);
            population[parentBidx].setShelter(false);
        } else if (randA <= probA && randB <= probB) {
            baby.setShelter(true);
            double coinToss = randomdouble(0, 1);
            if (coinToss <= (probA / (probA + probB))) {
                population[parentAidx].setShelter(false);
            } else {
                population[parentBidx].setShelter(false);
            }
        }
    }
    *parentA = parentAidx;
    *parentB = parentBidx;
    return baby;
}


// Make sure that a baby satisfies the empirical constraints on trait values

bool checkBaby(Individual baby, double lowerLim, double upperLim) {

    double lb = baby.getLB();
    double bp = baby.getBP();

    // Body size constraints are user-specified; bequeathal probability has to be in [0, 1]
    if ((lb < lowerLim || lb > upperLim) || (bp < 0.0 || bp > 1.0)) {
        return false;
    }
    
    return true;
}


// Step the population forward by 1 generation

void evolvePop(vector<Individual> &population, double target, double selStrength, double lowerLim, double upperLim, int popSize, const std::vector<std::vector<double> > mutList, int numShelters, double calamFreq, double calamStrength, bool evolveLB, bool evolveBP, bool showShelterStats, int *uniqueParents) {
    std::vector<int> parentIndices;
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
    double totalFitness = getTotalFitness(reproductivePop, target, selStrength, fitnessArr);
    
    // Loop for length of population, make babies
    int counter = 0;
    int stallCount = 0;
    Individual ind = reproductivePop[0];
    std::vector<Individual> newPop(popSize, ind);
    while(counter < popSize) {
        int parentA, parentB;
        // Pick parents and make a baby
        Individual baby = pickAndMateParents(reproductivePop, totalFitness, target, fitnessArr, mutList, evolveLB, evolveBP, &parentA, &parentB);

        if (!checkBaby(baby, lowerLim, upperLim)) {
            stallCount++;
            if (stallCount > 500) {
                cout << "problem" << endl;
                stallCount = 0;
            }
        } else {
            newPop[counter] = baby;
            parentIndices.push_back(parentA);
     		parentIndices.push_back(parentB);
            // Add baby to newPop
            counter++;
        }
    }
    *uniqueParents = std::set<int>( parentIndices.begin(), parentIndices.end() ).size();
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
        cout << "   Unique parents: " << std::set<int>( parentIndices.begin(), parentIndices.end() ).size() << endl;
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

double getTotalFitness(std::vector<Individual> &population, double target, double selStrength, std::vector<double> &fitnessArr) {
    double totalFitness = 0;
    for(int i = 0; i < population.size(); i++) {
        double currentFitness = population[i].getFitness(target, selStrength);
        fitnessArr[i] = currentFitness;
        totalFitness += currentFitness;
    }
    return totalFitness;
}


// Return the average log(body size)

double getTraitMean(std::vector<Individual> &population, string trait) {
    double traitSum = 0;
    if (trait == "bodySize") {
        for(int i = 0; i < population.size(); i++) {
            traitSum += population[i].getLB();
        }
    } else if (trait == "beqProb") {
        for(int i = 0; i < population.size(); i++) {
            traitSum += population[i].getBP();
        }
    } else {
        throw std::invalid_argument("Argument 'trait' must be one of the following: bodySize, beqProb");
    }
    return traitSum/population.size();
}


// Get the total variation (sum of variances) of a trait in the population

double getTraitVar(std::vector<Individual> &population, string trait) {

    int popSize = population.size();
    // creating a vector to hold the trait values
    std::vector<double> popVect(popSize, 0.0);
  
    // enter the trait values into the vector
    if (trait == "bodySize") {
        for(int i = 0; i < popSize; i++) {
            popVect[i] = population[i].getLB();
        }
    } else if (trait == "beqProb") {
        for(int i = 0; i < popSize; i++) {
            popVect[i] = population[i].getBP();
        }
    } else {
        throw std::invalid_argument("Argument 'trait' must be one of the following: bodySize, beqProb");
    }

    double traitMean = (std::accumulate(popVect.begin(), popVect.end(), 0.0)) / popSize; // body size mean
    double traitVar = variance(popVect, traitMean);                                      // body size variance

    return traitVar;
}


// Get the population means and variances of individual traits

string getData(std::vector<Individual> &population, double target, double selStrength) {

    std::stringstream s;
    int popSize = population.size();
    // creating a vector of vectors to hold all the phenotype values
    std::vector<double> rowinit(popSize, 0.0);
    std::vector<std::vector<double> > popVects(2, rowinit);
  
    for(int i = 0; i < popSize; i++) {
        popVects[0][i] = population[i].getLB();
        popVects[1][i] = population[i].getBP();
    }

    double lbMean = (std::accumulate(popVects[0].begin(), popVects[0].end(), 0.0)) / popSize; // body size mean
    double lbVar = variance(popVects[0], lbMean);                                             // body size variance
    double bpMean = (std::accumulate(popVects[1].begin(), popVects[1].end(), 0.0)) / popSize; // bequeathal prob. mean
    double bpVar = variance(popVects[1], bpMean);                                             // bequeathal prob. variance
    
    // getting average fitness
    std::vector<double> fitnessArr(popSize, 0.0);
    double avgFitness = getTotalFitness(population, target, selStrength, fitnessArr) / popSize;
    
    s << lbMean << "\t" << lbVar << "\t" << bpMean << "\t" << bpVar << "\t" << avgFitness << "\t";
    string data = s.str();
    return data;
}
