#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>


class Individual{


public:
    double logBodySize;
    double bequeathalProb;
    std::vector<double> lbMutations; // genotype underlying log body size
    std::vector<double> bpMutations; // genotype underlying bequeathal probability
    double mutStrength = 100;
    bool shelter;

    double getLB(void) {

        double lbsize = logBodySize;
        for(int i = 0; i < lbMutations.size(); i++) {
            lbsize += lbMutations[i] * mutStrength;
        }

        return lbsize;
    }
    
    double getBP(void) {
        
        double bprob = bequeathalProb;
        for(int i = 0; i < bpMutations.size(); i++) {
            bprob += bpMutations[i] * mutStrength;
        }
        
        return bprob;
    }
    
    double getFitness(double target) {
        double distance = pow((getLB() - target), 2) * -.5;
        return exp(distance);
    }
    
    bool hasShelter(void) {
        return shelter;
    }
    
    void setShelter(bool sShelter) {
        shelter = sShelter;
    }

    // Constructor sets shelter = false by default
    Individual(double sPhenotypes[]) {
        // Set the phenotype of this individual equal to the passed-in values
        logBodySize = sPhenotypes[0];
        bequeathalProb = sPhenotypes[1];
        lbMutations = std::vector<double> (1000, 0.0);
        bpMutations = std::vector<double> (1000, 0.0);
        shelter = false;
    };

};

     
 #endif // INDIVIDUAL_H
