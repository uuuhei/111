#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>

//需要重新制造一个head！
class Individual{


public:
      double phenotypeVal; //注意重新define过了 不再是arr了
      std::vector<double> lbMutations;
      double mutStrength = 100;

//得到body size的第一个trait value，定义了function

        double getLB(void)
        {

            double lbsize = phenotypeVal;
            for(int i =0; i < lbMutations.size(); i++){
                    lbsize+=lbMutations[i]*mutStrength;
                }

            return lbsize;
        };

//通过预计的trait的数值，得到这个个体的fitness
        double getFitness(double target)
        {
            double distance = pow((getLB() - target), 2) * -.5;
            return exp(distance);
        }


        Individual(double sPhenotypes)
        {
            //setting the phenotypes of this individual equal to the passed in phenotype vals
            phenotypeVal = sPhenotypes; //这里去掉了s！
            lbMutations = std::vector<double> (1000, 0.0);

        };

};

     
 #endif // INDIVIDUAL_H
