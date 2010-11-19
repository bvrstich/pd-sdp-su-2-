/**
 * @mainpage 
 * This is an implementation of a spinsymmetrical primal dual interior point method
 * for optimizing the second order density matrix using the P Q G T1 and T2 N-representability conditions.
 * The method used is a path following algorithm with predictor corrector steps.
 * At compile time you can decide which condtions will be active compile with make PQ, PQG, PQGT1, PQGT2 or PQGT=(for all conditions).
 * @author Brecht Verstichel
 * @date 14-04-2010
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>

using std::cout;
using std::endl;
using std::ofstream;

#include "include.h"

/**
 * 
 * In the main the actual program is run.\n 
 * Part 1: An easy initial point is taken and then centered to the required precision (flag == 0)\n
 * Part 2: When the primal dual point is sufficiently centered steps are taken to reduce the
 * primal dual gap and take a large step in that direction (predictor) (flag == 1)\n
 * After each step a correcting step (flag == 2) is taken that brings the primal dual point closer to
 * the central path.\n
 * Part 3: When the primal dual gap is smaller that the required accuracy exit the while. (flag == 3)\n
 * For more information on the actual method, see primal_dual.pdf
 */

int main(int argc,char *argv[]){

   //initialize the random nr generator
   srand(time(NULL));

   cout.precision(10);

   int M = 32;//dim sp hilbert space
   int N = 16;//nr of particles

   char input_filename[100];
//   char brol[500];

   for(double g = 0.1;g < 5.001;g+= 0.1){

      sprintf(input_filename,"/home/bright/bestanden/results/sp_pairing/equidistant/8_16/PQG/DM_out/rdm_g=%f.out",g);

      //make the ifstream object
      ifstream input(input_filename);

      TPM tpm(M,N);

      tpm.in(input);

      TPM ham(M,N);
      ham.sp_pairing(g*M/2);

      cout << g << "\t" << tpm.ddot(ham) << endl;

   }

/* THIS IS FOR COUNTING NR OF LINES IN OUTPUT FILE
      //make the filename
      sprintf(input_filename,"/home/bright/bestanden/results/sp_pairing/equidistant/8_16/PQG/output/g=%f.out",g);

      //make the ifstream object
      ifstream input(input_filename);

      int iter = 0;

      while(input){
         input.getline(brol,500);
         ++iter;
      }

      cout << g << "\t" << iter << endl;

   }
*/

   return 0;

}
