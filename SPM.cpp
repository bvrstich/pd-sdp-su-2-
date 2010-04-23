#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * constructor, makes matrix of dimension M/2, there will be two degenerate blocks +1/2,-1/2. So
 * only one matrix is needed of dimension M/2 is needed to represent the SPM.
 * @param M dimension of single particle space
 * @param N Nr of particles
 */
SPM::SPM(int M,int N) : Matrix(M/2) {

   this->M = M;
   this->N = N;

}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(SPM &spm_copy) : Matrix(spm_copy) {

   this->M = spm_copy.gM();
   this->N = spm_copy.gN();

}

/**
 * TPM constructor: Creates a SPM initialized on the "bar" of the TPM.
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be initiated.
 */
SPM::SPM(double scale,TPM &tpm) : Matrix(tpm.gM()/2) {

   this->M = tpm.gM();
   this->N = tpm.gN();

   this->bar(scale,tpm);

}

/**
 * PHM constructor: Creates a SPM initialized on the "bar" of the PHM.
 * @param scale the factor u want the SPM to be scaled with
 * @param phm the PHM out of which the SPM will be initiated.
 */
SPM::SPM(double scale,PHM &phm) : Matrix(phm.gM()/2) {

   this->M = phm.gM();
   this->N = phm.gN();

   this->bar(scale,phm);

}

/**
 * destructor
 */
SPM::~SPM(){

}

/**
 * @return nr of particles
 */
int SPM::gN(){

   return N;

}

/**
 * @return dimension of sp space
 */
int SPM::gM(){

   return M;

}

ostream &operator<<(ostream &output,SPM &spm_p){

   for(int i = 0;i < spm_p.gn();++i)
      for(int j = 0;j < spm_p.gn();++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,TPM &tpm){

   //hulpvariabele
   double ward;

   for(int a = 0;a < M/2;++a)
      for(int c = a;c < M/2;++c){

         (*this)(a,c) = 0.0;

         for(int b = 0;b < M/2;++b){

            //S = 0 stuk
            ward = tpm(0,a,b,c,b);

            if(a == b)
               ward *= std::sqrt(2.0);

            if(c == b)
               ward *= std::sqrt(2.0);

            (*this)(a,c) += ward;

            //S = 1 stuk: hier kan nooit a = b en c = d wegens antisymmetrie
            (*this)(a,c) += 3.0*tpm(1,a,b,c,b);

         }

         //nog schalen
         (*this)(a,c) *= 0.5*scale;

      }

   this->symmetrize();

}

/**
 * Trace out a set of indices to create the "bar" matrix of a PHM, slight difference from the bar(TPM) function (normalization of the tp basisset).
 * @param scale the factor u want the SPM to be scaled with
 * @param phm the PHM out of which the SPM will be filled
 */
void SPM::bar(double scale,PHM &phm){

   for(int a = 0;a < M/2;++a)
      for(int c = a;c < M/2;++c){

         (*this)(a,c) = 0.0;

         for(int b = 0;b < M/2;++b){

            //S = 0 stuk
            for(int S = 0;S < 2;++S)
               (*this)(a,c) += phm.gdeg(S)*phm(S,a,b,c,b);

         }

         //nog schalen
         (*this)(a,c) *= 0.5*scale;

      }

   this->symmetrize();

}
