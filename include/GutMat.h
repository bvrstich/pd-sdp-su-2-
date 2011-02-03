#ifndef GUTMAT_H
#define GUTMAT_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 03-02-2011\n\n
 * This class GutMat was written for the gutzwiller constraint matrices on singly occupied space. (see notes on hubbard).
 */

class GutMat : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param gm_p de GutMat you want to print
    */
   friend ostream &operator<<(ostream &output,GutMat &gm_p);

   public:
      
      //constructor
      GutMat(int M,int N);

      //copy constructor
      GutMat(GutMat &);

      //destructor
      virtual ~GutMat();

      using Matrix::operator=;

      int gN();

      int gM();

      void p(TPM &);

      void q(TPM &);

      void printeig();

   private:

      //!dimension of single particle space
      int M;

      //!nr of particles
      int N;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
