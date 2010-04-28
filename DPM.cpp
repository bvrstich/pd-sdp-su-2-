#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::endl;

#include "include.h"

int DPM::counter = 0;

int ***DPM::dp2s;
int *****DPM::s2dp;

double **DPM::_6j;

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 1/2 and 3/2.
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
DPM::DPM(int M,int N) : BlockMatrix(2) {

   this->N = N;
   this->M = M;

   //set the dimension and the degeneracies of the blocks
   this->setMatrixDim(0,M/2*(M/2 - 1) + M/2*(M/2 - 1)*(M/2 - 2)/3,2);
   this->setMatrixDim(1,M/2*(M/2 - 1)*(M/2 - 2)/6,4);

   if(counter == 0)//make the lists
      this->construct_lists();

   ++counter;

}

/**
 * copy constructor: constructs BlockMatrix object with two blocks, on for S=1/2 and one for S=3/2, and copies the content of the dpm_c blocks into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param dpm_c DPM to be copied into (*this)
 */
DPM::DPM(DPM &dpm_c) : BlockMatrix(dpm_c) {

   this->N = dpm_c.gN();
   this->M = dpm_c.gM();

   if(counter == 0)
      this->construct_lists();

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists dp2s en s2dp will be deleted.
 */
DPM::~DPM(){

   if(counter == 1){

      //first delete S = 1/2 part
      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < M/2;++a){

            for(int b = 0;b < M/2;++b)
               delete [] s2dp[0][S_ab][a][b];

            delete [] s2dp[0][S_ab][a];

         }

         delete [] s2dp[0][S_ab];

      }

      //then the S = 3/2 part
      for(int a = 0;a < M/2;++a){

         for(int b = 0;b < M/2;++b)
            delete [] s2dp[1][0][a][b];

         delete [] s2dp[1][0][a];

      }

      delete [] s2dp[1][0];

      for(int S = 0;S < 2;++S)
         delete [] s2dp[S];

      delete [] s2dp;

      //now delete dp2s 
      for(int S = 0;S < 2;++S){

         for(int i = 0;i < this->gdim(S);++i)
            delete [] dp2s[S][i];

               delete [] dp2s[S];

      }

      delete [] dp2s;

      for(int S = 0;S < 2;++S)
         delete [] _6j[S];

      delete [] _6j;

   }

   --counter;

}

/** 
 * Function that allocates and initializes the lists needed in the program, called when the first DPM object is constructed,
 */
void DPM::construct_lists(){

   //first allocation
   dp2s = new int ** [2];//two total spinblocks

   for(int S = 0;S < 2;++S){

      dp2s[S] = new int * [this->gdim(S)];//dimension of the blocks

      for(int i = 0;i < this->gdim(S);++i)
         dp2s[S][i] = new int [4];//amount of information stored for an index: (S_ab,a,b,c)

   }

   s2dp = new int **** [2];//two spinblocks

   s2dp[0] = new int *** [2];//for the S = 1/2, we have that S_ab can be 0 or 1

   for(int S_ab = 0;S_ab < 2;++S_ab){//loop and allocate

      s2dp[0][S_ab] = new int ** [M/2];

      for(int a = 0;a < M/2;++a){

         s2dp[0][S_ab][a] = new int * [M/2];

         for(int b = 0;b < M/2;++b)
            s2dp[0][S_ab][a][b] = new int [M/2];

      }

   }

   s2dp[1] = new int *** [1];//for the S = 3/2, we have that S_ab can be only 1

   s2dp[1][0] = new int ** [M/2];

   for(int a = 0;a < M/2;++a){//loop and allocate

      s2dp[1][0][a] = new int * [M/2];

      for(int b = 0;b < M/2;++b)
         s2dp[1][0][a][b] = new int [M/2];

   }

   //initialize the lists
   int teller = 0;

   //first S == 1/2, S_ab == 0 and a == b != c
   for(int a = 0;a < M/2;++a){

      for(int c = 0;c < a;++c){

         s2dp[0][0][a][a][c] = teller;

         dp2s[0][teller][0] = 0;//S_ab

         dp2s[0][teller][1] = a;
         dp2s[0][teller][2] = a;
         dp2s[0][teller][3] = c;

         ++teller;

      }

      for(int c = a + 1;c < M/2;++c){

         s2dp[0][0][a][a][c] = teller;

         dp2s[0][teller][0] = 0;//S_ab

         dp2s[0][teller][1] = a;
         dp2s[0][teller][2] = a;
         dp2s[0][teller][3] = c;

         ++teller;

      }

   }

   //S and S_ab the same but a != b != c
   for(int a = 0;a < M/2;++a)
      for(int b = a + 1;b < M/2;++b)
         for(int c = b + 1;c < M/2;++c){

            s2dp[0][0][a][b][c] = teller;

            dp2s[0][teller][0] = 0;//S_ab

            dp2s[0][teller][1] = a;
            dp2s[0][teller][2] = b;
            dp2s[0][teller][3] = c;

            ++teller;

         }

   //S == 0, S_ab == 1, a != b != c
   for(int a = 0;a < M/2;++a)
      for(int b = a + 1;b < M/2;++b)
         for(int c = b + 1;c < M/2;++c){

            s2dp[0][1][a][b][c] = teller;

            dp2s[0][teller][0] = 1;//S_ab

            dp2s[0][teller][1] = a;
            dp2s[0][teller][2] = b;
            dp2s[0][teller][3] = c;

            ++teller;

         }

   //re-init teller
   teller = 0;

   //S == 1, S_ab == 1, a != b != c
   for(int a = 0;a < M/2;++a)
      for(int b = a + 1;b < M/2;++b)
         for(int c = b + 1;c < M/2;++c){

            s2dp[1][0][a][b][c] = teller;

            dp2s[1][teller][0] = 1;//S_ab

            dp2s[1][teller][1] = a;
            dp2s[1][teller][2] = b;
            dp2s[1][teller][3] = c;

            ++teller;

         }

   //allocate
   _6j = new double * [2];

   for(int S = 0;S < 2;++S)
      _6j[S] = new double [2];

   //initialize
   _6j[0][0] = -0.5;
   _6j[0][1] = 0.5;
   _6j[1][0] = 0.5;
   _6j[1][1] = 1.0/6.0;

}

/**
 * @return number of particles
 */
int DPM::gN(){

   return N;

}

/**
 * @return number of single particle oribals
 */
int DPM::gM(){

   return M;

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymemtry relations are automatically accounted for:\n\n
 * DPM(S,S_ab,a,b,c,S_de,d,e,f) = sum_S_ac (some terms dependent on spin) DPM(S,S_ac,a,c,b,S_de,d,e,f) etc...
 * @param S The block index, when == 0 then access the block S = 1/2, for block == 1 we access the S = 3/2.
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param a first sp index that forms the dp row index i together with b, c and S_ab in block S
 * @param b second sp index that forms the dp row index i together with a, c and S_ab in block S
 * @param c third sp index that forms the dp row index i together with a, b and S_ab in block S
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param d first sp index that forms the dp column index j together with e, z and S_de in block S
 * @param e second sp index that forms the dp column index j together with d, z and S_de in block S
 * @param z third sp index that forms the dp column index j together with d, e and S_de in block S
 * @return the number on place DPM(S,i,j) with the right phase and forefactor.
 */
double DPM::operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   int *i = new int [2];
   double *coef_i = new double [2];

   int dim_i = get_inco(S,S_ab,a,b,c,i,coef_i);

   if(dim_i == 0){

      delete [] i;
      delete [] coef_i;

      return 0.0;

   }

   int *j = new int [2];
   double *coef_j = new double [2];

   int dim_j = get_inco(S,S_de,d,e,z,j,coef_j);

   if(dim_j == 0){

      delete [] i;
      delete [] j;

      delete [] coef_i;
      delete [] coef_j;

      return 0.0;

   }

   double ward = 0.0;

   for(int I = 0;I < dim_i;++I)
      for(int J = 0;J < dim_j;++J)
         ward += coef_i[I] * coef_j[J] * (*this)(S,i[I],j[J]);

   delete [] i;
   delete [] j;

   delete [] coef_i;
   delete [] coef_j;

   return ward;

}

/** 
 * Static member function that gets the dp-indices and their coefficients of the sp indices S,S_ab,a,b,c.
 * @param S block index of the state
 * @param S_ab intermediate spincoupling of a and b.
 * @param a first sp orbital
 * @param b second sp orbital
 * @param c third sp orbital
 * @param i pointer of dim 1 or 2 containing the indices occuring in the expansion of this particular dp state in the normal basis (a==b,c a < b < c).
 * @param coef pointer of dim 1 or 2 containing the coefficients occuring in the expansion.
 * @return the number of terms in the expansion (1 or 2), also the dim of pointers i and coef. When zero is returned this is not a valid element.
 */
int DPM::get_inco(int S,int S_ab,int a,int b,int c,int *i,double *coef){

   //they cannot all be equal
   if(a == b && b == c)
      return 0;

   if(S == 0){//spin 1/2 block:

      //if normal basis:
      if(a == b){

         if(S_ab == 1)//spin has to be zero for a == b
            return 0;

         i[0] = s2dp[0][0][a][b][c];
         coef[0] = 1;

         return 1;

      }
      else if (a < b && b < c){

         i[0] = s2dp[0][S_ab][a][b][c];
         coef[0] = 1;

         return 1;

      }
      else{//anomal basis:

         int min,max,phase;

         //first order a and b for code saving reasons
         if(a < b){

            min = a;
            max = b;

            phase = 1;

         }
         else{

            min = b;
            max = a;

            phase = 1 - 2*S_ab;

            if(c > max){//we still have one simple dim = 1 term left: b < a < c

               i[0] = s2dp[0][S_ab][b][a][c];
               coef[0] = phase;

               return 1;

            }

         }

         //now we have four possibilities left:
         //don't forget to multiply every result by phase to get the right a and b for min and max!
         // 1) c < min < max
         // 2) c == min < max
         // 3) min < c < max
         // 4) min < max == c
         if(c < min){//c < min < max

            //the S_ca == 0 part:
            i[0] = s2dp[0][0][c][min][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * _6j[0][S_ab];

            //the S_ca == 1 part:
            i[1] = s2dp[0][1][c][min][max];
            coef[1] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * std::sqrt(3.0) * _6j[1][S_ab];

            return 2;

         }
         else if(c == min){//c == min < max: this will also be a 1 dim list, because S_ac can only be 0 if a == c.

            i[0] = s2dp[0][0][c][min][max];
            coef[0] = std::sqrt(2.0) * phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * _6j[0][S_ab];

            return 1;

         }
         else if(c < max){//min < c < max

            //S_ac == 0 part:
            i[0] = s2dp[0][0][min][c][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * _6j[0][S_ab];

            //S_ac == 1 part:
            i[1] = s2dp[0][1][min][c][max];
            coef[1] = - phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * std::sqrt(3.0) * _6j[1][S_ab];

            return 2;

         }
         else{// min < c == max: also a 1 dim list, S_bc can only be 0 if b == c

            i[0] = s2dp[0][0][max][c][min];
            coef[0] = phase * std::sqrt(2.0) * std::sqrt(2.0*S_ab + 1.0) *_6j[0][S_ab];

            return 1;

         }

      }

   }
   else{//spin 3/2 block, totally antisymmetrical in the spatial sp orbs.

      //only S_ab == 1 can couple to 3/2's.
      if(S_ab == 0)
         return 0;

      //if any of the sp orbs are equal, antisymmetry leads to zero:
      if(a == b || b == c || c == a)
         return 0;

      if(a < b){

         if(b < c){//a < b < c

            i[0] = s2dp[1][0][a][b][c];
            coef[0] = 1;

         }
         else if(c < a){//c < a < b

            i[0] = s2dp[1][0][c][a][b];
            coef[0] = 1;

         }
         else{//a < c < b

            i[0] = s2dp[1][0][a][c][b];
            coef[0] = -1;

         }

      }
      else{//b < a

         if(a < c){//b < a < c

            i[0] = s2dp[1][0][b][a][c];
            coef[0] = -1;

         }
         else if(c < b){//c < b < a

            i[0] = s2dp[1][0][c][b][a];
            coef[0] = -1;

         }
         else{//b < c < a

            i[0] = s2dp[1][0][b][c][a];
            coef[0] = 1;

         }

      }

      return 1;

   }

}

/**
 * Print the uncoupled version of the DPM. Really only needed for debugging purposes, so can be inefficient.
 */
void DPM::uncouple(const char *filename){

   ofstream output(filename);

   output.precision(10);

   //first make table of cg-coefs needed:
   double cg[2][3][2][4][4];

   //init
   for(int i = 0;i < 2;++i)
      for(int j = 0;j < 3;++j)
         for(int k = 0;k < 2;++k)
            for(int l = 0;l < 4;++l)
               for(int m = 0;m < 4;++m)
                  cg[i][j][k][l][m] = 0.0;

   //(1/2 +1/2 1/2 -/1/2 | 0 0)
   cg[0][1][0][0][1] = 1.0/std::sqrt(2.0);

   //(1/2 -1/2 1/2 +/1/2 | 0 0)
   cg[0][0][1][0][1] = -1.0/std::sqrt(2.0);

   //(1/2 +1/2 1/2 -/1/2 | 1 0)
   cg[0][1][0][2][1] = 1.0/std::sqrt(2.0);

   //(1/2 -1/2 1/2 +/1/2 | 1 0)
   cg[0][0][1][2][1] = 1.0/std::sqrt(2.0);

   //(1/2 +1/2 1/2 +1/2 | 1 +1)
   cg[0][1][1][2][2] = 1.0;

   //(1/2 -1/2 1/2 -1/2 | 1 -1)
   cg[0][0][0][2][0] = 1.0;

   //(1 +1 1/2 +1/2 | 3/2 + 3/2)
   cg[1][2][1][3][3] = 1.0;

   //(1 -1 1/2 -1/2 | 3/2 - 3/2)
   cg[1][0][0][3][0] = 1.0;

   //(1 +1 1/2 -1/2 | 3/2 +1/2) and (1 -1 1/2 +1/2 | 3/2 -1/2)
   cg[1][2][0][3][2] = cg[1][0][1][3][1] = 1.0/std::sqrt(3.0);

   //(1 +1 1/2 -1/2 | 3/2 +1/2) and (1 -1 1/2 +1/2 | 3/2 -1/2)
   cg[1][2][0][3][2] = cg[1][0][1][3][1] = 1.0/std::sqrt(3.0);

   //(1 0 1/2 +1/2 | 3/2 +1/2 ) and (1 0 1/2 -1/2 | 3/2 -1/2)
   cg[1][1][1][3][2] = cg[1][1][0][3][1] = std::sqrt(2.0/3.0);

   //(1 +1 1/2 -1/2 | 1/2 +1/2) and (1 -1 1/2 +1/2 | 1/2 -1/2)
   cg[1][2][0][1][2] = cg[1][0][1][1][1] = std::sqrt(2.0/3.0);

   //(1 0 1/2 +1/2 | 1/2 +1/2) and (1 0 1/2 -1/2 | 1/2 -1/2)
   cg[1][1][1][1][2] = cg[1][1][0][1][1] = -1.0/std::sqrt(3.0);

   int a,b,c,d,e,z;
   int s_a,s_b,s_c,s_d,s_e,s_z;

   int M_ab,M_de;
   int M_abc,M_dez;

   double ward;

   //three row indices of the DPM
   for(int alpha = 0;alpha < M;++alpha)
      for(int beta = alpha + 1;beta < M;++beta)
         for(int gamma = beta + 1;gamma < M;++gamma){

            //watch it now, 0 is down, 1 is up.
            a = alpha/2;
            s_a = alpha%2;

            b = beta/2;
            s_b = beta%2;

            c = gamma/2;
            s_c = gamma%2;

            M_ab = s_a + s_b;
            M_abc = M_ab + s_c;

            //and the three column indices of the DPM
            for(int delta = alpha;delta < M;++delta)
               for(int epsilon = delta + 1;epsilon < M;++epsilon)
                  for(int zeta = epsilon + 1;zeta < M;++zeta){

                     d = delta/2;
                     s_d = delta%2;

                     e = epsilon/2;
                     s_e = epsilon%2;

                     z = zeta/2;
                     s_z = zeta%2;

                     M_de = s_d + s_e;
                     M_dez = M_de + s_z;

                     ward = 0.0;

                     //all can go through
                     if(M_abc == M_dez){

                        //for S = 1/2
                        if(M_abc != 0 && M_abc != 3){//if M = -3/2 or 3/2, then the S = 1/2 term is not there

                           //4 terms left: S_ab and S_de = 0 or 1 

                           //first S_ab == 0 and S_de == 0
                           if(M_ab == 1 && M_de == 1)//only contribution from the M_ab == M_de == 0 term (0 == 1 here because -1/2 == 0 and +1/2 == 1)
                              ward += cg[0][s_a][s_b][0][1] * cg[0][s_d][s_e][0][1] * (*this)(0,0,a,b,c,0,d,e,z);

                           //S_ab == 0 and S_de == 1
                           if(M_ab == 1)//no limit on M_de this time
                              ward += cg[0][s_a][s_b][0][1] * cg[0][s_d][s_e][2][M_de] * cg[1][M_de][s_z][1][M_dez] * (*this)(0,0,a,b,c,1,d,e,z);

                           //S_ab == 1 and S_de == 0
                           if(M_de == 1)//no limit on M_ab this time
                              ward += cg[0][s_a][s_b][2][M_ab] * cg[1][M_ab][s_c][1][M_abc] * cg[0][s_d][s_e][0][1] * (*this)(0,1,a,b,c,0,d,e,z);

                           //S_ab == 1 and S_de == 1, no constraints
                           ward += cg[0][s_a][s_b][2][M_ab] * cg[1][M_ab][s_c][1][M_abc] * cg[0][s_d][s_e][2][M_de] * cg[1][M_de][s_z][1][M_dez] 
                           
                              * (*this)(0,1,a,b,c,1,d,e,z);

                        }

                        //for S = 3/2: no constraints, but only the S_ab = S_de = 1 term will contribute
                        ward += cg[0][s_a][s_b][2][M_ab] * cg[1][M_ab][s_c][3][M_abc] * cg[0][s_d][s_e][2][M_de] * cg[1][M_de][s_z][3][M_dez] * (*this)(1,1,a,b,c,1,d,e,z);

                        //norms:
                        if(a == b)
                           ward *= std::sqrt(2.0);

                        if(d == e)
                           ward *= std::sqrt(2.0);

                     }

                     output << alpha << "\t" << beta << "\t" << gamma << "\t" << delta << "\t" << epsilon << "\t" << zeta << "\t" << ward << endl;

                  }

         }

}
