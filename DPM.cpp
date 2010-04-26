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

   }

   --counter;

}

/** 
 * Function that allocates and initializes the lists needed in the program, called when the first DPM object is constructed,
 * I think that this should actually be a static function, but for purely esthetic reasons.
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

   //if all the indices on either side are equal, the result is zero due to antisymmetry
   if( ( a == b && b == c ) || ( d == e && e == z ) )
      return 0;

   //if spin == 1/2
   if(S == 0){

      if( (a == b) || (a < b && b < c) ){//normal row index ordering

         int i = s2dp[0][S_ab][a][b][c];

         if( (d == e) || (d < e && e < z) ){//normal column ordering

            int j = s2dp[0][S_de][d][e][z];

            return (*this)(0,i,j);

         }
         else{//anomal column ordering in the morning

            return 0.0;

         }

      }
      else{//anomal row ordening in the morning

         return 0.0;

      }

   }
   else{//the spin == 3/2: this is a totally antisymmetrical (in the sp spatial orbitals) block

      //eerst kijken of er geen indices gelijk zijn:
      if(a == b || a == c || b == c)
         return 0;

      if(d == e || d == z || e == z)
         return 0;

      if(S_ab == 0 || S_de == 0)
         return 0;

      //dan kijken wel dp index met welke fase moet genomen worden:
      //eerst voor de i
      int i;

      int phase = 1;

      if(a < b){

         if(b < c)
            i = s2dp[1][0][a][b][c];
         else if(c < a)
            i = s2dp[1][0][c][a][b];
         else{

            i = s2dp[1][0][a][c][b];
            phase *= -1;

         }

      }
      else{

         if(a < c){

            i = s2dp[1][0][b][a][c];
            phase *= -1;

         }
         else if(c < b){

            i = s2dp[1][0][c][b][a];
            phase *= -1;

         }
         else
            i = s2dp[1][0][b][c][a];

      }

      //idem voor j maar met d e z
      int j;

      if(d < e){

         if(e < z)
            j = s2dp[1][0][d][e][z];
         else if(z < d)
            j = s2dp[1][0][z][d][e];
         else{

            j = s2dp[1][0][d][z][e];
            phase *= -1;

         }

      }
      else{

         if(d < z){

            j = s2dp[1][0][e][d][z];
            phase *= -1;

         }
         else if(z < e){

            j = s2dp[1][0][z][e][d];
            phase *= -1;

         }
         else
            j = s2dp[1][0][e][z][d];

      }

      return phase*(*this)(1,i,j);

   }

}
