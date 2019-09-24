#include <simint/simint.h>

extern struct simint_multi_shellpair left_pair; //bra SIMINT shell pair structure
extern struct simint_multi_shellpair right_pair;  //ket SIMINT shell pair structure

#define MAX_CONTRACTION 12

extern struct shell {
  long long int atom_id;

  double exponents[MAX_CONTRACTION];
  double coefficients[2*MAX_CONTRACTION];

  double atom_center[3];

  long long int am;
  long long int nbas;
  long long int nprim;
  long long int pos;
};

#undef MAX_CONTRACTION

/*    
      The following serves as the documentation example for the SIMGMS routine.
      
      We will use the example from the ATOMS subroutine in inputa.src -
      CSiH with the 6-311G(D,P) basis. For just those three atoms, in that order,
      the NSHELL variables will look like this: 
              S  L  L  L  D    S  L  L  L  L  D    S  S  S  P
      KATOM   1  1  1  1  1    2  2  2  2  2  2    3  3  3  3
      KNG     6  3  1  1  1    6  6  3  1  1  1    3  1  1  1
      KTYPE   1  2  2  2  3    1  2  2  2  2  3    1  1  1  2
      KMIN    1  1  1  1  5    1  1  1  1  1  5    1  1  1  2
      KMAX    1  4  4  4 10    1  4  4  4  4 10    1  1  1  4
      KSTART  1  7 10 11 12   13 19 25 28 29 30   31 34 35 36 (SUM KNG)
      KLOC    1  2  6 10 14   20 21 25 29 33 37   43 44 45 46

      As SIMINT does not have L shells, SIMGMS "splits" each L shell into
      its S and P components. The basis set and SIMGMS arrays will then look 
      like this:
                      S | S  P | S  P | S  P | D    S | S  P | S  P | S  P | S  P | D    S  S  S  P
      ksize           1     4      4      4    6    1     4      4      4      4    6    1  1  1  3
      ksize_simint    1   1  3   1  3   1  3   6    1   1  3   1  3   1  3   1  3   6    1  1  1  3
      kstart_simint   1   2  3   6  7   10 11  14   20  21 22  25 26  29 30  33 34  37   43 44 45 46
      sp_shell        0     1      1      1    0    0     1      1      1      1    0    0  0  0  0
*/





//All code below this line is automatically generated
//-------------------------------------------------------------------//

void simgms_retrieve_eris_c_0000(int ii, int jj, int kk, int ll, double* eri); 

//void simgms_retrieve_eris_c_L(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes, int* L_);
