//--uses minimal algorithm for SIMGMS implementation--//

#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "simint.h"

static struct simint_shell* shells = NULL; //array of basis set shells for SIMINT
static struct simint_multi_shellpair* shell_pair_data = NULL; //array of shell pair data 

static int* sp_shell = NULL; //array telling if given shell is L shell or not
static int* ksize = NULL; //array telling if given shell is L shell or not
static int* ksize_simint = NULL; //array telling if given shell is L shell or not
static int* kstart_simint = NULL; //array telling if given shell is L shell or not

int iold, jold;
int nshells, ishell, ishell_base;
int nshell_simint;

//--------------------------------//
//--temporary decompose function--//
//--------------------------------//
long long int decompose(long long int input) {
  long long int test1 = 1+8*input;
  double d_test1 = (double)test1;
  double test2 = sqrt(d_test1);
  double test3 = -1.0 + test2;
  double ret = ceil(test3/2.0);
  return (long long int)ret;
}

//---------------------//
//--initialize SIMINT--//
//---------------------//
void initialize_c()
{
  simint_init();

  shells = malloc(1*sizeof(struct simint_shell));
  shell_pair_data = malloc(1*sizeof(struct simint_multi_shellpair));
  
  sp_shell = malloc(1*sizeof(int));
  ksize = malloc(1*sizeof(int));
  ksize_simint = malloc(1*sizeof(int));
  kstart_simint = malloc(1*sizeof(int));

  iold = -1; jold = -1;
  ishell = 0; ishell_base = 0;
}

//-------------------//
//--SIMINT clean-up--//
//-------------------//
void finalize_c()
{
    //--Free remaining memory--//
    free(kstart_simint);
    free(ksize_simint);
    free(ksize);
    free(sp_shell);

    free(shell_pair_data);
    for (int i = 0; i != nshell_simint; ++i) simint_free_shell(&shells[i]);
    free(shells);

    //--Finalize the SIMINT library--//
    simint_finalize();
}

//--------------------------//
//--reset SIMINT variables--//
//--------------------------//
void reset_c()
{
  //--reset necessary variables--//
  iold = -1; jold = -1;
  ishell = 0; ishell_base = 0;
}

//--------------------------------//
//--Get info on basis set shells--//
//--------------------------------//
void get_julia_shell_info_c(struct shell* p_input)
{
  struct shell input = (*p_input);

  printf("ATOM ID: %lld\n", input.atom_id);
  printf("ATOM AM: %lld\n", input.am);
  printf("ATOM NBAS: %lld\n", input.nbas);
  printf("ATOM NPRIM: %lld\n", input.nprim);
  printf("ATOM POS: %lld\n", input.pos);

  printf("ATOM EXPONENTS:\n");
  for (int i = 0; i != input.nprim; ++i) {
    printf("%f\n",input.exponents[i]);
  }
  printf("\n");

  int coeff_count = input.nbas == 4 ? 2*input.nprim : input.nprim;
  printf("ATOM COEFFICIENTS:\n");
  for (int i = 0; i != coeff_count; ++i) {
    printf("%f\n",input.coefficients[i]);
  }
  printf("\n");
}

void get_simint_shell_info_c(long long int shell_num)
{
  struct simint_shell input = shells[shell_num];

  printf("ATOM AM: %lld\n", input.am);
  printf("ATOM NPRIM: %lld\n", input.nprim);
  printf("ATOM POS: %lld\n", shell_num);

  printf("ATOM COORD: %f, %f, %f\n:", input.x, input.y, input.z);

  printf("ATOM EXPONENTS:\n");
  for (int i = 0; i != input.nprim; ++i) {
    printf("%f\n",input.alpha[i]);
  }
  printf("\n");

  printf("ATOM COEFFICIENTS:\n");
  for (int i = 0; i != input.nprim; ++i) {
    printf("%f\n",input.coef[i]);
  }
  printf("\n");
}

//---------------------------------------------//
//--Translate JuliaChem shell to simint_shell--//
//---------------------------------------------//
void allocate_shell_array_c(long long int nshell,
  long long int t_nshell_simint)
{
  //--Some variable set-up for basis set translation--//
  //int sizeof_simint_shell = 2*sizeof(int); //account for am and nprim
  //sizeof_simint_shell += 3*sizeof(double); //account for x,y and z
  //sizeof_simint_shell += 2*sizeof(double*); //account for alpha and coeff
  //sizeof_simint_shell += sizeof(size_t); //account for memsize
  //sizeof_simint_shell += sizeof(void*); //account for ptr

  nshells = (int)nshell;
  int nshell_simint = (int)t_nshell_simint;

  //--Resize arrays--//
  shells = realloc(shells, nshell_simint*sizeof(struct simint_shell));
  shell_pair_data = realloc(shell_pair_data,(nshells*(nshells+1)/2)*sizeof(struct simint_multi_shellpair));
  //printf("SHELL PAIR DATA SIZE: %d\n", (nshells*(nshells+1)/2)*sizeof(struct simint_multi_shellpair));

  sp_shell = realloc(sp_shell, nshells*sizeof(int));
  ksize = realloc(ksize, nshells*sizeof(int));
  ksize_simint = realloc(ksize_simint, nshell_simint*sizeof(int));
  kstart_simint = realloc(kstart_simint, nshell_simint*sizeof(int));
}

void add_shell_c(struct shell* p_input)
{
  struct shell input = (*p_input);

  int sp = input.sp;

  sp_shell[ishell_base] = sp ? 1 : 0;
  ksize[ishell_base] = input.nbas;

  for (int isp = 0; isp < sp+1; ++isp) { //two iterations for L shells to split them into s and p components
    simint_initialize_shell(&shells[ishell]);

    shells[ishell].x = input.atom_center[0];
    shells[ishell].y = input.atom_center[1];
    shells[ishell].z = input.atom_center[2];

    shells[ishell].am = sp ? isp : (int)(input.am-1);
    shells[ishell].nprim = (int)input.nprim;

    simint_allocate_shell(shells[ishell].nprim, &shells[ishell]);

    int nprim = (int)input.nprim;
    for (int iprim = 0; iprim != nprim; ++iprim) {
      shells[ishell].alpha[iprim] = input.exponents[iprim];
      shells[ishell].coef[iprim] = input.coefficients[iprim+nprim*isp];
    };

    if (sp)
      ksize_simint[ishell] = isp ? 3 : 1;
    else
      ksize_simint[ishell] = input.nbas;

    kstart_simint[ishell] = ishell == 0 ? 0 : kstart_simint[ishell-1] + ksize_simint[ishell-1];

    ++ishell;
  }

  ++ishell_base;
}

//------------------------//
//--Normalize shell list--//
//------------------------//
void normalize_shells_c()
{
  simint_normalize_shells(nshell_simint, shells);
}

//----------------------------//
//-- One-electron integrals --//
//----------------------------//
void compute_overlap_c(long long int ash, long long int bsh, double* ovr) {
  int ncomputed = 0;
  ncomputed = simint_compute_overlap(&shells[ash-1], 
    &shells[bsh-1], ovr);
}

void compute_ke_c(long long int ash, long long int bsh, double* ke) {
  int ncomputed = 0;
  ncomputed = simint_compute_ke(&shells[ash-1], 
    &shells[bsh-1], ke);
}

void compute_nah_c(long long int ncenter, double* Z, double* x, 
  double* y, double* z, long long int ash, long long int bsh, double* ke) {

  int ncomputed = 0;
  ncomputed = simint_compute_potential(ncenter, Z, x, y, z, &shells[ash-1], 
    &shells[bsh-1], ke);
}

//------------------------------//
//--Precompute shell pair data--//
//------------------------------//
void precompute_shell_pair_data_c() {
  simint_initialize_multi_shellpairs(nshells*(nshells+1)/2, shell_pair_data);

  for (int sha = 0; sha != nshells; ++sha) {
    for (int shb = 0; shb <= sha; ++shb) { 
      int sh_idx = (sha*(sha+1)/2) + shb;
      //printf("ENTER PRECOMPUTE\n");
      //printf("%d, %d\n",sh_idx, nshells*(nshells+1)/2);
      //printf("%p\n",&shell_pair_data[sh_idx]);
   
      //simint_create_multi_shellpair(1, &shells[sha], 1, &shells[shb],
      //  &left_pair, 0);
      //printf("SHELL DATA SIZE: %d\n", sizeof(left_pair));
     
      simint_create_multi_shellpair(1, &shells[sha], 1, &shells[shb],
        &shell_pair_data[sh_idx], 0);
      //printf("EXIT PRECOMPUTE\n");
    }
  }
}

//-------------------------------//
//--Create and fill shell pairs--//
//-------------------------------//
/*
void create_ij_shell_pair_c(long long int ish, long long int jsh) {
  simint_create_multi_shellpair(1, &shells[ish-1], 1, &shells[jsh-1],
    &left_pair, 0);
}

void allocate_kl_shell_pair_c(long long int ksh, long long int lsh) {
  simint_allocate_multi_shellpair(1, &shells[ksh-1], 1, &shells[lsh-1],
    &right_pair, 0);
}

void create_kl_shell_pair_c(long long int ksh, long long int lsh) {
  simint_create_multi_shellpair(1, &shells[ksh-1], 1, &shells[lsh-1],
    &right_pair, 0);
}

void fill_kl_shell_pair_c(long long int ksh, long long int lsh) {
  simint_fill_multi_shellpair(1, &shells[ksh-1], 1, &shells[lsh-1],
    &right_pair, 0);
}
*/

//----------------//
//--Compute ERIs--//
//----------------//
void compute_eris_c(long long int ish, long long int jsh, long long int ksh,
  long long int lsh, double* eri, double* work) {
  
  //-- set up work buffers --//
  int ij_idx = (ish*(ish-1)/2) + jsh - 1;
  int kl_idx = (ksh*(ksh-1)/2) + lsh - 1;

  struct simint_multi_shellpair left_pair = shell_pair_data[ij_idx]; 
  struct simint_multi_shellpair right_pair = shell_pair_data[kl_idx]; 
 
#if 0 
  printf("IJ %d, %d, %d, %d:\n", ish, jsh, ksh, lsh);
  printf("%d, %d, %d\n",left_pair_.am1, left_pair_.am2, left_pair_.nprim);
  printf("%d, %d, %d\n",left_pair_.nshell12, left_pair_.nshell12_clip, *(left_pair_.nprim12));
  printf("%f, %f, %f\n",*(left_pair_.x), *(left_pair_.y), *(left_pair_.z));

  printf("KL %d, %d, %d, %d:\n", ish, jsh, ksh, lsh);
  printf("%d, %d, %d\n",right_pair_.am1, right_pair_.am2, right_pair_.nprim);
  printf("%d, %d, %d\n",right_pair_.nshell12, right_pair_.nshell12_clip, *(right_pair_.nprim12));
  printf("%f, %f, %f\n",*(right_pair_.x), *(right_pair_.y), *(right_pair_.z));
 #endif

  int ncomputed = 0;
  ncomputed = simint_compute_eri(&left_pair,
    &right_pair, 0.0, work, eri);
  
  //SIMINT_FREE(work);
}
