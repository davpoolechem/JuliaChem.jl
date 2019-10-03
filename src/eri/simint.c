//--uses minimal algorithm for SIMGMS implementation--//

#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <assert.h>

#include "simint.h"

static double* buffer = NULL; //shared workspace for SIMINT ERI computations
static double* work = NULL; //shared workspace for SIMINT ERI computations
static struct simint_shell* shells = NULL; //array of basis set shells for SIMINT

static int* sp_shell = NULL; //array telling if given shell is L shell or not
static int* ksize = NULL; //array telling if given shell is L shell or not
static int* ksize_simint = NULL; //array telling if given shell is L shell or not
static int* kstart_simint = NULL; //array telling if given shell is L shell or not

struct simint_multi_shellpair left_pair; //bra SIMINT shell pair structure
struct simint_multi_shellpair right_pair; //ket SIMINT shell pair structure

int nshells, ishell, ishell_base;
int nshell_simint;

//---------------------//
//--initialize SIMINT--//
//---------------------//
void initialize_c()
{
  simint_init();

  simint_initialize_multi_shellpair(&left_pair);
  simint_initialize_multi_shellpair(&right_pair);

  shells = malloc(1*sizeof(struct simint_shell));

  sp_shell = malloc(1*sizeof(int));
  ksize = malloc(1*sizeof(int));
  ksize_simint = malloc(1*sizeof(int));
  kstart_simint = malloc(1*sizeof(int));

  work = malloc(simint_ostei_workmem(0,1)*sizeof(double));
  buffer = malloc(81*sizeof(double));

  ishell = 0; ishell_base = 0;
}

//-------------------//
//--SIMINT clean-up--//
//-------------------//
void finalize_c()
{
    //--Free remaining memory--//
    SIMINT_FREE(buffer);
    SIMINT_FREE(work);

    free(kstart_simint);
    free(ksize_simint);
    free(ksize);
    free(sp_shell);

    for (int i = 0; i != nshell_simint; ++i) simint_free_shell(&shells[i]);
    free(shells);

    simint_free_multi_shellpair(&right_pair);
    simint_free_multi_shellpair(&left_pair);

    //--Finalize the SIMINT library--//
    simint_finalize();
}

//--------------------------//
//--reset SIMINT variables--//
//--------------------------//
void reset_c()
{
    //--reset necessary variables--//
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

  sp_shell = realloc(sp_shell, nshells*sizeof(int));
  ksize = realloc(ksize, nshells*sizeof(int));
  ksize_simint = realloc(ksize_simint, nshell_simint*sizeof(int));
  kstart_simint = realloc(kstart_simint, nshell_simint*sizeof(int));

  work = realloc(work, simint_ostei_workmem(0,1)*sizeof(double));
  buffer = realloc(buffer, 1296*sizeof(double));
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

//----------------------//
//--Create shell pairs--//
//----------------------//
void create_ij_shell_pair_c(long long int ish, long long int jsh) {
  simint_create_multi_shellpair(1, &shells[ish-1], 1, &shells[jsh-1],
    &left_pair, 0);
}

void create_kl_shell_pair_c(long long int ksh, long long int lsh) {
  simint_create_multi_shellpair(1, &shells[ksh-1], 1, &shells[lsh-1],
    &right_pair, 0);
}

//----------------//
//--Compute ERIs--//
//----------------//
void compute_eris_c(double* eri) {
  int ncomputed = 0;
  ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, eri);
}

//-------------------------------------------------------------//
//--Copy list of ERIs of the form (ish jsh|ksh lsh) to eri--//
//-------------------------------------------------------------//
/*
double* retrieve_eris_c(long long int ish, long long int jsh,
  long long int ksh, long long int lsh, double* eri)
{
  //--initialize some variables--//
  //int ii = (int)ish-1, jj = (int)jsh-1, kk = (int)ksh-1, ll = (int)lsh-1; //account for 1-indexing of ish,jsh

  //--start ERI computation--//
  bool new_ij = ish != iold || jsh != jold;
  if (new_ij) {
    simint_create_multi_shellpair(1, &shells[ish-1], 1, &shells[jsh-1], &left_pair, 0);
    iold = ish; jold = jsh;
  }

  simint_create_multi_shellpair(1, &shells[ksh-1], 1, &shells[lsh-1], &right_pair, 0);
  //printf("IJ %d, %d, %d, %d:\n", ii, jj, kk, ll);
  //printf("%d, %d, %d\n",left_pair.am1, left_pair.am2, left_pair.nprim);
  //printf("%d, %d, %d\n",left_pair.nshell12, left_pair.nshell12_clip, *(left_pair.nprim12));
  //printf("%f, %f, %f\n",*(left_pair.x), *(left_pair.y), *(left_pair.z));

  //printf("KL %d, %d, %d, %d:\n", ii, jj, kk, ll);
  //printf("%d, %d, %d\n",right_pair.am1, right_pair.am2, right_pair.nprim);
  //printf("%d, %d, %d\n",right_pair.nshell12, right_pair.nshell12_clip, *(right_pair.nprim12));
  //printf("%f, %f, %f\n",*(right_pair.x), *(right_pair.y), *(right_pair.z));
  int ncomputed = 0;
  ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, eri);
  //printf("%f\n",eri[0]);

  //for (int i = 0; i != ish-1; ++i) ii += sp_shell[i]; //account for splitting of GAMESS L shells into s and p SIMINT shells
  //for (int j = 0; j != jsh-1; ++j) jj += sp_shell[j];
  //for (int k = 0; k != ksh-1; ++k) kk += sp_shell[k];
  //for (int l = 0; l != lsh-1; ++l) ll += sp_shell[l];

  //int fullsizes[4] = { ksize[ish-1], ksize[jsh-1], ksize[ksh-1], ksize[lsh-1] };
  //int L_[4] = { sp_shell[ish-1], sp_shell[jsh-1], sp_shell[ksh-1], sp_shell[lsh-1] };

  //printf("START: %d, %d, %d, %d\n", kstart_simint[0], kstart_simint[1], kstart_simint[2], kstart_simint[3]);
  //--start ERI computation--//
  //if (L_[0] == 0 && L_[1] == 0 && L_[2] == 0 && L_[3] == 0)
  //  simgms_retrieve_eris_c_0000(ii, jj, kk, ll, eri);
  //else
  //  simgms_retrieve_eris_c_L(ii, jj, kk, ll, eri, fullsizes, L_);

  return eri;
}

//--------------------------------------------------------------------------------------------------//
//--Codes for calculating ERIs for shell quartets with different numbers and positions of L shells--//
//--------------------------------------------------------------------------------------------------//
//For fxn simgms_retrieve_eris_c_ijkl, any of shells i,j,k,l = 1 means that shell is an L shell
//Ex: simgms_retrieve_eris_c_0010 handles shell quartets where the third shell is an L shell

/*
void simgms_retrieve_eris_c_0000(int ii, int jj, int kk, int ll, double* eri) {

  int ncomputed = 0;

  //--start ERI computation--//
  bool new_ij = ii != iold || jj != jold;
  if (new_ij) {
    simint_create_multi_shellpair(1, &shells[ii], 1, &shells[jj], &left_pair, 0);
    iold = ii; jold = jj;
  }

  simint_create_multi_shellpair(1, &shells[kk], 1, &shells[ll], &right_pair, 0);
  //printf("IJ %d, %d, %d, %d:\n", ii, jj, kk, ll);
  //printf("%d, %d, %d\n",left_pair.am1, left_pair.am2, left_pair.nprim);
  //printf("%d, %d, %d\n",left_pair.nshell12, left_pair.nshell12_clip, *(left_pair.nprim12));
  //printf("%f, %f, %f\n",*(left_pair.x), *(left_pair.y), *(left_pair.z));

  //printf("KL %d, %d, %d, %d:\n", ii, jj, kk, ll);
  //printf("%d, %d, %d\n",right_pair.am1, right_pair.am2, right_pair.nprim);
  //printf("%d, %d, %d\n",right_pair.nshell12, right_pair.nshell12_clip, *(right_pair.nprim12));
  //printf("%f, %f, %f\n",*(right_pair.x), *(right_pair.y), *(right_pair.z));

  ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, eri);
  //printf("%f\n",eri[0]);
}

void simgms_retrieve_eris_c_L(int ii, int jj, int kk, int ll, double* eri, int* fullsizes, int* L_) {
    int key = 8*L_[0] + 4*L_[1] + 2*L_[2] + 1*L_[3];
    switch(key) {
        case 1: simgms_retrieve_eris_c_L_0001(ii, jj, kk, ll, eri, fullsizes); break;
        case 2: simgms_retrieve_eris_c_L_0010(ii, jj, kk, ll, eri, fullsizes); break;
        case 3: simgms_retrieve_eris_c_L_0011(ii, jj, kk, ll, eri, fullsizes); break;
        case 4: simgms_retrieve_eris_c_L_0100(ii, jj, kk, ll, eri, fullsizes); break;
        case 5: simgms_retrieve_eris_c_L_0101(ii, jj, kk, ll, eri, fullsizes); break;
        case 6: simgms_retrieve_eris_c_L_0110(ii, jj, kk, ll, eri, fullsizes); break;
        case 7: simgms_retrieve_eris_c_L_0111(ii, jj, kk, ll, eri, fullsizes); break;
        case 8: simgms_retrieve_eris_c_L_1000(ii, jj, kk, ll, eri, fullsizes); break;
        case 9: simgms_retrieve_eris_c_L_1001(ii, jj, kk, ll, eri, fullsizes); break;
        case 10: simgms_retrieve_eris_c_L_1010(ii, jj, kk, ll, eri, fullsizes); break;
        case 11: simgms_retrieve_eris_c_L_1011(ii, jj, kk, ll, eri, fullsizes); break;
        case 12: simgms_retrieve_eris_c_L_1100(ii, jj, kk, ll, eri, fullsizes); break;
        case 13: simgms_retrieve_eris_c_L_1101(ii, jj, kk, ll, eri, fullsizes); break;
        case 14: simgms_retrieve_eris_c_L_1110(ii, jj, kk, ll, eri, fullsizes); break;
        case 15: simgms_retrieve_eris_c_L_1111(ii, jj, kk, ll, eri, fullsizes); break;
    }
}

#define L_LOOP_BOUND(X) ((X) ? (3) : (1))

void simgms_retrieve_eris_c_L_0001(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 0; ++isp) {
    for (int jsp = 0; jsp <= 0; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 0; ++ksp) {
        for (int lsp = 0; lsp <= 1; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          int m_loop_bound = ksize_simint[ii+isp];
          int n_loop_bound = ksize_simint[jj+jsp];
          int o_loop_bound = ksize_simint[kk+ksp];
          //int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < m_loop_bound; ++m)
          for(int n = 0; n < n_loop_bound; ++n)
          for(int o = 0; o < o_loop_bound; ++o)
          for(int p = 0; p < L_LOOP_BOUND(lsp); ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_0010(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 0; ++isp) {
    for (int jsp = 0; jsp <= 0; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 1; ++ksp) {
        for (int lsp = 0; lsp <= 0; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          int m_loop_bound = ksize_simint[ii+isp];
          int n_loop_bound = ksize_simint[jj+jsp];
          //int o_loop_bound = ksize_simint[kk+ksp];
          int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < m_loop_bound; ++m)
          for(int n = 0; n < n_loop_bound; ++n)
          for(int o = 0; o < L_LOOP_BOUND(ksp); ++o)
          for(int p = 0; p < p_loop_bound; ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_0011(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 0; ++isp) {
    for (int jsp = 0; jsp <= 0; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 1; ++ksp) {
        for (int lsp = 0; lsp <= 1; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          int m_loop_bound = ksize_simint[ii+isp];
          int n_loop_bound = ksize_simint[jj+jsp];
          //int o_loop_bound = ksize_simint[kk+ksp];
          //int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < m_loop_bound; ++m)
          for(int n = 0; n < n_loop_bound; ++n)
          for(int o = 0; o < L_LOOP_BOUND(ksp); ++o)
          for(int p = 0; p < L_LOOP_BOUND(lsp); ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_0100(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 0; ++isp) {
    for (int jsp = 0; jsp <= 1; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 0; ++ksp) {
        for (int lsp = 0; lsp <= 0; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          int m_loop_bound = ksize_simint[ii+isp];
          //int n_loop_bound = ksize_simint[jj+jsp];
          int o_loop_bound = ksize_simint[kk+ksp];
          int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < m_loop_bound; ++m)
          for(int n = 0; n < L_LOOP_BOUND(jsp); ++n)
          for(int o = 0; o < o_loop_bound; ++o)
          for(int p = 0; p < p_loop_bound; ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_0101(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 0; ++isp) {
    for (int jsp = 0; jsp <= 1; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 0; ++ksp) {
        for (int lsp = 0; lsp <= 1; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          int m_loop_bound = ksize_simint[ii+isp];
          //int n_loop_bound = ksize_simint[jj+jsp];
          int o_loop_bound = ksize_simint[kk+ksp];
          //int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < m_loop_bound; ++m)
          for(int n = 0; n < L_LOOP_BOUND(jsp); ++n)
          for(int o = 0; o < o_loop_bound; ++o)
          for(int p = 0; p < L_LOOP_BOUND(lsp); ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_0110(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 0; ++isp) {
    for (int jsp = 0; jsp <= 1; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 1; ++ksp) {
        for (int lsp = 0; lsp <= 0; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          int m_loop_bound = ksize_simint[ii+isp];
          //int n_loop_bound = ksize_simint[jj+jsp];
          //int o_loop_bound = ksize_simint[kk+ksp];
          int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < m_loop_bound; ++m)
          for(int n = 0; n < L_LOOP_BOUND(jsp); ++n)
          for(int o = 0; o < L_LOOP_BOUND(ksp); ++o)
          for(int p = 0; p < p_loop_bound; ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_0111(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 0; ++isp) {
    for (int jsp = 0; jsp <= 1; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 1; ++ksp) {
        for (int lsp = 0; lsp <= 1; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          int m_loop_bound = ksize_simint[ii+isp];
          //int n_loop_bound = ksize_simint[jj+jsp];
          //int o_loop_bound = ksize_simint[kk+ksp];
          //int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < m_loop_bound; ++m)
          for(int n = 0; n < L_LOOP_BOUND(jsp); ++n)
          for(int o = 0; o < L_LOOP_BOUND(ksp); ++o)
          for(int p = 0; p < L_LOOP_BOUND(lsp); ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_1000(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 1; ++isp) {
    for (int jsp = 0; jsp <= 0; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 0; ++ksp) {
        for (int lsp = 0; lsp <= 0; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          //int m_loop_bound = ksize_simint[ii+isp];
          int n_loop_bound = ksize_simint[jj+jsp];
          int o_loop_bound = ksize_simint[kk+ksp];
          int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < L_LOOP_BOUND(isp); ++m)
          for(int n = 0; n < n_loop_bound; ++n)
          for(int o = 0; o < o_loop_bound; ++o)
          for(int p = 0; p < p_loop_bound; ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //printf("%s\n", ntotal-ncomputed+buffer_idx == eri_idx ? "true" : "false");
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_1001(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 1; ++isp) {
    for (int jsp = 0; jsp <= 0; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 0; ++ksp) {
        for (int lsp = 0; lsp <= 1; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0

          //int m_loop_bound = ksize_simint[ii+isp];
          int n_loop_bound = ksize_simint[jj+jsp];
          int o_loop_bound = ksize_simint[kk+ksp];
          //int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < L_LOOP_BOUND(isp); ++m)
          for(int n = 0; n < n_loop_bound; ++n)
          for(int o = 0; o < o_loop_bound; ++o)
          for(int p = 0; p < L_LOOP_BOUND(lsp); ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //printf("%s\n", ntotal-ncomputed+buffer_idx == eri_idx ? "true" : "false");
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_1010(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 1; ++isp) {
    for (int jsp = 0; jsp <= 0; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 1; ++ksp) {
        for (int lsp = 0; lsp <= 0; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          //int m_loop_bound = ksize_simint[ii+isp];
          int n_loop_bound = ksize_simint[jj+jsp];
          //int o_loop_bound = ksize_simint[kk+ksp];
          int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < L_LOOP_BOUND(isp); ++m)
          for(int n = 0; n < n_loop_bound; ++n)
          for(int o = 0; o < L_LOOP_BOUND(ksp); ++o)
          for(int p = 0; p < p_loop_bound; ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //printf("%s\n", ntotal-ncomputed+buffer_idx == eri_idx ? "true" : "false");
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_1011(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 1; ++isp) {
    for (int jsp = 0; jsp <= 0; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 1; ++ksp) {
        for (int lsp = 0; lsp <= 1; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          //int m_loop_bound = ksize_simint[ii+isp];
          int n_loop_bound = ksize_simint[jj+jsp];
          //int o_loop_bound = ksize_simint[kk+ksp];
          //int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < L_LOOP_BOUND(isp); ++m)
          for(int n = 0; n < n_loop_bound; ++n)
          for(int o = 0; o < L_LOOP_BOUND(ksp); ++o)
          for(int p = 0; p < L_LOOP_BOUND(lsp); ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //printf("%s\n", ntotal-ncomputed+buffer_idx == eri_idx ? "true" : "false");
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_1100(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 1; ++isp) {
    for (int jsp = 0; jsp <= 1; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 0; ++ksp) {
        for (int lsp = 0; lsp <= 0; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          //int m_loop_bound = ksize_simint[ii+isp];
          //int n_loop_bound = ksize_simint[jj+jsp];
          int o_loop_bound = ksize_simint[kk+ksp];
          int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < L_LOOP_BOUND(isp); ++m)
          for(int n = 0; n < L_LOOP_BOUND(jsp); ++n)
          for(int o = 0; o < o_loop_bound; ++o)
          for(int p = 0; p < p_loop_bound; ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //printf("%s\n", ntotal-ncomputed+buffer_idx == eri_idx ? "true" : "false");
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_1101(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 1; ++isp) {
    for (int jsp = 0; jsp <= 1; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 0; ++ksp) {
        for (int lsp = 0; lsp <= 1; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          //int m_loop_bound = ksize_simint[ii+isp];
          //int n_loop_bound = ksize_simint[jj+jsp];
          int o_loop_bound = ksize_simint[kk+ksp];
          //int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < L_LOOP_BOUND(isp); ++m)
          for(int n = 0; n < L_LOOP_BOUND(jsp); ++n)
          for(int o = 0; o < o_loop_bound; ++o)
          for(int p = 0; p < L_LOOP_BOUND(lsp); ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //printf("%s\n", ntotal-ncomputed+buffer_idx == eri_idx ? "true" : "false");
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_1110(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 1; ++isp) {
    for (int jsp = 0; jsp <= 1; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 1; ++ksp) {
        for (int lsp = 0; lsp <= 0; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          //int m_loop_bound = ksize_simint[ii+isp];
          //int n_loop_bound = ksize_simint[jj+jsp];
          //int o_loop_bound = ksize_simint[kk+ksp];
          int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < L_LOOP_BOUND(isp); ++m)
          for(int n = 0; n < L_LOOP_BOUND(jsp); ++n)
          for(int o = 0; o < L_LOOP_BOUND(ksp); ++o)
          for(int p = 0; p < p_loop_bound; ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //eri[mnop_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);

            //printf("%s\n", ntotal-ncomputed+buffer_idx == eri_idx ? "true" : "false");
            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}

void simgms_retrieve_eris_c_L_1111(int ii, int jj, int kk, int ll, double* eri, int* fullsizes) {

  int ncomputed = 0, ntotal = 0;

  //--start ERI computation--//
  for (int isp = 0; isp <= 1; ++isp) {
    for (int jsp = 0; jsp <= 1; ++jsp) {

      bool new_ij = ii+isp != iold || jj+jsp != jold;
      if (new_ij) {
        simint_create_multi_shellpair(1, &shells[ii+isp], 1, &shells[jj+jsp], &left_pair, 0);
        iold = ii+isp; jold = jj+jsp;
      }

      for (int ksp = 0; ksp <= 1; ++ksp) {
        for (int lsp = 0; lsp <= 1; ++lsp) {
          simint_create_multi_shellpair(1, &shells[kk+ksp], 1, &shells[ll+lsp], &right_pair, 0);
          ncomputed = simint_compute_eri(&left_pair, &right_pair, 0.0, work, buffer);

          int sizes[4] = { ksize_simint[ii+isp], ksize_simint[jj+jsp], ksize_simint[kk+ksp], ksize_simint[ll+lsp] };
          ncomputed *= sizes[0]*sizes[1]*sizes[2]*sizes[3];
          ntotal += ncomputed;

          //--sort separated L shells into proper JuliaChem L shell order--//

          int buffer_idx = 0;
          //#if 0
          //int m_loop_bound = ksize_simint[ii+isp];
          //int n_loop_bound = ksize_simint[jj+jsp];
          //int o_loop_bound = ksize_simint[kk+ksp];
          //int p_loop_bound = ksize_simint[ll+lsp];

          for(int m = 0; m < L_LOOP_BOUND(isp); ++m)
          for(int n = 0; n < L_LOOP_BOUND(jsp); ++n)
          for(int o = 0; o < L_LOOP_BOUND(ksp); ++o)
          for(int p = 0; p < L_LOOP_BOUND(lsp); ++p)
          {
            //int m_idx = kstart_simint[ii+isp] + m;
            //int n_idx = kstart_simint[jj+jsp] + n;
            //int o_idx = kstart_simint[kk+ksp] + o;
            //int p_idx = kstart_simint[ll+lsp] + p;

            eri[ntotal-ncomputed+buffer_idx] = buffer[buffer_idx];
            //printf("%d, %d, %d, %d, %lf\n", m_idx+1, n_idx+1, o_idx+1, p_idx+1, eri[eri_idx]);
            //printf("%s\n", ntotal-ncomputed+buffer_idx == eri_idx ? "true" : "false");
            ++buffer_idx;
          }
          //#endif
        }
      }
    }
  }
}
*/
#undef L_LOOP_BOUND
