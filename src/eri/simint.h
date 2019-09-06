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

/*
void simgms_retrieve_eris_c_L(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes, int* L_);

void simgms_retrieve_eris_c_0000(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_0001(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_0010(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_0011(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_0100(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_0101(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_0110(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_0111(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_1000(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_1001(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_1010(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_1011(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_1100(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_1101(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_1110(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void simgms_retrieve_eris_c_1111(int ii, int jj, int kk, int ll, double* ghondo, int* fullsizes);

void sort_L_shells_base(double* ghondo, double* target, int* sizes, int* L_,
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssss_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssp_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssd_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssf_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sssg_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssps_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspp_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspd_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspf_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sspg_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssds_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssds_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssdg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfs_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfs_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssff_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssfg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgs_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgs_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ssgg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spss_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsp_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsd_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsf_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spsg_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spps_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppp_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppd_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppf_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sppg_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spds_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spds_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spdg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfs_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfs_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spff_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spfg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgs_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgs_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_spgg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdsg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdpg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sddp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sddp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sddp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sddd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sddf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sddg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdfs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdfs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdfp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdfp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdfd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdfg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdgs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdgs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdgp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdgp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdgd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdgf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sdgg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfsg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfpg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfdp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfdp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfdd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfdf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfdg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sffs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sffs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sffs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sffp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sffp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sffp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sffd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sffg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfgs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfgs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfgp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfgp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfgd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfgf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sfgg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgsg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgpg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgdp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgdp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgdd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgdf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgdg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgfs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgfs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgfp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgfp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgfd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sgfg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sggs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sggs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sggs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sggp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sggp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sggp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sggd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sggf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_sggg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psss_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssp_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssd_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssf_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pssg_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psps_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspp_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspd_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspf_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pspg_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psds_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psds_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psdg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfs_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfs_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psff_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psfg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgs_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgs_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_psgg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppss_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsp_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsd_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsf_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppsg_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppps_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppp_1111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppd_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppf_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pppg_1110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppds_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppds_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppdg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfs_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfs_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppff_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppfg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgs_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgs_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgp_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgp_1101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgd_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgf_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ppgg_1100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdsg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdpg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pddp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pddp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pddp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pddd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pddf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pddg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdfs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdfs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdfp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdfp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdfd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdfg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdgs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdgs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdgp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdgp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdgd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdgf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pdgg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfsg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfpg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfdp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfdp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfdd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfdf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfdg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pffs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pffs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pffs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pffp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pffp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pffp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pffd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pffg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfgs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfgs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfgp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfgp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfgd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfgf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pfgg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgss_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgss_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgss_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgss_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgsg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgps_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgps_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgps_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgps_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpp_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpp_1011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpd_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpf_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgpg_1010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgds_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgds_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgdp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgdp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgdd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgdf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgdg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgfs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgfs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgfp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgfp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgfd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgff_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pgfg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pggs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pggs_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pggs_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pggp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pggp_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pggp_1001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pggd_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pggf_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_pggg_1000(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dssg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dspg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dsgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpsg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dppg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dpgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dddp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ddgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dffs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dffp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dfgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dgfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dggs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_dggp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fssg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fspg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fsgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpsg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fppg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fpgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fddp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fdgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fffs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fffp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ffgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fgfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fggs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_fggp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gssg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gspg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gsgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpss_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpss_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpss_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpss_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpsg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpps_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpps_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpps_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpps_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppp_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppp_0111(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppd_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppf_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gppg_0110(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpds_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpds_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpdp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpdp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpdd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpdf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpdg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpfs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpfs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpfp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpfp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpfd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpff_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpfg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpgs_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpgs_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpgp_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpgp_0101(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpgd_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpgf_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gpgg_0100(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gddp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gdgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gffs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gffp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfgs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gfgp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggss_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggss_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggss_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggsp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggsp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggsp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggsd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggsf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggsg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggps_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggps_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggps_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggpp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggpp_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggpp_0011(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggpd_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggpf_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggpg_0010(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggds_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggdp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggfs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_ggfp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gggs_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old); 

void sort_L_shells_gggp_0001(double* ghondo, double* target, 
        int* target_idx, int* buffer_idx, int* buffer_old);
*/ 
