////////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////////
#ifndef __TD1C__
#define __TD1C__
/////////////////////////////////////////////////////////////// external
//nyi #include <mpi.h>
#include <ctime>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
// Orimo_ECS
// prototype declaration
class Surff;
// Orimo_ECS
////////////////////////////////////////////////////////////////// local
typedef std::complex<double> dcomplex;
typedef std::vector<dcomplex* > type_cic;
typedef std::vector<std::vector<dcomplex* > > type_orb;
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// MPI
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clmpi
{
public:
  int nrank;                            // number of MPI processes
  int nthr;                             // number of OMP threads
  int myrank;                           // my rank
//nyi  int hostname_len;                      // hostname
//nyi  char hostname[MPI_MAX_PROCESSOR_NAME]; // hostname

  // basic
  clmpi();
  ~clmpi();
  void init();
  void final();
  void clear();
  void print() const;
  void gen();

  // functions
  void mpi_barrier() const;
  void omp_divide(int ithr, int ll0, int ul0, int& ll1, int& ul1) const;

private:
  void mpi_init();
  void mpi_final();
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Input and output
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clio
{
public:
  std::string name;
  std::string inp;
  std::string orb;
  std::string cic;
  std::string ene;
  std::string op1, opx, op0, opd, atrP, atrQ;
  std::string dipn, veln, accn;
  std::string ipx, ipd;
  std::string dipipx, velipx, accipx;
  std::string dipipd, velipd, accipd;
  std::string orbp;
  std::string rrad;
  std::string torb, tcic;
  std::string rhok, rhokz;
  FILE *fp_ene;
  FILE *fp_op1, *fp_opx, *fp_op0, *fp_opd, *fp_atrQ, *fp_atrP;
  FILE *fp_dipn, *fp_veln, *fp_accn;
  FILE *fp_ipx, *fp_ipd;
  FILE *fp_dipipx, *fp_velipx, *fp_accipx;
  FILE *fp_dipipd, *fp_velipd, *fp_accipd;
  FILE *fp_rrad;
  FILE *fp_torb, *fp_tcic;

  int iprint;
  std::string job_type;
  int nprint_ene;
  int nprint_op1, nprint_opx, nprint_opn, nprint_op1tr;
  int nprint_ipx, nprint_ipd, nprint_opipx, nprint_opipd;
  int nprint_norm;
  int nprint_rrad;
  int nprint_full;
  int nread_full;
// tdcis-teramura
  int nprint_every;
  std::string ci0;
  FILE *fp_ci0;
  int nprint_ci0;
// tdcis-teramura

  // basic
  clio();
  clio(std::string);
  clio(int, char**);
  ~clio();
  void print() const;
  void gen();
  void gen(std::string);
  void gen(int, char**);

  // functions
  void read_info(std::string, std::string&) const;
  void read_info(std::string, bool&) const;
  void read_info(std::string, int&) const;
  void read_info(std::string, double&) const;
  void read_info(std::string, dcomplex&) const;
  void read_info(std::string, std::vector<int>&) const;
  void read_info(std::string, std::vector<double>&) const;
  void read_info(std::string, std::vector<dcomplex>&) const;

  void read_info(std::string, std::string, std::string&) const;
  void read_info(std::string, bool, bool&) const;
  void read_info(std::string, int, int&) const;
  void read_info(std::string, double, double&) const;
  void read_info(std::string, dcomplex, dcomplex&) const;

private:
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Job control parameters
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clcontrol
{
public:
  static bool fedvr_normalized;
  static int icomp, igauge;
  static int oorot_type, split_type;
  static int iprojfc, type_fcx, type_dcx;
  static bool jfc_implicit, xfc_implicit;
  static int h1rat_maxcyc;
  static double h1rat_thresh;
  static int maxipx,maxipd;
  static double radipx;
  static bool docs1, docs2;
  static bool sae, psp;
  static int psp_type;
  static int dft_type;
  static int reg_type;
  static double throcc1;
  static double throcc2;
  static double throcc3;
  static int xact2_type;
  static int xact2_maxitr;
  static double xact2_thresh;
  static int ncut_occ3;
  static bool exact3j;
  static bool cionly;
// tdcis-teramura
  static bool istdcis, tdcis_rvg; //teramura tdcis 
  static bool print_xrad;         //teramura print xrad
// tdcis-teramura
  static int rrad_type;
  static double rrad_rion;
  static double rrad_r2in;
  static double rrad_r2out;
  static bool tsurff;

  // basic
  clcontrol();
  clcontrol(const clio&);
  ~clcontrol();
  void gen(const clio&);

private:
  static int num_control;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Laser field
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clfield1
{
public:
  int field_num;
  std::string env_type;
  double cyc1, cyc2, cyc3, numcyc, FWHM, sigma;
  double period, fint, famp, wlen, freq, dfreq, cep, tau, delay;
  double tleft, tcenter, tright;

  // basic
  clfield1();
  ~clfield1();
  void gen1(const clmpi&, const clio&, int field_num);
  void gen2(const clmpi&, const clio&, double tshift);

  // functions
  double ncyc() const;
  void get_value(double ttmp, double* lfield) const;
  void get_evalue(double ttmp, double* lfield) const;
  void get_avalue(double ttmp, double* lfield) const;
  void get_value(double t1, double t2, double c1, double c2, double* lfield) const;
  void get_der(double ttmp, double* lfield) const;

private:
  std::string gauge;
  void get_evalue_explicit(double ttmp, double* lfield) const;
  void get_evalue_implicit(double ttmp, double* lfield) const;
  void get_avalue_explicit(double ttmp, double* lfield) const;
  void get_avalue_implicit(double ttmp, double* lfield) const;
  void get_evalue_numder(double ttmp, double* lfield) const;
  void get_avalue_numint(double ttmp, double* lfield) const;
  void get_eder_explicit(double ttmp, double* lfield) const;
  void get_eder_implicit(double ttmp, double* lfield) const;
  void get_ader_explicit(double ttmp, double* lfield) const;
  void get_ader_implicit(double ttmp, double* lfield) const;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Laser field
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clfield
{
public:
  int num_field;
  int td_type;
  std::string gauge;
  bool lgauge, vgauge;
  int igauge;
  int nstep;
  double cyctot;
  double cycinit;
  int step;
  double time;
  double dtime;

  // basic
  clfield();
  ~clfield();
  clfield(const clmpi&, const clio&);
  void gen(const clmpi&, const clio&);

  // functions
  void init();
  double ncyc() const;
  bool finished() const;
  void get_value(double* lfield) const;
  void get_value(double ttmp, double* lfield) const;
  void get_value(double t1, double t2, double c1, double c2, double* lfield) const;
  void get_evalue(double* lfield) const;
  void get_evalue(double ttmp, double* lfield) const;
  void get_avalue(double* lfield) const;
  void get_avalue(double ttmp, double* lfield) const;
  void get_der(double* lfield) const;
  void get_der(double ttmp, double* lfield) const;
  std::vector<clfield1> Fields;
private:
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Physical observables
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
struct typhys
{
  double ene, ene1, ene2;
  double dip[3], vel[3], acc[3];
  double sz, s2;
  double lz, l2;
  double tot;
  double ipx[10];
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Discrete variable representation
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class cldvr
{
public:

  int ndvr;
  int mmin;
  int mmax;
  double xmin;
  double xmax;
  std::vector<double> xpt;
  std::vector<double> wpt;
  std::vector<double> dshape;
// Orimo_ECS
  double exp_alpha;
// Orimo_ECS

  // basic
  cldvr();
  cldvr(const int, const int, const int, const double, const double);
  ~cldvr();
  void init();
  void final();
  void clear();
  void print() const;
  void gen(const int, const int, const int, const double, const double);
  double get_val(int m, double rval) const;
// Orimo_ECS
  void gen_radau(int n, int m0, int m1, double x0, double x1, double d_exp_factor);
// Orimo_ECS

private:
  void gen_xpt();
  void gen_wpt();
  void gen_dshape();
  void get_legendre(double&, double&);
// Orimo_ECS
  void gen_xw_radau();
  void gen_dshape_radau();
// Orimo_ECS
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Finite element Discrete variable representation
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clfedvr
{
public:
  int nfe;                 // Number of finite elements
  int nrad;                // Number of radial grid points
  int nradfc;              // Number of radial grid points for frozen-core orbitals
  int nmax;                // Maximum numver of DVR functions
  double rmask;             // Edge of the boundary mask
  std::string mask_type;    // Type of mask function (cos/sawada)
  int mask_cos_order;      // Order of cos mask function

// Orimo_ECS
  int pot_type;     // Type of potential (atom, harmonic)
  int trunc_irad;
  std::string abc_type;     // Type of absorbing boundary condition(mask/ecs)
  double rab;               // r of Absorbing boundary
  int ecs_flag;            // Flag of ecs to give ecs info to fortran (ecs = 1 , other = 0)
  double theta;             // Angle/PI of the Exterior Complex Scaling, namely input coefficinet of PI
  double recs;              // Edge of the Exterior Complex Scaling
  int irad_ecs;            // boundary of ecs on irad
  bool inf_range;
  int irad_inf;
  double exp_factor;
  int type_mkint1_sph;
  int type_mkint2_sph;
  int type_mkv2mf;
  int type_mkxmat_aa;
  int switchoff;
  int irad_sw;
  int retreive_flag;       // Flag of retreiving wave funcion
// Orimo_ECS
// tdcis-teramura
  int nradgs;
// tdcis-teramura

  std::vector<int> mapf;   // Map from i to first \mu
  std::vector<int> mapb;   // Map from \mu to i
  std::vector<double> mask; // Mask function
  std::vector<double> xrad; // Global radial grid points
  std::vector<double> wrad; // Global radial quad weights

// Orimo_ECS
  std::vector<dcomplex> cxrad; //Complex global radial grid points
  std::vector<dcomplex> cwrad; // Complex global radial quad weights
  std::vector<double> bra_wrad; //for integralatiom in ecs
  std::vector<dcomplex> rdr; // R(r) / conjg(R(r))
  std::vector<dcomplex> wdw; // sqrt(cwrad / conjg(cwrad))
// Orimo_ECS

  //old  std::vector<double> kloc;      // Local kinetic energy matrix
  //old  std::vector<double> nabla_new; // Global first derivative
  std::vector<double> radk;             // radial kinetic energy
  std::vector<double> radk0;            // radial kinetic energy
  std::vector<double> radp;             // radial first derivative, transposed
  dcomplex radkI_ecs;                   // Rradial kinetic energy of ecs at R_0(ecs boundary)

  // basic
  clfedvr();
  clfedvr(const clmpi&, const clio&);
  ~clfedvr();
  void init();
  void final();
  void clear();
  void print() const;
  void gen(const clmpi&, const clio&);
  int get_ndvr(int ife) const;
  int get_ife(double rval) const;
  int get_irad(double rval) const;
  int get_irad(int ife, int m) const;
  double get_x0(int ife) const;
  double get_x1(int ife) const;
  double get_val(int irad, double rval) const;
  double get_val(int ife, int m, double rval) const;
  double get_val0(int ife, double rval) const;

  //  private: Changed by ORIMO to use these in class Basis
  double rll;               // Lower limit of r coordinate
  double rul;               // Upper limit of r coordinate
  std::vector<double> x0;   // Coordinates of left edge of elements {i}
  std::vector<double> x1;   // Coordinates of right edge of elements {i}
  std::vector<double> dfe;  // Size of elements {i}
  std::vector<int> ndvr;   // Number of DVR functions in elements {i}
  std::vector<cldvr> DVR;   // Local DVR, with (n+1) functions
  
private:
  void gen_x();
  void gen_m();
  void gen_map();
  void gen_grid();
  void gen_mask();
// Orimo_ECS
  void gen_abc(const clio& IO);
  void set_potential(const clio& IO);
// Orimo_ECS

  //old  void gen_kinetic();
  //old  void gen_nabla();
  //old  void gen_nabla_new();
  void gen_radk();
  void gen_radp();
  void read_x(const clio&);
};
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//// Spherical harmonics transformation
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
class clsph
{
public:
  int lmax1, mmax1, nsph1;
  int lmax2, mmax2, nsph2;
  int nlat, nphi, nang;
  double wgt_phi;
  std::vector<double> wgt_lat;
  std::vector<double> wgt_ang;
  std::vector<double> cost;
  std::vector<double> sint;
  std::vector<double> legf1, legf2;  // forward  Legendre transformation coefficients
  std::vector<double> legb1, legb2;  // backward Legendre transformation coefficients
// tdcis-teramura
  int lmax1gs, lmax2gs;
// tdcis-teramura

  clsph();
  clsph(const clmpi&, const clio&);
  ~clsph();
  void print(const clio&) const;
  void gen(const clmpi&, const clio&);

private:
  int read_mmax1(const clio&);
//shtns public:
//shtns   int sph_type;
//shtns   int lmax, mmax, nsph;
//shtns   int lmax2, mmax2;
//shtns   int nlat, nphi, nang;
//shtns   int nsph_shtns, nphi_shtns, nang_shtns;
//shtns   int const_mres;
//shtns   double wgt_phi;
//shtns   std::vector<double> cost;
//shtns   std::vector<double> sint;
//shtns   std::vector<double> wgt_lat;
//shtns   std::vector<double> wgt_ang;
//shtns   std::vector<double> legf;    // forward  Legendre transformation coefficients
//shtns   std::vector<double> legb;    // backward Legendre transformation coefficients
//shtns 
//shtns   // basic
//shtns   clsph();
//shtns   clsph(const clmpi&, const clio&, int type_sph);
//shtns   ~clsph();
//shtns   void init();
//shtns   void final();
//shtns   void clear();
//shtns   void print(const clio&) const;
//shtns   void gen(const clmpi&, const clio&, int type_sph);
//shtns 
//shtns   // extension
//shtns   void sph2ang(dcomplex*, dcomplex*);
//shtns   void ang2sph(dcomplex*, dcomplex*);
//shtns   void ang_sph_ang(const clmpi&, const int&);
//shtns   int read_mmax1(const clio&);
//shtns 
//shtns   // scratches
//shtns   //old  std::vector<dcomplex* > wfn_ang;
//shtns   //old  std::vector<dcomplex* > wfn_sph;
//shtns 
//shtns private:
};
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Determinants
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clormas
{
public:
  int nfcore2, nfcore1, nfcore, ndcore, ncore, nact, nocc, nvir, nfun;
  int nelcore[3], nelact[3], neltot[3];
  int nblock, nsub; 
  std::vector<int> mval, froz;
  std::vector<int> type_block, nfun_block;
  std::vector<int> norb_sub, lorb_sub, min_sub, max_sub;
  int nstr_alph, nstr_beta, ndet, lcic, ndetx;
  double thradfc;
  bool fab_den2, den2_abonly, donly;
  bool tdcc;
  int cc_type;
  std::string cc_code;

// tdcis-teramura
  int nfcore_tdcis;
// tdcis-teramura

  bool ras;
  bool dplus;
  int nact1,act1_ll,act1_ul;

  clormas();
  clormas(const clmpi&, const clio&);
  ~clormas();
  void gen(const clmpi&, const clio&);
  void cic0(dcomplex*) const;
  void check_orb(const clio&);
  void read_mval(const clio&);

private:
};
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Determinants
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class cldet
{
public:
  int nfcore, ndcore, ncore, nact, nvir, nocc, nfun;
  int neltot[3], nelact[3];
  int nstra, nstrb;

  std::vector<int> wgta;
  std::vector<int> onva;
  std::vector<int> orba;
  std::vector<int> n1xa;
  std::vector<int> p1xa;
  std::vector<int> h1xa;
  std::vector<int> eq1xa;
  std::vector<int> sgn1xa;
  std::vector<int> n1xra;
  std::vector<int> r1xra;
  std::vector<int> l1xra;
  std::vector<int> sgn1xra;

  std::vector<int> wgtb;
  std::vector<int> onvb;
  std::vector<int> orbb;
  std::vector<int> n1xb;
  std::vector<int> p1xb;
  std::vector<int> h1xb;
  std::vector<int> eq1xb;
  std::vector<int> sgn1xb;
  std::vector<int> n1xrb;
  std::vector<int> r1xrb;
  std::vector<int> l1xrb;
  std::vector<int> sgn1xrb;

  // basic
  cldet();
  cldet(const clmpi&, const clio&);
  ~cldet();
  void init();
  void final();
  void clear();
  void alloc();
  void print() const;
  void read_info(const clio&);
  void read_info_orb(const clio&);
  void gen(const clmpi&, const clio&);

  // functions
private:
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Basis
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clbas
{
public:
  int znuc;     // nuclear charge
  int smul;     // 2S*1 of S(S+1)
  int ltot;     // L of L(L+1)
  int mtot;     // L_z
  int nbas;     // number of 1e basis functions for orbitals
  int nbas2;    // number of 1e basis functions for Coulomb potential
  int ngrid;    // number of spatial grid points
  bool lconst;   // fix l of each orbital
  bool lconst_core;
  int psp_label; // pseudopotential label

  clfedvr GRad;  // radial grid information
  clsph GAng;    // angular grid information
  clsph Sph1;    // angular grid information
  clsph Sph2;    // angular grid information
  //old  cldet Det;     // Determinant basis informatin
  clormas ORMAS; // Determinant-based CI information

  std::vector<int> nval, lval, mval;
  std::vector<double> grid, wgt;
  std::vector<double> alph_lm, alph_rlm;
  std::vector<dcomplex> tmat, kmat, pmat; 
  std::vector<double> d2ll;
  std::vector<dcomplex> bas_zfac, bas_pzfac1, bas_pzfac2, bas_azfac;
  std::vector<double> bas_d2fac1, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1;
// Orimo_ECS
  std::vector<dcomplex> bas_d2crpl1;
  std::vector<dcomplex> d1mat, d2mat; // ecs contour redial first, second derivative
  std::vector<double> confd2ll;
  std::vector<dcomplex> d2ll_ecs; // for ecs
  std::vector<int> ipiv_ecs;         // for ecs
// Orimo_ECS

  // basic
  clbas();
  clbas(const clmpi&, const clio&);
  ~clbas();
  void init();
  void final();
  void clear();
  void gen(const clmpi&, const clio&);
  void read_info(const clio&);
  void print(const clio&) const;
  void ang2sph1(const clmpi&, const std::vector<dcomplex>& fang, std::vector<dcomplex>& fsph) const;
  void sph2ang1(const clmpi&, const std::vector<dcomplex>& fsph, std::vector<dcomplex>& fang) const;
  void ang2sph2(const clmpi&, const std::vector<dcomplex>& fang, std::vector<dcomplex>& fsph) const;
  void sph2ang2(const clmpi&, const std::vector<dcomplex>& fsph, std::vector<dcomplex>& fang) const;
  void proj(const clbas&, const std::vector<dcomplex>& fun1, std::vector<dcomplex>& fun2) const;
  void ppgenkb(const std::vector<dcomplex>& orb) const;
private:
  static int num_bas;
  void gen_grid();
  void gen_wgt();
  void read_lmval(const clio&);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Wavefunction
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clwfn
{
public:
  int size;  // total size of wavefunction
  int size1; // orbital-part size (angular momentum-space)
  int sizeg; // orbital-part size (real-space)
  int size2; // ci-part sizes
  std::vector<int> typefc; // angular momentum of frozen-core orbitals

  std::vector<dcomplex > wfn;
  std::vector<dcomplex > wfng;
//old  type_cic cic;  // pointer to ci part: cic[nstrb][nstra]
//old  type_orb orb;  // pointer to orbital part: orb[nfun][nlm][nrad+1]
//old  type_orb orbg; // temporary space: orbg[nfun][nang][nrad+1]

// Sato_tSURFF
// Orimo_ECS
  // for t-SURFF
  bool dosurff;
  Surff *surff;
// Orimo_ECS
// Sato_tSURFF

  // basic
  clwfn();
  clwfn(bool);
  clwfn(const clbas&, const clwfn&);
  clwfn(const clmpi&, const clio&, const clbas&);
  clwfn(const clmpi&, const clio&, const clbas&, bool);
  ~clwfn();
  void init();
  void final();
  void clear(const clbas&);
  void clearo(const clbas&);
  void clearc(const clbas&);
  void alloc(const clbas&);
  void copy(const clbas&, const clwfn&);
  void adjust(const clbas&, const clwfn&);
  void gen(const clmpi&, const clio&, const clbas&);
  void get_nradfc(const clmpi&, const clio&, const clbas&);

  void print() const;
  void print(const clbas&) const;
  void printc(const clbas&) const;
  void printo(const clbas&) const;
  void read(const clmpi&, const clio&, const clbas&);
  void read_orb(const clmpi&, const clio&, const clbas&);
  void read_cic(const clmpi&, const clio&, const clbas&);
  void write(const clmpi&, const clio&, const clbas&) const;
  void read_info(const clio&);

  // functions
  void orth(const clbas&);
  void ortho(const clbas&);
  void orthc(const clbas&);
  void proj(const clbas&, const clwfn&);
  void projg(const clmpi&, const clbas&, const clwfn&);
  void mask(const clbas&);
  void ladapt(const clmpi& Proc, const clbas& Bas);
  void madapt(const clmpi& Proc, const clbas& Bas);
  void print_cic(const clbas&) const;
  void print_cic(std::string, const clbas&) const;
  void print_orb(const clbas&) const;
  void print_orb(int, const clio&, const clbas&) const;
  void print_orb(int, std::string, const clbas&) const;
  void print_orbg(const clbas&) const;
  void print_wang(const clbas&) const;
// tdcis-teramura
  void read_nradgs(const clmpi& Proc, const clio& IO, const clbas& Bas);
// tdcis-teramura

private:
  std::string orth_type;    // Type of mask function (cos/sawada)
  bool ci_normalize;        // normalize ci vector in real time propagation
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Hamiltonian times wavefunction
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clhprod
{
public:
  clhprod();
  ~clhprod();
  clhprod(const clmpi&, const clio&, const clbas&, const clfield&);
  void gen(const clmpi&, const clio&, const clbas&, const clfield&);

  bool projhigh;
  double projhigh_cutoff;
  double dip_exp, vel_exp, acc_exp;
  double ene_fcore, ene_dcore, ene_core, ene_act, ene_tot;

  void clear(const clmpi&, const clbas&, clwfn& Wfn);
  void clearo(const clmpi&, const clbas&, clwfn& Wfn);
  void clearc(const clmpi&, const clbas&, clwfn& Wfn);
  void clearog(const clmpi&, const clbas&, clwfn& Wfn);
  void copy(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void copyo(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void copyc(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void copyog(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void scal(const clmpi&, const clbas&, dcomplex fac, clwfn& Wfn);
  void scalo(const clmpi&, const clbas&, dcomplex fac, clwfn& Wfn);
  void scalc(const clmpi&, const clbas&, dcomplex fac, clwfn& Wfn);
  void scalog(const clmpi&, const clbas&, dcomplex fac, clwfn& Wfn);
  void xpy(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void xpyo(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void xpyc(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void xpyog(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void xpyz(const clmpi&, const clbas&, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void xpyzo(const clmpi&, const clbas&, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void xpyzc(const clmpi&, const clbas&, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void xpyzog(const clmpi&, const clbas&, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void xmy(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void xmyo(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void xmyc(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void xmyog(const clmpi&, const clbas&, const clwfn& Wfn1, clwfn& Wfn2);
  void xmyz(const clmpi&, const clbas&, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void xmyzo(const clmpi&, const clbas&, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void xmyzc(const clmpi&, const clbas&, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void xmyzog(const clmpi&, const clbas&, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void axpy(const clmpi&, const clbas&, dcomplex fac, const clwfn& Wfn1, clwfn& Wfn2);
  void axpyo(const clmpi&, const clbas&, dcomplex fac, const clwfn& Wfn1, clwfn& Wfn2);
  void axpyc(const clmpi&, const clbas&, dcomplex fac, const clwfn& Wfn1, clwfn& Wfn2);
  void axpyog(const clmpi&, const clbas&, dcomplex fac, const clwfn& Wfn1, clwfn& Wfn2);
  void axpyz(const clmpi&, const clbas&, dcomplex fac, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void axpyzo(const clmpi&, const clbas&, dcomplex fac, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void axpyzc(const clmpi&, const clbas&, dcomplex fac, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void axpyzog(const clmpi&, const clbas&, dcomplex fac, const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3);
  void axpbyz(const clmpi&, const clbas&, dcomplex fac1, const clwfn& Wfn1, dcomplex fac2, const clwfn& Wfn2, clwfn& Wfn3);
  void axpbyzo(const clmpi&, const clbas&, dcomplex fac1, const clwfn& Wfn1, dcomplex fac2, const clwfn& Wfn2, clwfn& Wfn3);
  void axpbyzc(const clmpi&, const clbas&, dcomplex fac1, const clwfn& Wfn1, dcomplex fac2, const clwfn& Wfn2, clwfn& Wfn3);
  void axpbyzog(const clmpi&, const clbas&, dcomplex fac1, const clwfn& Wfn1, dcomplex fac2, const clwfn& Wfn2, clwfn& Wfn3);

  void mkrrad(const clwfn& Wfn, std::vector<dcomplex>& rrad, std::vector<dcomplex>& rradpw);
  void mkrradx(const clwfn& Wfn, std::vector<dcomplex>& rrad, std::vector<dcomplex>& rradpw);
  void mkrrad0(const clwfn& Wfn, std::vector<dcomplex>& rrad, std::vector<dcomplex>& rradpw);
  void mkrrad1(const clwfn& Wfn, std::vector<dcomplex>& rrad, std::vector<dcomplex>& rradpw);
  void spin(const clmpi&, const clbas&, const clwfn& Wfn, typhys&);
  void oang(const clmpi&, const clbas&, const clwfn& Wfn, typhys&);
  void ladapt(const clmpi&, const clbas&, clwfn& Wfn);
  void dipole(const clmpi&, const clbas&, const double* lfield, const clwfn& Wfn, typhys&);
  void dipole_ipx(const clmpi&, const clbas&, const double* lfield, const clwfn& Wfn, typhys*);
  void dipole_ipd(const clmpi&, const clbas&, const double* lfield, const clwfn& Wfn, typhys*);
  void energy(const clmpi&, const clbas&, const double* lfield, const clwfn& Wfn, typhys&);
  void chkconv(const double* lfield, const clwfn& Wfn, typhys&);
  void cidiag(const clmpi&, const clbas&, clwfn& Wfn);
  void orbene(const clmpi&, const clbas&, const double* lfield, const clwfn& Wfn, std::vector<double>& Eig);
  void fockdiag(const clmpi&, const clbas&, const double* lfield,       clwfn& Wfn, std::vector<double>& Eig);
  void enepole(const clmpi&, const clbas&, const double* lfield, const clwfn& Wfn, typhys&);
  void normx(const clmpi&, const clbas&, const double* lfield, const clwfn&, typhys&);
  void ionpx(const clmpi&, const clbas&, const double* lfield, const clwfn&, typhys&);
  void ionpd(const clmpi&, const clbas&, const double* lfield, const clwfn&, typhys&);

  void set_wfn0(const clmpi&, const clio&, const clbas&, const clwfn& Wfn) const;
  void ProjHigh(clwfn& Wfn) const;
  void htot(const clmpi&, const clbas&, double dt, const double* lfield, const clwfn& Wfn, clwfn& hWfn);
  void htotx(const clmpi&, const clbas&, double dt, const double* lfield, const clwfn& Wfn);
  void htoto(const clmpi&, const clbas&, double dt, const double* lfield, const clwfn& Wfn, clwfn& hWfn);
  void htot0(const clmpi&, const clbas&, double dt, const double* lfield, const clwfn& Wfn, clwfn& hWfn);
  void h1tot(const clmpi&, const clbas&, double dt, const double* lfield, const clwfn& Wfn, clwfn& hWfn);
  void v1tot(const clmpi&, const clbas&, double dt, const double* lfield, const clwfn& Wfn, clwfn& hWfn);
  void h1add(const clmpi&, const clbas&, dcomplex zfac, const double* lfield, const clwfn& Wfn, clwfn& hWfn);
  void v1ext(const clmpi&, const clbas&, dcomplex zfac, const double* lfield, const clwfn& Wfn, clwfn& hWfn);
  void fulldiag(const clmpi&, const clbas&, clwfn& wfn);
  void getno(const clmpi&, const clio&, const clbas&, const clwfn&, clwfn&);
  void kickk(const clmpi&, const clio&, const clbas&, clwfn&, double knorm) const;
  void getden1(std::vector<dcomplex>&) const;
  void getden2(std::vector<dcomplex>&) const;
  void getint1(std::vector<dcomplex>&) const;
  void getint2(std::vector<dcomplex>&) const;

// tdcis-teramura
  void phys_tdcis(const clmpi&, const clbas&, const double* lfield, const clwfn& Wfn, typhys&, const double);
  void orbin_tdcis_init();
// tdcis-teramura

  void op1tr(const double* lfield, typhys&) const;
  void op1tr_printp_tag(const clio&, const clbas&) const;
  void op1tr_printq_tag(const clio&, const clbas&) const;
  void op1tr_print(const clio&, const clbas&, const clfield&, const double*, const typhys&) const;
  void op1tr_init(const clio&, const clbas&, const clwfn&);
  void op1tr_final();
  int op1tr_nfun,op1tr_nrad;
  std::vector<dcomplex> vdP,vvP,vaP;
  std::vector<dcomplex*> dP,vP,aP;
  std::vector<dcomplex> dQ,vQ,aQ;
private:
  static int num_hprod;
//old  void hcic(const clmpi&, const clbas&);
//old  void mkden(const clmpi&, const clbas&);
//old  void mkden1(const clmpi&, const clbas&);
//old  void mkden2(const clmpi&, const clbas&);
//old  void invden(const clmpi&, const clbas&);
//old  void d2rden(const clmpi&, const clbas&);
//old  void mkint(const clmpi&, const clbas&);
//old  void mkint1(const clmpi&, const clbas&);
//old  void mkint2(const clmpi&, const clbas&);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Gaussian basis function
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clgbas
{
public:

  std::string fck;
  int ngbas;
  int ngfun;
  int glmax;
  int nshell;
  int nshelp;
  std::vector<int> type;
  std::vector<int> nprm;
  std::vector<double> alph;
  std::vector<double> cont;
  std::vector<double> cmo;

  clgbas(const clio&);
  ~clgbas();
  void read_fck_info(const std::string, std::string&) const;
  void read_fck_info(const std::string, bool&) const;
  void read_fck_info(const std::string, int&) const;
  void read_fck_info(const std::string, double&) const;
  void read_fck_array(const std::string, int, std::vector<int>&) const;
  void read_fck_array(const std::string, int, std::vector<double>&) const;

private:
  void overlap();
  void normalize1();
  void normalize2();
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Rational functions of one-electron operators
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clh1rat
{
public:
  clh1rat();
  ~clh1rat();
  void gen(const clmpi&, const clio&, const clbas&, double, int, int, int,
	   std::vector<dcomplex>&, std::vector<dcomplex>&, std::vector<dcomplex>&);
  void prod(const clmpi&, const clbas&, double time, const clfield&, clhprod&, clwfn&) const;
  void prod(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, clwfn&) const;
  void prod2(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, clwfn&, clwfn&) const;   
  void prod_numer(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, const clwfn&, clwfn&) const;
  void prod_denom(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, clwfn&, clwfn&) const;
  void prod2(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, 
	     const std::vector<dcomplex>&, const std::vector<dcomplex>&, const std::vector<dcomplex>&, clwfn&, clwfn&) const;
  void prod_numer(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, 
		  const std::vector<dcomplex>&, const std::vector<dcomplex>&, const std::vector<dcomplex>&, const clwfn&, clwfn&) const;
  void prod_denom(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, 
		  const std::vector<dcomplex>&, const std::vector<dcomplex>&, const std::vector<dcomplex>&, clwfn&, clwfn&) const;
  void prod2(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, double, 
	     const std::vector<dcomplex>&, const std::vector<dcomplex>&, const std::vector<dcomplex>&, clwfn&, clwfn&) const;
  void prod_numer(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, double, 
		  const std::vector<dcomplex>&, const std::vector<dcomplex>&, const std::vector<dcomplex>&, const clwfn&, clwfn&) const;
  void prod_denom(const clmpi&, const clbas&, double time, dcomplex zfac, const clfield&, clhprod&, double, 
		  const std::vector<dcomplex>&, const std::vector<dcomplex>&, const std::vector<dcomplex>&, clwfn&, clwfn&) const;
  dcomplex get_numer_limit() const;
  dcomplex get_denom_limit() const;
  dcomplex get_coeff_limit() const;
  void print() const;

private:
  double dtime;
  int prodci_type;
  int h1ci_type;
  int dim_numer;
  int dim_denom;
  dcomplex numer_limit;
  dcomplex denom_limit;
  dcomplex coeff_limit;
  std::vector<dcomplex> ncoeff;
  std::vector<dcomplex> dcoeff0;
  std::vector<dcomplex> dcoeff1;
  std::vector<int>     tpiv;
  std::vector<dcomplex> tinv;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Generalized denominator of Crank-Nicolson operator
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class cldpade
{
public:
  cldpade();
  ~cldpade();
  cldpade(const clmpi&, const clio&, const clbas&, double dt, dcomplex alpha, int icomp, int isplit);
  void gen(const clmpi&, const clio&, const clbas&, double dt, dcomplex alpha, int icomp, int isplit);
  void gen2(const clmpi&, const clio&, const clbas&, dcomplex d0, dcomplex d1);
  void prod(const clmpi&, const clbas&, double time, const clfield&, clhprod&, clwfn&) const;
  void prod2(const clmpi&, const clbas&, double time, const clfield&, clhprod&, clwfn&) const;
  void cnic(const clmpi&, const clbas&, double time, const clfield&, clhprod&, clwfn&) const;

  bool dpade_midpt;

private:
  int icomp;
  int dpade_maxcyc;
  double dpade_thresh;
  double dpade_dtime;
  dcomplex dpade_alpha;
  dcomplex dpade_coeff0, dpade_coeff1;

  std::vector<int> tpiv;
  std::vector<dcomplex> tinv;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Basic class for h1 propagation subclasses
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clh1prop
{
public:
  clh1prop();
  virtual ~clh1prop();
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&) = 0;
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, double dtime, const clhprod&) = 0;
  virtual void prop(const clmpi&, const clbas&, const clfield&, int STEP, clhprod&, clwfn&) = 0;
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&) = 0;

protected:
  bool projfc;
  void gen_basic(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// alternate direction implicit propagator
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clh1adi : public clh1prop
{
public:
  clh1adi();
  virtual ~clh1adi();
  clh1adi(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, double dtime, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, int STEP, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  bool doext;
  std::vector<int> tpiv2;
  std::vector<dcomplex> tadi1, tadi2;
  void propl(const clmpi&, const clbas&, const clfield&, int STEP, clhprod&, clwfn&);
  void propv(const clmpi&, const clbas&, const clfield&, int STEP, clhprod&, clwfn&);
  void t_explicit(const clmpi&, const clbas&, clhprod&, clwfn&);
  void t_implicit(const clmpi&, const clbas&, clhprod&, clwfn&);
  void laser_lgauge(const clmpi&, const clbas&, int istag, double time, double dt, const clfield&, clwfn&);
  void laser_vgauge(const clmpi&, const clbas&, int istag, double time, double dt, const clfield&, clwfn&);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Matrix iteration method
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clh1itr : public clh1prop
{
public:
  clh1itr();
  virtual ~clh1itr();
  clh1itr(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  clh1itr(const clmpi&, const clio&, const clbas&, const clfield&, double dtime, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, double dtime, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, int STEP, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);

private:
  int h1itr_maxcyc;
  double h1itr_thresh;
  std::vector<int> cnpiv;
  std::vector<dcomplex> cninv;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Lancosz exponentialization of one-electron hamiltonian
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clh1exp : public clh1prop
{
public:
  clh1exp();
  virtual ~clh1exp();
  clh1exp(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, double dtime, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, int STEP, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);

private:
  int h1exp_maxcyc;
  double h1exp_thresh;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Basic class for (h1 +) h2 propagation subclasses
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clh2prop
{
public:
  clh2prop();
  virtual ~clh2prop();
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&) = 0;
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&) = 0;
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&) = 0;

protected:
  void gen_basic(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Constant time-step Runge-Kuta propagator
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clcrk : public clh2prop
{
public:
  clcrk();
  virtual ~clcrk();
  clcrk(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  int rk_order;
  int rk_level;
  int rk_nstage;
  int rk_max_stage;
  std::string rk_type;
  std::string rk_formula;
  std::vector<double> rk_nodes;
  std::vector<std::vector<double> > rk_wghts;
  std::vector<std::vector<double> > rk_coeff;

  clwfn Wfn0;
  clwfn tWfn;
  std::vector<clwfn> kWfn;

  void set_coeff();
  void set_coeff_01();
  void set_coeff_02();
  void set_coeff_03();
  void set_coeff_04();
  void set_coeff_04_kutta();
  void set_coeff_54_fehlberg();
  void set_coeff_54_dormand();
  void set_coeff_65_prince();
  void set_coeff_65_calvo();
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Constant time-step implicit Runge-Kuta propagator
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clirk : public clh2prop
{
public:
  clirk();
  virtual ~clirk();
  clirk(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  int irk_order;
  int irk_maxcyc;
  double irk_thresh;
  std::string irk_formula;
  clwfn Wfn0, Wfn1, hWfn0, hWfn1;
  void prop1(const clmpi&, const clbas&, double, double, const clfield&, clhprod&, clwfn&);
  void prop2m(const clmpi&, const clbas&, double, double, const clfield&, clhprod&, clwfn&);
  void prop2t(const clmpi&, const clbas&, double, double, const clfield&, clhprod&, clwfn&);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Exponential Runge-Kuta with Pade approximation for phi functions
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class cletdrb : public clh2prop
{
public:
  cletdrb();
  virtual ~cletdrb();
  cletdrb(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  double etd_dt1, etd_dt2;
  clwfn  Wfn0,  Wfn1,  Wfn2,  Wfn3;
  clwfn hWfn0, hWfn1, hWfn2, hWfn3, tWfn, dWfn;
  double Eref;
  std::vector<dcomplex> Den1, Int1e, Int2e, hDiag;

  int time_of_h;
  bool pade_equiv;
  bool rosenbrock;
  static int max_dim_numer;
  static int max_dim_denom;
  int rk_order;
  int dim_numer1, dim_denom1;
  int dim_numer2, dim_denom2;
  int dim_numer3, dim_denom3;
  int cisplit;

  dcomplex phi1rk, phid1rk, phi2rk, phi3rk;
  clh1rat phi1dt1, phid1dt1, phi2dt1, phi3dt1;
  clh1rat phi1dt2, phid1dt2, phi2dt2, phi3dt2;

  void gen_phi1(const clmpi&, const clio&, const clbas&);
  void gen_phid1(const clmpi&, const clio&, const clbas&);
  void gen_phi2(const clmpi&, const clio&, const clbas&);
  void gen_phi3(const clmpi&, const clio&, const clbas&);
  void prop1(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop2(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop3(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4_orb(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4_orbci(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4v1(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4v2(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4v3(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4v4(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prod1(const clmpi&, const clbas&, clhprod&, double, const clfield&, const clh1rat&, dcomplex, const clwfn& WIn, clwfn& WOut);
  void prod2(const clmpi&, const clbas&, clhprod&, double, const clfield&, const clh1rat&, dcomplex, clwfn& WIn, clwfn& WOut);
  dcomplex test_pade(const std::vector<dcomplex>&, const std::vector<dcomplex>&, const std::vector<dcomplex>&, dcomplex) const;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// ETD Adams-Bashforth method with Pade approximation for phi functions
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clexpab : public clh2prop
{
public:
  clexpab();
  virtual ~clexpab();
  clexpab(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:

  double expab_dt;
  dcomplex expab_cdt;
  std::vector<double> timeP;
  clwfn tWfn;
  std::vector<clwfn> WfnB, WfnP, hWfnP;
  static int max_dim_numer;
  static int max_dim_denom;
  int num_step;
  bool expab_pc;
  int expab_order;
  int expab_npade;
  int expab_dpade;
  int expab_ndpade;
  int expab_h1type;
  clh1rat expab_phi;

  //for debug
  clh1rat expab_phi0;
  clh1rat expab_phi1;
  void prop1(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop2(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop2_0(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop2_1(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void gen_phi2_phi0(const clmpi&, const clio&, const clbas&);
  void gen_phi2_phi1(const clmpi&, const clio&, const clbas&);
  //for debug
  void prop_general(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop_general_corrector(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void gen_phi1(const clmpi&, const clio&, const clbas&);
  void gen_phi2(const clmpi&, const clio&, const clbas&);
  void gen_phi3(const clmpi&, const clio&, const clbas&);
  void gen_phi4(const clmpi&, const clio&, const clbas&);
  void gen_phi5(const clmpi&, const clio&, const clbas&);
  void gen_phi6(const clmpi&, const clio&, const clbas&);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Lawson Adams-Bashforth method with Pade approximation for phi functions
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class cllawab : public clh2prop
{
public:
  cllawab();
  virtual ~cllawab();
  cllawab(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:

  double lawab_dt;
  dcomplex lawab_cdt;
  std::vector<double> timeP;
  clwfn tWfn;
  std::vector<clwfn> WfnP, hWfnP;
  static int max_dim_numer;
  static int max_dim_denom;
  static int max_order;
  static dcomplex coeffAB[10][10];
  static dcomplex coeffAM[10][10];

  int num_step;
  int lawab_order;
  int lawab_npade;
  int lawab_dpade;
  int lawab_ndpade;
  int lawab_h1type;
  int lawab_ncorr;
  std::vector<clh1rat> lawab_exp;
  void prop1(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop2(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop_general(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop_general(const clmpi&, const clbas&, int, const clfield&, clhprod&, clwfn&);
  void gen_exp(const clmpi&, const clio&, const clbas&);
  void get_coeff(dcomplex*, dcomplex*) const;
  void get_coeff(int, dcomplex*, dcomplex*);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Exponential Runge-Kuta with Pade approximation for phi functions
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class cletdrk : public clh2prop
{
public:
  cletdrk();
  virtual ~cletdrk();
  cletdrk(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  double etd_dt1, etd_dt2;
  clwfn  Wfn0,  Wfn1,  Wfn2,  Wfn3;
  clwfn hWfn0, hWfn1, hWfn2, hWfn3, tWfn;

  bool pade_equiv;
  bool rosenbrock;
  static int max_dim_numer;
  static int max_dim_denom;
  int rk_order, dim_numer, dim_denom;

  dcomplex phi1rk, phi2rk, phi3rk;
  clh1rat phi1dt1, phi2dt1, phi3dt1;
  clh1rat phi1dt2, phi2dt2, phi3dt2;

  void gen_phi1(const clmpi&, const clio&, const clbas&);
  void gen_phi2(const clmpi&, const clio&, const clbas&);
  void gen_phi3(const clmpi&, const clio&, const clbas&);
  void prop1(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop2(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop3(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prod1(const clmpi&, const clbas&, clhprod&, double, const clfield&, const clh1rat&, dcomplex, const clwfn& WIn, clwfn& WOut);
  void prod2(const clmpi&, const clbas&, clhprod&, double, const clfield&, const clh1rat&, dcomplex, clwfn& WIn, clwfn& WOut);
  dcomplex test_pade(const std::vector<dcomplex>&, const std::vector<dcomplex>&, const std::vector<dcomplex>&, dcomplex) const;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Exponential time differencing second-order Runge-Kuta with Pade(1,1)
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class cletdrk2 : public clh2prop
{
public:
  cletdrk2();
  virtual ~cletdrk2();
  cletdrk2(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  cldpade dPade1;
  clwfn Wfn0, Wfn1, hWfn0, hWfn1, tWfn;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Exponential time differencing 4th-order Runge-Kuta with Pade(2,2)
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class cletdrk4 : public clh2prop
{
public:
  cletdrk4();
  virtual ~cletdrk4();
  cletdrk4(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  double etd_dt1, etd_dt2;
  clwfn  Wfn0,  Wfn1,  Wfn2,  Wfn3;
  clwfn hWfn0, hWfn1, hWfn2, hWfn3, tWfn;

  static int max_dim_numer;
  static int max_dim_denom;
  int dim_numer, dim_denom;
  bool pade_old;
//  std::vector<dcomplex> exp_n, gfun_n, f0fun_n, f1fun_n, f2fun_n;
//  std::vector<dcomplex> exp_d0, gfun_d0, f0fun_d0, f1fun_d0, f2fun_d0;
//  std::vector<dcomplex> exp_d1, gfun_d1, f0fun_d1, f1fun_d1, f2fun_d1;

  clh1rat exp1;
  clh1rat exp2;
  clh1rat gfun;
  clh1rat f0fun;
  clh1rat f1fun;
  clh1rat f2fun;

  dcomplex exp1rk;
  dcomplex exp2rk;
  dcomplex gfunrk;
  dcomplex f0funrk;
  dcomplex f1funrk;
  dcomplex f2funrk;

  void gen_exp(const clmpi&, const clio&, const clbas&);
  void gen_gfun(const clmpi&, const clio&, const clbas&);
  void gen_f0fun(const clmpi&, const clio&, const clbas&);
  void gen_f1fun(const clmpi&, const clio&, const clbas&);
  void gen_f2fun(const clmpi&, const clio&, const clbas&);
  void prod1(const clmpi&, const clbas&, clhprod&, double, const clfield&, const clh1rat&, dcomplex, const clwfn& WIn, clwfn& WOut);
  void prod2(const clmpi&, const clbas&, clhprod&, double, const clfield&, const clh1rat&, dcomplex, clwfn& WIn, clwfn& WOut);
  dcomplex test_pade(const std::vector<dcomplex>&, const std::vector<dcomplex>&, const std::vector<dcomplex>&, dcomplex) const;
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Exponential OIFS propagator with 4th-order Runge-Kuta time stepping
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class cllawrk : public clh2prop
{
public:
  cllawrk();
  virtual ~cllawrk();
  cllawrk(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  int lawrk_order, lawrk_npade;
  double lawrk_dt1, lawrk_dt2;
  cldpade dPade2p, dPade2m;
  clwfn Wfn0, Wfn1, hWfn0, hWfn1, v0Wfn, vWfn, xWfn;
  void prop1(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop2(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop3(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4_old(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void expwfn(const clmpi&, const clbas&, const clfield&, clhprod&, const clwfn&, clwfn&);
  void expwfn11(const clmpi&, const clbas&, const clfield&, clhprod&, const clwfn&, clwfn&);
  void expwfn12(const clmpi&, const clbas&, const clfield&, clhprod&, const clwfn&, clwfn&);
  void expwfn22(const clmpi&, const clbas&, const clfield&, clhprod&, const clwfn&, clwfn&);
};
////////////////////////////////////////////////////////////////////////
class cllawmul : public clh2prop
{
public:
  cllawmul();
  virtual ~cllawmul();
  cllawmul(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  int lawrk_order, lawrk_npade;
  double lawrk_dt1, lawrk_dt2;
  cldpade dPade2p, dPade2m;
  clwfn Wfn0, Wfn1, hWfn0, hWfn1, v0Wfn, vWfn, xWfn;
  void prop1(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop2(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop3(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void prop4_old(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  void expwfn(const clmpi&, const clbas&, const clfield&, clhprod&, const clwfn&, clwfn&);
  void expwfn11(const clmpi&, const clbas&, const clfield&, clhprod&, const clwfn&, clwfn&);
  void expwfn12(const clmpi&, const clbas&, const clfield&, clhprod&, const clwfn&, clwfn&);
  void expwfn22(const clmpi&, const clbas&, const clfield&, clhprod&, const clwfn&, clwfn&);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Newton iterated Crank Nicolson with approximate Jacobian
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clcrnic : public clh2prop
{
public:
  clcrnic();
  virtual ~clcrnic();
  clcrnic(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void gen(const clmpi&, const clio&, const clbas&, const clfield&, const clhprod&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, clhprod&, clwfn&);
  virtual void prop(const clmpi&, const clbas&, const clfield&, double time, double dtime, clhprod&, clwfn&);
private:
  int crnic_maxcyc;
  cldpade dPade;
  clwfn Wfn0, Wfn1, hWfn0, hWfn1;
  double get_res(const clmpi&, const clbas&, const clwfn&);
};
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Determinants
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class clpes
{
public:
  double pes_r_min;  // lower limit of r, should be sufficiently large to remove bound states.
  double pes_k_min;  // lower limit of k
  double pes_k_max;  // upper limit of k
  double pes_k_step; // k grid step
  int pes_numk; // number of k grids
  int pes_llr;  // lower limit of r grid index
  int pes_ulr;  // upper limit of r grid index

  std::vector<dcomplex> pes_psik; // k space wavefunctions:    psik[i][l][k];
  std::vector<dcomplex> pes_rhok; // k space orbital products: rhok[i][j][k];
  std::vector<double> pes_bess; // spherical bessel function of the first kind

  clpes();
  ~clpes();
  clpes(const clmpi&, const clio&, const clbas&, const clhprod&);
  void gen(const clmpi&, const clio&, const clbas&, const clhprod&);
  void spec1_k(const clmpi&, const clio&, const clbas&, const clhprod&, const clwfn&) const;
  void spec1_kz(const clmpi&, const clio&, const clbas&, const clhprod&, const clwfn&) const;
private:
  void gen_bess(const clmpi&, const clbas&, const clhprod&);
};
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Main routines and utilities
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void print_date();

int get_abs(int);
double get_abs(double);
double get_abs(dcomplex);

void guess(const clmpi&, const clio&, const clbas&, clwfn&);
double guess_rfun_hlike(int n, int l, const int& znuc, const double& xrad, const double& wrad);
dcomplex guess_rfun_hlike_ecs(int n, int l, const int& znuc, const double& xrad, const dcomplex& cwrad, const clbas& Bas);
double guess_rfun_slater(double neff, double xi, double xrad, double wrad);
dcomplex guess_rfun_slater_ecs(double neff, double xi, double xrad, dcomplex cwrad, const clbas& Bas);
void guess_g09(const clmpi&, const clio&, const clbas&, clwfn&);
void guess_card(const clmpi&, const clio&, const clbas&, clwfn&);
void guess_card_getnlm(const clio&, int ifun, int& n, int& l, int& m);
void guess_card1(const clmpi&, const clio&, const clbas&, clwfn&);
void guess_card2(const clmpi&, const clio&, const clbas&, clwfn&);
void guess_card2_getnlm(const clio&, int ifun, int& n, int& l, int& m);
void guess_card3(const clmpi&, const clio&, const clbas&, clwfn&);
void guess_card4(const clmpi&, const clio&, const clbas&, clwfn&);
void guess_hforb(const clmpi&, const clio&, const clbas&, clwfn&);
void guess_hforb_read(const int& n, const int& l, const int& m, const int& nrad, std::vector<dcomplex>& chirad);
void guess_aufbau(const clmpi&, const clio&, const clbas&, clwfn&);
void guess_aufbau_getnlm(int ifun, int& n, int& l, int& m);
void guess_aufbau_getlm(int lm, int& l, int& m);
void guess_aufbau_getneff(int n, int l, int m, int nel, double znuc, double& neff, double& xi);

void init(const clmpi&, const clio&, const clbas&, clwfn&);
void init_prop12(const clmpi&, const clio&, const clbas&, clfield&, clh2prop&, clhprod&, clwfn&);	      
void init_split2(const clmpi&, const clio&, const clbas&, clfield&, clh1prop&, clh2prop&, clhprod&, clwfn&);
void init_lawrk(const clmpi&, const clio&, const clbas&, clfield&, cllawrk&, clhprod&, clwfn&);

void tdse(const clmpi&, const clio&, const clbas&, clwfn&);
void tdse_prop12(const clmpi&, const clio&, const clbas&, clfield&, clh2prop&, clhprod&, clwfn&);
void tdse_split2(const clmpi&, const clio&, const clbas&, clfield&, clh1prop&, clh2prop&, clhprod&, clwfn&);
void tdse_split4(const clmpi&, const clio&, const clbas&, clfield&, clhprod&, clwfn&);
void tdse_split4x(const clmpi&, const clio&, const clbas&, clfield&, 
		  clh1prop&, clh1prop&, clh2prop&, double dt1, double dt2, clhprod&, clwfn&);
void tdse_read(const clmpi&, const clio&, const clbas&, clfield&, clhprod&, clwfn&);
void tdse_print(const clmpi&, const clio&, const clbas&, const clfield&, clhprod&, const clwfn&);
//void rtdisp_print(const clmpi&, const clio&, const clbas&, const clfield&, clhprod&, const clwfn&);
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Global scope
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
extern double aval_prev[10];
extern double time_prev[10];
const bool CTRUE = true;
const bool CFALSE = false;
const int LZERO  = 0;
const int LONE   = 1;
const int LTWO   = 2;
const int LTHREE = 3;
const int LFOUR  = 4;
const int LFIVE  = 5;
const double ZERO  = 0.0;
const double ONE   = 1.0;
const double TWO   = 2.0;
const double THREE = 3.0;
const double FOUR  = 4.0;
const double FIVE  = 5.0;
const double SIX   = 6.0;
const double EIGHT = 8.0;
const double HALF  = 0.5;
const double THIRD = ONE / THREE;
const double SIXTH = ONE / SIX;
const double PI = 3.14159265358979323846264338327950288;
const dcomplex CZERO(ZERO, ZERO);
const dcomplex RUNIT(ONE, ZERO);
const dcomplex CTWO(TWO, ZERO);
const dcomplex CTHREE(THREE, ZERO);
const dcomplex CFOUR(FOUR, ZERO);
const dcomplex CFIVE(FIVE, ZERO);
const dcomplex CHALF(HALF, ZERO);
const dcomplex CTHIRD(THIRD, ZERO);
const dcomplex IUNIT(ZERO, ONE);
const dcomplex ITWO(ZERO, TWO);
////////////////////////////////////////////////////////////////////////
#endif
