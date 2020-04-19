#ifndef SURFF_INCLUDED
#define SURFF_INCLUDED

//prototype declaration
#include <string>
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
#include <stdarg.h>

class clmpi;
class clio;
class clbas;
class clwfn;
class clhprod;
class clfield;
class clormas;
typedef std::complex<double> dcomplex;
typedef std::vector<dcomplex* > type_cic;
typedef std::vector<std::vector<dcomplex* > > type_orb;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// MKL lapack
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
#define MKL_Complex16 std::complex<double>
// Sato_TSURFF
// #include <mkl_lapack.h>
// //#include <mkl_lapacke.h>
// Sato_TSURFF


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// benri
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
static std::string strsprintf(const char* format,...);
std::string strsprintf(const char* format,...){
  va_list ap;
  va_start(ap, format);
  char* alloc;
  if(vasprintf(&alloc,format,ap) == -1) {
    return std::string("");
  }
  va_end(ap);
  std::string retStr = std::string(alloc);
  free(alloc);
  return retStr;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// integral of surface flux (t-surff)
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class Integ
{
protected:
  int ksize;
  int norb;

public:
  Integ(){};
  Integ(int vir_ksize, int vir_norb){
    std::cout << "# Integ()" << std::endl;;
    ksize = vir_ksize;
    norb = vir_norb;
  };

  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
		    std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes) = 0;
  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
		    const std::vector<dcomplex> &den1, std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes){};
  
  //for EI
  virtual void set_renormarize(int iren){};

};
/////////////////////////////////////
class Integ_trap : public Integ
{
private:
  std::vector<dcomplex> tmp_opes;

public:
  Integ_trap(int vir_ksize, int vir_norb);
  ~Integ_trap();
  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
	    std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes);
  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
		    const std::vector<dcomplex> &den1, std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes);

};
/////////////////////////////////////
class Integ_simpson : public Integ
{
public:
  Integ_simpson();
  ~Integ_simpson();
  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
	    std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes);
  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
		    const std::vector<dcomplex> &den1, std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes);

};
/////////////////////////////////////
class Integ_CN : public Integ
{
private:
  std::vector<dcomplex> tmp_opes;
  std::vector<dcomplex> exp_vmat;
  std::vector<dcomplex> imp_vmat;

  std::vector<dcomplex> old_v2xmat;
  std::vector<dcomplex> old_opes_dt;

public:
  Integ_CN(int vir_ksize, int vir_norb);
  ~Integ_CN();
  void rec_old_state(const clmpi&, const std::vector<dcomplex> &v2xmat, std::vector<dcomplex> &opes_dt);
  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
		    std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes);
  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
		    const std::vector<dcomplex> &den1, std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes);
  virtual void prop_old(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
	    std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes);

};
/////////////////////////////////////
class Integ_EI : public Integ
{
private:
  std::vector<dcomplex> tmp_opes;
  int renormarize;
  

public:
  Integ_EI(int vir_ksize, int vir_norb);
  ~Integ_EI();

  void diag_xmat();
  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
		    std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes);
  virtual void prop(const clmpi&, const clbas&, const clfield&, const std::vector<dcomplex> &v2xmat,
		    const std::vector<dcomplex> &den1, std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes);
  void eigen_decompose(double dt, const std::vector<dcomplex> &v2xmat, std::vector<dcomplex> &eval,  
		       std::vector<dcomplex> &evl,  std::vector<dcomplex> &evr);
  void inv_mat(std::vector<dcomplex> &inv_mat);
  void calc_pade1(std::vector<dcomplex> &eval, std::vector<dcomplex> &pade);
  void calc_pade2(std::vector<dcomplex> &eval, std::vector<dcomplex> &pade);
  
  virtual void set_renormarize(int iren);

  
};
/////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// time dependent surface flux (t-surff)
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
class Surff
{
    
public:
  int nprint;
  std::string integ_type;
  int num_krad;
  int num_kang;
  int ksize;
  int norb;
  int lnum;
  int mnum;
  double norm;
  std::vector<double> krad;
  std::vector<double> kang;
  std::vector<double> wang;
  std::vector<dcomplex> pes1e_surff;
  std::vector<double> sph_ang;
  std::vector<dcomplex> sbessel;
  std::vector<dcomplex> sbessel_dr;
  std::vector<double> alph_lm; //<Y_lm|cos(theta)|Y_l'm'>
  std::vector<double> vphase; // volkov phase
  std::vector<double> tdfac; // time dependent factor of volkov phase
  std::vector<dcomplex> v2xmat;
  std::vector<dcomplex> orb_surf;
  std::vector<dcomplex> dorb_surf;
  std::vector<dcomplex> orbpes;
  std::vector<dcomplex> orbpes_dt;
  std::vector<dcomplex> torbpes;

  std::vector<dcomplex> rhok;
  std::vector<dcomplex> ang_spec;
  std::vector<double> mom_spec;
  std::vector<double> ene_spec;
    
  Integ* integ;
  Surff(const clmpi&, const clio&, const clbas&);
  ~Surff();
  void prop(const clmpi&, const clwfn&, const clbas&, const clfield&, const clhprod&);
  void calc_me_spectrum(const clmpi&, const clwfn&, const clbas&, clhprod&);    
  void rec_v2xmat(double dtime, const std::vector<dcomplex>& xmat);
  void rec_v2xmat(double dtime);
  void print(const clio&);
  void print_file(const clio&);
  void print_amplitude(const clio&);
    
private:
  int irad_surf;
  double rrad_surf;
  void init(const clmpi&, const clio&, const clbas&);
  void init_kang(const clmpi&, const clio&);
  void calc_dt(const clmpi&, const clwfn&, const clbas&, const clfield&, std::vector<dcomplex> &tmppes_dt);
  void prop_trap(const clmpi&, const clbas&, const clfield&, std::vector<dcomplex> &tmppes_dt, std::vector<dcomplex> &tmppes);
  void prop_trap_remove_Q(const clmpi&, const clbas&, const clfield&, std::vector<dcomplex> &tmppes_dt, std::vector<dcomplex> &tmppes);
  void prop_trap_vphase(const clmpi& MPIP, double afield, double time, double dtime);
  void calc_orb(const clbas &Bas, const clwfn &Wfn, 
	  	int srad, std::vector<dcomplex> &norb, std::vector<dcomplex> &dorb);

  void calc_momentum_spectrum(const clmpi&, const dcomplex*, const clbas&, clhprod &HPW);
  void calc_momentum_spectrum_1e(const clmpi&, const dcomplex*, const clbas&, clhprod &HPW);
  void convert_energy_spectrum();
  void calc_t_surff(const clmpi&, const clwfn& wfn, const clbas&, const clfield&);
};

#endif
