////////////////////////////////////////////////////////////////////////
// Laser field
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
double time_prev[10]; // This is for clfield1::get_avalue_numint. DO NOT touch 
double aval_prev[10]; // these arrays anywhere except in clfield1::get_avalue_numint.
////////////////////////////////////////////////////////////////////////
clfield::clfield()
{
}
////////////////////////////////////////////////////////////////////////
clfield::~clfield()
{
}
////////////////////////////////////////////////////////////////////////
clfield::clfield(const clmpi& MPIP, const clio& IO)
{
  gen(MPIP, IO);
}
//////////////////////////////////////////////////////////////////////////
void clfield::gen(const clmpi& MPIP, const clio& IO)
{
  if (IO.job_type.compare("td") == 0) {
    td_type = 1;
  } else {
    td_type = 0;
  }

  IO.read_info("gauge", "length1", gauge);
  if (gauge.compare("length1") == 0 ||
      gauge.compare("length2") == 0) {
    lgauge = true;
    vgauge = false;
    igauge = 0;
  } else if (gauge.compare("velocity1") == 0 ||
	     gauge.compare("velocity2") == 0) {
    lgauge = false;
    vgauge = true;
    igauge = 1;
  } else {
    std::cout << "clfield::read_info: bad gauge" << std::endl;
    abort();
  }

  IO.read_info("num_field", LONE, num_field);
  IO.read_info("nstep", 10000, nstep);
  IO.read_info("cyctot", ONE, cyctot);
  IO.read_info("cycinit", ZERO, cycinit);

  Fields.resize(num_field);
  double ttmpl, ttmpr;
  double tinit = ZERO;
  double tlast = ZERO;
  for (long ifield = 0; ifield < num_field; ifield++) {
    //    std::cout << "clfield::gen (1)" << std::endl;
    Fields[ifield].gen1(MPIP, IO, ifield);
    //    std::cout << "clfield::gen (2)" << std::endl;
    ttmpl = Fields[ifield].delay - Fields[ifield].tau;
    ttmpr = Fields[ifield].delay + Fields[ifield].tau;
    if (ttmpl < tinit) tinit = ttmpl;
    if (ttmpr > tlast) tlast = ttmpr;
    time_prev[ifield] = ZERO;
    aval_prev[ifield] = ZERO;
  }
  dtime = Fields[0].period / nstep;

  double tshift = -tinit;
  tlast += tshift;
  printf("clfield::gen: t_last is %20.10f cycles of pulse #0.\n", tlast/Fields[0].period);
  printf("clfield::gen: cyctot is %20.10f cycles of pulse #0.\n", cyctot);

  for (long ifield = 0; ifield < num_field; ifield++) {
    Fields[ifield].gen2(MPIP, IO, tshift);
  }

  if (tlast > cyctot * Fields[0].period) {
    std::cout << "clfield::gen: tlast > total simulation time." << std::endl;
    //    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::init()
{
  step = LZERO;
  time = cycinit * Fields[0].period;
}
////////////////////////////////////////////////////////////////////////
double clfield::ncyc() const
{
  return time / Fields[0].period;
}
////////////////////////////////////////////////////////////////////////
bool clfield::finished() const
{
    return time > cyctot * Fields[0].period + dtime / TWO;
}
////////////////////////////////////////////////////////////////////////
void clfield::get_value(double* lfield) const
{
  //DEBUG
  //  printf("clfield::get_value: time = %20.10f\n", time);
  //DEBUG  
  get_value(time, lfield);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_evalue(double* efield) const
{
  get_evalue(time, efield);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_avalue(double* afield) const
{
  get_avalue(time, afield);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_der(double* dfield) const
{
  get_der(time, dfield);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_value(double ttmp, double* lfield) const
{
  //DEBUG
  //  printf("clfield::get_value: ttmp = %20.10f\n", ttmp);
  //DEBUG  
  for (long i = 0; i < 9; i++) {
    lfield[i] = ZERO;
  }
  for (long ifield = 0; ifield < num_field; ifield++) {
    if (Fields[ifield].env_type.compare("pulse") == 0) {
      //if (ttmp < dtime + 1.E-10) {
      if (step == 0) {
	lfield[2] += Fields[ifield].famp / dtime;
	lfield[5] += Fields[ifield].famp / dtime;
      }
    } else {
      Fields[ifield].get_value(ttmp, lfield);
    }
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_der(double ttmp, double* dfield) const
{
  for (long i = 0; i < 9; i++) {
    dfield[i] = ZERO;
  }
  for (long ifield = 0; ifield < num_field; ifield++) {
    if (Fields[ifield].env_type.compare("pulse") == 0) {
    } else {
      Fields[ifield].get_der(ttmp, dfield);
    }
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_evalue(double ttmp, double* efield) const
{
  for (long i = 0; i < 3; i++) {
    efield[i] = ZERO;
  }
  for (long ifield = 0; ifield < num_field; ifield++) {
    if (Fields[ifield].env_type.compare("pulse") == 0) {
      //if (ttmp < dtime + 1.E-10) {
      if (step == 0) {
	efield[2] += Fields[ifield].famp / dtime;
      }
    } else {
      Fields[ifield].get_evalue(ttmp, efield);
    }
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_avalue(double ttmp, double* afield) const
{
  for (long i = 0; i < 3; i++) {
    afield[i] = ZERO;
  }
  for (long ifield = 0; ifield < num_field; ifield++) {
    if (Fields[ifield].env_type.compare("pulse") == 0) {
    } else {
      Fields[ifield].get_avalue(ttmp, afield);
    }
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_value(double t1, double t2, double c1, double c2, double* lfield) const
{
  for (long i = 0; i < 9; i++) {
    lfield[i] = ZERO;
  }
  for (long ifield = 0; ifield < num_field; ifield++) {
    if (Fields[ifield].env_type.compare("pulse") == 0) {
      //if (t1 < dtime + 1.E-10) {
      if (step == 0) {
	lfield[2] += Fields[ifield].famp * c1 / dtime;
	lfield[5] += Fields[ifield].famp * c1 / dtime;
      }
      //if (t2 < dtime + 1.E-10) {
      if (step == 0) {
	lfield[2] += Fields[ifield].famp * c2 / dtime;
	lfield[5] += Fields[ifield].famp * c2 / dtime;
      }
    } else {
      Fields[ifield].get_value(t1, t2, c1, c2, lfield);
    }
  }
}
//////////////////////////////////////////////////////////////////////////
