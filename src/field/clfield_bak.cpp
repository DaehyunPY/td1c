////////////////////////////////////////////////////////////////////////
// Laser field
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
double aval_pre = 0;
double time_pre = 0;
////////////////////////////////////////////////////////////////////////
extern "C"
{
}
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

  IO.read_info("cep", ZERO, cep);
  IO.read_info("fint", ZERO, fint);
  IO.read_info("wlen", 800.0, wlen);
  IO.read_info("env_type", "sin2", env_type);

  if (env_type.compare("zero") == 0) {
  } else if (env_type.compare("trape") == 0) {
    double F13 = ONE / THREE;
    IO.read_info("cyc1", F13, cyc1);
    IO.read_info("cyc2", F13*TWO, cyc2);
    IO.read_info("cyc3", F13*THREE, cyc3);
  } else if (env_type.compare("sin1") == 0) {
    IO.read_info("sin1cyc", ONE, cyc3);
  } else if (env_type.compare("sin2") == 0) {
    IO.read_info("sin2cyc", ONE, cyc3);
  } else if (env_type.compare("sin2flat") == 0) {
    double F13 = ONE / THREE;
    IO.read_info("sin2cyc1", F13, cyc1);
    IO.read_info("sin2cyc2", F13*TWO, cyc2);
    IO.read_info("sin2cyc3", F13*THREE, cyc3);
  } else if (env_type.compare("gauss") == 0) {
    IO.read_info("FWHM", ONE, FWHM_as);
    FWHM = FWHM_as / 24.188843265;
    sigma = FWHM / 2.35482;
    //too long...     tcenter = sigma * 8;
    //too short...    tcenter = sigma * 4;
    tcenter = sigma * 6;
  } else {
    std::cout << "clfield::read_info: bad env_type." << std::endl;
    abort();
  }

  cep = cep / 180.0 * PI;
  famp = sqrt(fint / 3.50944506E+16);
  freq = 1239.84190 / (wlen * 27.2113845);
  period = TWO * PI / freq;

  IO.read_info("nstep", LONE*10000, nstep);
  IO.read_info("cyctot", ONE, cyctot);
  dtime = period / nstep;

}
////////////////////////////////////////////////////////////////////////
double clfield::ncyc() const
{
  return time / period;
}
////////////////////////////////////////////////////////////////////////
bool clfield::finished() const
{
    return time > cyctot * period + dtime / TWO;
}
////////////////////////////////////////////////////////////////////////
void clfield::get_value(double* lfield) const
{
  get_value(time, lfield);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_value(double t1, double t2, double c1, double c2, double* lfield) const
{
  double lfield1[9];
  double lfield2[9];
  get_value(t1, lfield1);
  get_value(t2, lfield2);
  for (int i = 0; i < 9; i++) {
    lfield[i] = c1 * lfield1[i] + c2 * lfield2[i];
  }
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
void clfield::get_value(double ttmp, double* lfield) const
{
  if (gauge.compare("length1") == 0) {
    get_evalue_explicit(ttmp, lfield);
  } else if (gauge.compare("length2") == 0) {
    get_evalue_implicit(ttmp, lfield);
  } else if (gauge.compare("velocity1") == 0) {
    get_avalue_explicit(ttmp, lfield);
  } else if (gauge.compare("velocity2") == 0) {
    get_avalue_implicit(ttmp, lfield);
  } else {
    std::cout << "clfield::get_value. bad gauge" << std::endl;
  }

  get_evalue(ttmp, &lfield[3]);
  get_avalue(ttmp, &lfield[6]);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_der(double* dfield) const
{
  get_der(time, dfield);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_der(double ttmp, double* dfield) const
{
  if (gauge.compare("length1") == 0) {
    get_eder_explicit(ttmp, dfield);
  } else if (gauge.compare("length2") == 0) {
    get_eder_implicit(ttmp, dfield);
  } else if (gauge.compare("velocity1") == 0) {
    get_ader_explicit(ttmp, dfield);
  } else if (gauge.compare("velocity2") == 0) {
    get_ader_implicit(ttmp, dfield);
  } else {
    std::cout << "clfield::get_der. bad gauge" << std::endl;
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_evalue(double ttmp, double* efield) const
{
  if (gauge.compare("length1") == 0) {
    get_evalue_explicit(ttmp, efield);
  } else if (gauge.compare("length2") == 0) {
    get_evalue_implicit(ttmp, efield);
  } else if (gauge.compare("velocity1") == 0) {
    get_evalue_implicit(ttmp, efield);
  } else if (gauge.compare("velocity2") == 0) {
    get_evalue_explicit(ttmp, efield);
  } else {
    std::cout << "clfield::get_evalue. bad gauge" << std::endl;
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_avalue(double ttmp, double* afield) const
{
  if (gauge.compare("length1") == 0) {
    get_avalue_implicit(ttmp, afield);
  } else if (gauge.compare("length2") == 0) {
    get_avalue_explicit(ttmp, afield);
  } else if (gauge.compare("velocity1") == 0) {
    get_avalue_explicit(ttmp, afield);
  } else if (gauge.compare("velocity2") == 0) {
    get_avalue_implicit(ttmp, afield);
  } else {
    std::cout << "clfield::get_avalue. bad gauge" << std::endl;
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_evalue_explicit(double ttmp, double* lfield) const
{
  double cyc, env;

  if (td_type == 0 || env_type.compare("zero") == 0) {
    env = ZERO;
  } else if (env_type.compare("trape") == 0) {
    cyc = ttmp / period;
    if (cyc < cyc1) {
  	env = cyc / cyc1;
    } else if (cyc < cyc2) {
  	env = ONE;
    } else if (cyc < cyc3) {
  	env = (cyc - cyc3) / (cyc2 - cyc3);
    } else {
  	env = ZERO;
    }
  } else if (env_type.compare("sin2") == 0) {
     cyc = ttmp / period;
     if (cyc < cyc3) {
  	env = sin(PI * ttmp / (cyc3 * period));
  	env = env * env;
    } else {
  	env = ZERO;
	}
	 // env = ONE;
  } else if (env_type.compare("sin2flat") == 0) {
    double ttime;
    double sin2cyc;
    cyc = ttmp / period;
    if (cyc < cyc1) {
  	ttime = ttmp;
  	sin2cyc = cyc1 * TWO;
  	env = sin(PI * ttime / (sin2cyc * period));
  	env = env * env;
    } else if (cyc < cyc2) {
  	env = ONE;
    } else if (cyc < cyc3) {
  	sin2cyc = (cyc3 - cyc2) * TWO;
  	ttime = ttmp - cyc2 * period + (cyc3 - cyc2) * period;
  	env = sin(PI * ttime / (sin2cyc * period));
  	env = env * env;
    } else {
  	env = ZERO;
    }
  } else if (env_type.compare("gauss") == 0) {
    if (ttmp < tcenter * 2){
      env = exp(-(ttmp-tcenter)*(ttmp-tcenter)/(2*sigma*sigma));
    } else {
      env = ZERO;
    }
  } else {
    std::cout << "clfield::get_evalue_explicit: bad env_type." << std::endl;
    abort();
  }
  lfield[0] = ZERO;
  lfield[1] = ZERO;
  lfield[2] = env * famp * sin(freq * ttmp + cep);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_eder_explicit(double ttmp, double* dfield) const
{
  double cyc, env, denv, fac;

  if (td_type == 0 || env_type.compare("zero") == 0) {
    env = ZERO;
    denv = ZERO;
  } else if (env_type.compare("trape") == 0) {
    std::cout << "get_eder_explicit: env_type = trape nyi." << std::endl;
    abort();
  } else if (env_type.compare("sin2") == 0) {
     cyc = ttmp / period;
     if (cyc < cyc3) {
       fac = PI / (cyc3 * period);
       env = sin(fac*ttmp);
       env = env * env;
       denv = TWO * fac * sin(fac*ttmp) * cos(fac*ttmp);
     } else {
       env = ZERO;
       denv = ZERO;
     }
  } else if (env_type.compare("sin2flat") == 0) {
    std::cout << "get_eder_explicit: env_type = sin2flat nyi." << std::endl;
    abort();
  } else if (env_type.compare("gauss") == 0) {
    std::cout << "get_eder_explicit: env_type = gauss nyi." << std::endl;
    abort();
  } else {
    std::cout << "clfield::get_eder_explicit: bad env_type." << std::endl;
    abort();
  }
  dfield[0] = ZERO;
  dfield[1] = ZERO;
  dfield[2] = famp * (freq * env * cos(freq * ttmp + cep) + denv * sin(freq * ttmp + cep));
}
////////////////////////////////////////////////////////////////////////
void clfield::get_eder_implicit(double ttmp, double* dfield) const
{
  double cyc, env, denv, fac;

  dfield[0] = ZERO;
  dfield[1] = ZERO;
  dfield[2] = ZERO;
  if (td_type == 0 || env_type.compare("zero") == 0) {
    env = ZERO;
    denv = ZERO;
  } else if (env_type.compare("gauss") == 0) {
    if (ttmp < tcenter * 2){
      env = exp(-(ttmp-tcenter)*(ttmp-tcenter)/(2*sigma*sigma));
    } else {
      env = ZERO;
    }
    double sig2 = sigma * sigma;
    double sig4 = sig2 * sig2;
    double tdiff = ttmp - tcenter;
    dfield[2] = env * famp * (cos(freq*ttmp + cep) * TWO * tdiff / sig2
			    + sin(freq*ttmp + cep) * (sig4*freq*freq + sig2 - tdiff*tdiff) / (sig4 * freq));
  } else {
    std::cout << "clfield::get_eder_implicit: bad env_type." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_ader_explicit(double ttmp, double* dfield) const
{
  get_evalue_implicit(ttmp, dfield);
  dfield[0] = -dfield[0];
  dfield[1] = -dfield[1];
  dfield[2] = -dfield[2];
}
////////////////////////////////////////////////////////////////////////
void clfield::get_ader_implicit(double ttmp, double* dfield) const
{
  get_evalue_explicit(ttmp, dfield);
  dfield[0] = -dfield[0];
  dfield[1] = -dfield[1];
  dfield[2] = -dfield[2];
}
////////////////////////////////////////////////////////////////////////
void clfield::get_evalue_implicit(double ttmp, double* lfield) const
{
  double cyc, env;
  lfield[0] = ZERO;
  lfield[1] = ZERO;
  lfield[2] = ZERO;

  if (td_type == 0 || env_type.compare("zero") == 0) {
  } else if (env_type.compare("trape") == 0) {
    double phase;
    double slope;
    cyc = ttmp / period;
    phase = freq * ttmp + cep;
    if (cyc < cyc1) {
      slope = ONE / (cyc1 * period);
      lfield[2] = famp * slope * (cos(phase) * ttmp + sin(phase) / freq);
    } else if (cyc < cyc2) {
      lfield[2] = famp * cos(phase);
    } else if (cyc < cyc3) {
      slope = ONE / ((cyc2 - cyc3) * period);
      lfield[2] = famp * slope * (cos(phase) * (ttmp - cyc3 * period) + sin(phase) / freq);
    } else {
      lfield[2] = ZERO;
    }
  } else if (env_type.compare("sin2") == 0) {
    double fac0, facp, facm;
    double wpls = (ONE + ONE / cyc3) * freq;
    double wmns = (ONE - ONE / cyc3) * freq;
    cyc = ttmp / period;
    if (cyc < cyc3) {
  	fac0 = cos(freq * ttmp + cep);
  	facp = cos(wpls * ttmp + cep) * wpls / freq * HALF;
  	facm = cos(wmns * ttmp + cep) * wmns / freq * HALF;
  	lfield[2] = - HALF * famp * (fac0 - facp - facm);
    } else {
  	lfield[2] = ZERO;
    }
  } else if (env_type.compare("gauss") == 0) {
    if (ttmp < tcenter * 2){
      env = exp(-(ttmp-tcenter)*(ttmp-tcenter)/(2*sigma*sigma));
    } else {
      env = ZERO;
    }
    lfield[2] = env * famp * (sin(freq*ttmp + cep) * (ttmp - tcenter)/(freq*sigma*sigma)
			    - cos(freq*ttmp + cep));
  } else {
    std::cout << "clfield::get_evalue_implicit: bad env_type." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_avalue_explicit(double ttmp, double* lfield) const
{
  double cyc, env;

  if (td_type == 0 || env_type.compare("zero") == 0) {
    env = ZERO;
  } else if (env_type.compare("trape") == 0) {
    cyc = ttmp / period;
    if (cyc < cyc1) {
  	env = cyc / cyc1;
    } else if (cyc < cyc2) {
  	env = ONE;
    } else if (cyc < cyc3) {
  	env = (cyc - cyc3) / (cyc2 - cyc3);
    } else {
  	env = ZERO;
    }
  } else if (env_type.compare("sin2") == 0) {
    cyc = ttmp / period;
    if (cyc < cyc3) {
  	env = sin(PI * ttmp / (cyc3 * period));
  	env = env * env;
    } else {
  	env = ZERO;
    }
  } else if (env_type.compare("sin2flat") == 0) {
    double ttime;
    double sin2cyc;
    cyc = ttmp / period;
    if (cyc < cyc1) {
  	ttime = ttmp;
  	sin2cyc = cyc1 * TWO;
  	env = sin(PI * ttime / (sin2cyc * period));
  	env = env * env;
    } else if (cyc < cyc2) {
  	env = ONE;
    } else if (cyc < cyc3) {
  	sin2cyc = (cyc3 - cyc2) * TWO;
  	ttime = ttmp - cyc2 * period + (cyc3 - cyc2) * period;
  	env = sin(PI * ttime / (sin2cyc * period));
  	env = env * env;
    } else {
  	env = ZERO;
    }
  } else if (env_type.compare("gauss") == 0) {
    if (ttmp < tcenter * 2){
      env = exp(-(ttmp-tcenter)*(ttmp-tcenter)/(2*sigma*sigma));
    } else {
      env = ZERO;
    }
  } else {
    std::cout << "clfield::get_avalue_explicit: bad env_type." << std::endl;
    abort();
  }
  lfield[0] = ZERO;
  lfield[1] = ZERO;
  lfield[2] = env * famp / freq * sin(freq * ttmp + cep);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_avalue_implicit(double ttmp, double* lfield) const
{
  double cyc, env;
  lfield[0] = ZERO;
  lfield[1] = ZERO;
  lfield[2] = ZERO;

  if (td_type == 0 || env_type.compare("zero") == 0) {
  } else if (env_type.compare("trape") == 0) {
    double slope, phase, p1, p2;
    cyc = ttmp / period;
    phase = freq * ttmp + cep;
    p1 = freq * cyc1 * period + cep;
    p2 = freq * cyc2 * period + cep;
    if (cyc < cyc1) {
      slope = ONE / (cyc1 * period);
      lfield[2] = famp / freq * slope * (cos(phase) * ttmp - (sin(phase) - sin(cep)) / freq);
    } else if (cyc < cyc2) {
      lfield[2] = famp / freq * (cos(phase) - (sin(p1) - sin(cep)) / freq);
    } else if (cyc < cyc3) {
      slope = ONE / ((cyc2 - cyc3) * period);
      lfield[2] = famp / freq * (
				 slope * (cos(phase) * (ttmp - cyc3 * period) - (sin(phase) - sin(p2)) / freq)
				 - (sin(p1) - sin(cep)) / freq
				);
    } else {
      lfield[2] = ZERO;
    }
  } else if (env_type.compare("sin2") == 0) {
    double fac0, facp, facm;
    double wpls = (ONE + ONE / cyc3) * freq;
    double wmns = (ONE - ONE / cyc3) * freq; 
    cyc = ttmp / period;
    if (cyc < cyc3) {
      fac0 = (cos(freq * ttmp + cep) - cos(cep)) / freq;
      facp = (cos(wpls * ttmp + cep) - cos(cep)) / wpls * HALF;
      facm = (cos(wmns * ttmp + cep) - cos(cep)) / wmns * HALF;
      lfield[2] = HALF * famp * (fac0 - facp - facm);
    } else {
      lfield[2] = ZERO;
    }
    //    get_avalue_implicit_numint(ttmp, lfield);
    //    lfield[2] = (cos(freq * ttmp + cep) - cos(cep)) / freq;
  } else if (env_type.compare("gauss") == 0) {
      get_avalue_implicit_numint(ttmp, lfield);
  } else {
    std::cout << "clfield::get_avalue_implicit: bad env_type." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clfield::get_avalue_implicit_numint(double ttmp, double* lfield) const
{

  double time_d = ttmp - time_pre ;
  double field_now[3];

  double time_now = ttmp - time_d;

  get_evalue_explicit(time_now, field_now);

  double k1 = - time_d * field_now[2];

  time_now = ttmp - time_d / 2.0;
  
  get_evalue_explicit(time_now, field_now);

  double k2 = - time_d * field_now[2];

  time_now = ttmp;

  get_evalue_explicit(time_now, field_now);

  double k4 = - time_d * field_now[2];

  lfield[2] = aval_pre + (k1 + 2 * k2 + 2 * k2 + k4) / 6.0;

  if(ttmp < dtime){
    lfield[2] = ZERO;
  }

  aval_pre = lfield[2];
  time_pre = ttmp;
}
//////////////////////////////////////////////////////////////////////////
