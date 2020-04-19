////////////////////////////////////////////////////////////////////////
// Laser field
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
double aval_pre = 0;
double time_pre = 0;
////////////////////////////////////////////////////////////////////////
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

  IO.read_info("env_type_XUV", "sin2", env_type);
  IO.read_info("delay_XUV", ZERO, delay);
  IO.read_info("cep_XUV", ZERO, cep);
  IO.read_info("fint_XUV", ZERO, fint);
  IO.read_info("wlen_XUV", 11.7, wlen);
  IO.read_info("env_type_NIR", "sin2", env_type2);
  IO.read_info("delay_NIR", ZERO, delay2);
  IO.read_info("cep_NIR", ZERO, cep2);
  IO.read_info("fint_NIR", ZERO, fint2);
  IO.read_info("wlen_NIR", 750.0, wlen2);

  if (env_type.compare("trape") == 0) {
    double F13 = ONE / THREE;
    IO.read_info("cyc1_XUV", F13, cyc1);
    IO.read_info("cyc2_XUV", F13*TWO, cyc2);
    IO.read_info("cyc3_XUV", F13*THREE, cyc3);
  } else if (env_type.compare("sin1") == 0) {
    IO.read_info("sin1cyc_XUV", ONE, cyc3);
  } else if (env_type.compare("sin2") == 0) {
    IO.read_info("sin2cyc_XUV", ONE, cyc3);
  } else if (env_type.compare("sin2flat") == 0) {
    double F13 = ONE / THREE;
    IO.read_info("sin2cyc1_XUV", F13, cyc1);
    IO.read_info("sin2cyc2_XUV", F13*TWO, cyc2);
    IO.read_info("sin2cyc3_XUV", F13*THREE, cyc3);
  } else if (env_type.compare("gauss") == 0) {
    IO.read_info("FWHM_XUV", 200.0, FWHM);
  } else if (env_type.compare("none") == 0) {

  } else {
    std::cout << "clfield::read_info: bad env_type." << std::endl;
    abort();
  }

  if (env_type2.compare("trape") == 0) {
    double F13 = ONE / THREE;
    IO.read_info("cyc1_NIR", F13, cyc12);
    IO.read_info("cyc2_NIR", F13*TWO, cyc22);
    IO.read_info("cyc3_NIR", F13*THREE, cyc32);
  } else if (env_type2.compare("sin1") == 0) {
    IO.read_info("sin1cyc_NIR", ONE, cyc32);
  } else if (env_type2.compare("sin2") == 0) {
    IO.read_info("sin2cyc_NIR", ONE, cyc32);
  } else if (env_type2.compare("sin2flat") == 0) {
    double F13 = ONE / THREE;
    IO.read_info("sin2cyc1_NIR", F13, cyc12);
    IO.read_info("sin2cyc2_NIR", F13*TWO, cyc22);
    IO.read_info("sin2cyc3_NIR", F13*THREE, cyc3);
  } else if (env_type2.compare("gauss") == 0) {
    IO.read_info("FWHM_NIR", 3300.0, FWHM2);
  } else if (env_type.compare("none") == 0) {

  } else {
    std::cout << "clfield::read_info: bad env_type." << std::endl;
    abort();
  }

  IO.read_info("nstep", LONE*10000, nstep);
  IO.read_info("cyctot", ONE, cyctot);

  env_XUV.get(env_type, fint, wlen, cep, delay, FWHM);
  env_NIR.get(env_type2, fint2, wlen2, cep2, delay2, FWHM2);
  env_XUV.tcenter = env_NIR.tcenter + env_NIR.delay;

  if (gauge.compare("length2") == 0 ||
      gauge.compare("velocity1") == 0) {
    env_XUV.famp = env_XUV.famp / env_XUV.freq ;
    env_NIR.famp = env_NIR.famp / env_NIR.freq ;
  }

  /*  if (env_type.compare("gauss") == 0 && env_type2.compare("gauss") == 0 ){
  } else {
    std::cout << "clenv:: Under construction!!" << std::endl;
    abort();
    }*/

  dtime = env_XUV.period / nstep;
  period = env_XUV.period;

  //estimate r_max
  /*    for(int i = 0; i < 100; i++){
    cyctot = 10 * (i + 1);
  double margin = 1.6;
  double gnu = freq / (2 * PI) ;
  double eps_2p = 2080.7 / 2625.49962;
  double v_e = sqrt(2 * (2 * PI * gnu - eps_2p));
  double r_max = v_e * (cyctot / gnu - tcenter) * margin ;
  std::cout << "When cyctot = " << cyctot << ", the estimated value of rmax is " << r_max << std::endl;
  }*/
    // abort();  
  
}

////////////////////////////////////////////////////////////////////////
double clfield::ncyc() const
{
  return time / env_XUV.period;
}
///////////////////////////////////////////////////////////////////////
bool clfield::finished() const
{
    return time > cyctot * env_XUV.period + dtime / TWO;
}
////////////////////////////////////////////////////////////////////////
void clfield::get_value(double* lfield) const
{
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
  /* double cyc, env;

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

  } else {
    std::cout << "clfield::get_evalue_explicit: bad env_type." << std::endl;
    abort();
  }
  */

  lfield[0] = ZERO;
  lfield[1] = ZERO;
  lfield[2] = env_NIR.field(ttmp) + env_XUV.field(ttmp);
}
////////////////////////////////////////////////////////////////////////
void clfield::get_evalue_implicit(double ttmp, double* lfield) const
{
  double cyc, env;
  lfield[0] = ZERO;
  lfield[1] = ZERO;
  lfield[2] = ZERO;

/*  if (td_type == 0 || env_type.compare("zero") == 0) {
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
  } else {
    std::cout << "clfield::get_evalue_implicit: bad env_type." << std::endl;
    abort();
  }
*/
    lfield[2] = - env_NIR.dfield_dt(ttmp) - env_XUV.dfield_dt(ttmp) ;
  
}
////////////////////////////////////////////////////////////////////////
void clfield::get_avalue_explicit(double ttmp, double* lfield) const
{
  double cyc, env;

/*  if (td_type == 0 || env_type.compare("zero") == 0) {
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
  } else {
    std::cout << "clfield::get_avalue_explicit: bad env_type." << std::endl;
    abort();
  }
*/ 
  lfield[0] = ZERO;
  lfield[1] = ZERO;
  lfield[2] = env_NIR.field(ttmp) + env_XUV.field(ttmp);
 
}
////////////////////////////////////////////////////////////////////////
void clfield::get_avalue_implicit(double ttmp, double* lfield) const
{
  double cyc, env;
  lfield[0] = ZERO;
  lfield[1] = ZERO;
  lfield[2] = ZERO;

  /*  if (td_type == 0 || env_type.compare("zero") == 0) {
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
  } else if (env_type.compare("gauss") == 0) {
      get_avalue_implicit_numint(ttmp, lfield);
  } else {
    std::cout << "clfield::get_avalue_implicit: bad env_type." << std::endl;
    abort();
    }*/

  get_avalue_implicit_numint(ttmp, lfield);

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
clenv::clenv(){}
//////////////////////////////////////////////////////////////////////////
void clenv::get(std::string env_type1, double fint1, double wlen1, double cep1, double delay1, double FWHM1)
{
    cep = cep1 / 180.0 * PI;
    famp = sqrt(fint1 / 3.50944506E+16);
    freq = 1239.84190 / (wlen1 * 27.2113845);
    period = TWO * PI / freq;
    FWHM = FWHM1 / 24.188843265;
    sigma = FWHM / 2.35482;
    tcenter = sigma * 4;
    delay = delay1 / 24.188843265;
    env_type = env_type1;
}
////////////////////////////////////////////////////////////////////////
double clenv::field(double ttmp)const
{
  double t2 = ttmp - delay;
  double result;
  double env;
  if (env_type.compare("gauss") == 0) {
    env = exp(- (t2 - tcenter) * (t2 - tcenter)/(2 * sigma * sigma));
    result = env * famp * sin(freq * t2 + cep); 
  } else if (env_type.compare("none") == 0) {
    result = 0.0;
  }
  return result;
}
////////////////////////////////////////////////////////////////////////
double clenv::dfield_dt(double ttmp)const
{
  double t2 = ttmp - delay;
  double result;
  double env;
  if (env_type.compare("gauss") == 0) {
    env = exp(- (t2 - tcenter) * (t2 - tcenter)/(2 * sigma * sigma));
    result = env * famp * ( freq * cos(freq * t2 + cep) - (t2 - tcenter) * sin(freq * t2 + cep) / (sigma * sigma)  ) ; 
  } else if (env_type.compare("none") == 0) {
    result = 0.0;
  }
  return result;
}
////////////////////////////////////////////////////////////////////////
