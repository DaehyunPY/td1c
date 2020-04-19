////////////////////////////////////////////////////////////////////////
// Laser field
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
clfield1::clfield1()
{
}
////////////////////////////////////////////////////////////////////////
clfield1::~clfield1()
{
}
//////////////////////////////////////////////////////////////////////////
void clfield1::gen1(const clmpi& MPIP, const clio& IO, long ifield)
{
  field_num = ifield;
  std::string key;
  std::string tag = ""; 
  //C++11  if (field_num > 0) tag = std::to_string(field_num);
  std::stringstream ss;
  ss << field_num;
  if (field_num > 0) tag = ss.str();

  key = "gauge";                IO.read_info(key, "length1", gauge);
  key = "env_type"; key += tag; IO.read_info(key, "sin2", env_type);
  key = "fint";     key += tag; IO.read_info(key, ZERO, fint);
  key = "wlen";     key += tag; IO.read_info(key, 800.0, wlen);
  key = "cep";      key += tag; IO.read_info(key, ZERO, cep);
  key = "delay";    key += tag; IO.read_info(key, ZERO, delay);
  key = "numcyc";   key += tag; IO.read_info(key, ONE, numcyc);
  key = "cyc1";  key += tag; IO.read_info(key, -ONE, cyc1);    // for trape/sin2flat
  key = "cyc2";  key += tag; IO.read_info(key, -ONE, cyc2);    // for trape/sin2flat
  key = "cyc3";  key += tag; IO.read_info(key, -ONE, cyc3);    // for trape/sin2flat

  double FWHM_fs, FWHM_as, FWHM_cyc;
  key = "FWHM_cyc"; key += tag; IO.read_info(key, -ONE, FWHM_cyc); // for gauss (cyc->au, highest priority)
  key = "FWHM_fs";  key += tag; IO.read_info(key, -ONE, FWHM_fs);  // for gauss (as->au, the next)
  key = "FWHM_as";  key += tag; IO.read_info(key, -ONE, FWHM_as);  // for gauss (as->au, lowest)

  cep = cep / 180.0 * PI;
  famp = sqrt(fint / 3.50944506E+16);
  freq = 1239.84190 / (wlen * 27.2113845);
  period = TWO * PI / freq;
  delay = delay / 24.188843265;

  if (FWHM_cyc > ZERO) {
    FWHM = FWHM_cyc * period;
  } else if (FWHM_fs > ZERO) {
    FWHM = FWHM_fs / 0.024188843265;
  } else if (FWHM_as > ZERO) {
    FWHM = FWHM_as / 24.188843265;
  }

  //  FWHM for E0
  //  sigma = FWHM / 2.35482;
  //  FWHM for I0 \propto E0*E0
  sigma = sqrt(2.0) * FWHM / 2.35482;

  if (env_type.compare("zero") == 0) {
    tau = ZERO;
  } else if (env_type.compare("pulse") == 0) {
    tau = numcyc * period * HALF;
  } else if (env_type.compare("flat") == 0) {
    tau = numcyc * period * HALF;
  } else if (env_type.compare("trape") == 0) {
    if (cyc1 < ONE || cyc2 < ZERO || cyc3 < ONE || cyc1 + cyc2 + cyc3 != numcyc) {
      std::cout << "bad cyc1-cyc3 and numcyc for env_type trape." << std::endl;
      abort();
    }
    cyc2 += cyc1;
    cyc3 += cyc2;
    tau = numcyc * period * HALF;
  } else if (env_type.compare("sin2") == 0) {
    if (numcyc < ZERO) {
      numcyc = PI * FWHM / sqrt(TWO) / period;
      numcyc = ceil(numcyc);
    }
    tau = numcyc * period * HALF;
  } else if (env_type.compare("sin2flat") == 0) {
    if (cyc1 < ONE || cyc2 < ZERO || cyc3 < ONE || cyc1 + cyc2 + cyc3 != numcyc) {
      std::cout << "bad cyc1-cyc3 and numcyc for env_type sin2flat." << std::endl;
      abort();
    }
    cyc2 += cyc1;
    cyc3 += cyc2;
    tau = numcyc * period * HALF;
  } else if (env_type.compare("sinc") == 0) {
    double freq0_ev, freq1_ev, freq0, freq1;
    if (freq < ZERO) {
      key = "freq0_ev"; key += tag; IO.read_info(key, freq0_ev);
      key = "freq1_ev"; key += tag; IO.read_info(key, freq1_ev);
      freq0 = freq0_ev / 27.2113845;
      freq1 = freq1_ev / 27.2113845;
    } else {
      freq0 = ZERO;
      freq1 = freq;
    }
    freq = (freq1 + freq0) * HALF;
    dfreq = (freq1 - freq0) * HALF;
    wlen = 1239.84190 / (freq * 27.2113845);
    period = TWO * PI / freq;
    tau = numcyc * period * HALF;
  } else if (env_type.compare("gauss") == 0) {
    if (FWHM < ZERO) {
      std::cout << "env_type gauss requires positive FWHM." << std::endl;
      abort();
    }
    // tau = sigma * 2; // far too short...
    tau = sigma * 4; // may be too short... (iwatsu-kun's sotsuron)
    // tau = sigma * 6; // seems appropriate
    // tau = sigma * 8; // may be too long...
  } else {
    std::cout << "clfield1::gen1: bad env_type." << std::endl;
    abort();
  }
}
//////////////////////////////////////////////////////////////////////////
void clfield1::gen2(const clmpi& MPIP, const clio& IO, double tshift)
{
  tcenter = tshift + delay;
  tleft = tcenter - tau;
  tright = tcenter + tau;
  // DEBUG
  printf("clfield1: # %5ld: tau = %20.10f.\n", field_num, tau);
  printf("clfield1: # %5ld: t_L = %20.10f.\n", field_num, tleft);
  printf("clfield1: # %5ld: t_0 = %20.10f.\n", field_num, tcenter);
  printf("clfield1: # %5ld: t_R = %20.10f.\n", field_num, tright);
  // DEBUG
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_value(double time, double* lfield) const
{
  if (env_type.compare("pulse") == 0) {
    if (std::abs(time) < 1.0E-10) {
      lfield[2] += famp;
      lfield[5] += famp;
      //DEBUG
      printf("### clfield1::get_value:  %20.10f%20.10f\n", time, lfield[2]);
      //DEBUG
    }
  } else if (time > tleft && time < tright || env_type.compare("gauss") == 0) {
     if (gauge.compare("length1") == 0) {
       get_evalue_explicit(time, lfield);
     } else if (gauge.compare("length2") == 0) {
       get_evalue_implicit(time, lfield);
     } else if (gauge.compare("velocity1") == 0) {
       get_avalue_explicit(time, lfield);
     } else if (gauge.compare("velocity2") == 0) {
       get_avalue_implicit(time, lfield);
     }
     get_evalue(time, &lfield[3]);
     get_avalue(time, &lfield[6]);
  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_evalue(double time, double* efield) const
{
  if (env_type.compare("pulse") == 0) {
    if (std::abs(time) < 1.0E-10) {
      efield[2] += famp;
      //DEBUG
      printf("### clfield1::get_evalue: %20.10f%20.10f\n", time, efield[2]);
      //DEBUG
    }
  } else if (time > tleft && time < tright || env_type.compare("gauss") == 0) {
     if (gauge.compare("length1") == 0) {
       get_evalue_explicit(time, efield);
     } else if (gauge.compare("length2") == 0) {
       get_evalue_implicit(time, efield);
     } else if (gauge.compare("velocity1") == 0) {
       get_evalue_implicit(time, efield);
     } else if (gauge.compare("velocity2") == 0) {
       get_evalue_explicit(time, efield);
     }
  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_avalue(double time, double* afield) const
{
  if (time > tleft && time < tright || env_type.compare("gauss") == 0) {
     if (gauge.compare("length1") == 0) {
       get_avalue_implicit(time, afield);
     } else if (gauge.compare("length2") == 0) {
       get_avalue_explicit(time, afield);
     } else if (gauge.compare("velocity1") == 0) {
       get_avalue_explicit(time, afield);
     } else if (gauge.compare("velocity2") == 0) {
       get_avalue_implicit(time, afield);
     }
  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_der(double time, double* dfield) const
{
  if (env_type.compare("pulse") == 0) {
    if (std::abs(time) < 1.0E-10) {
      //DEBUG
      printf("### clfield1::get_der: WARNING: derivative of delta-function?\n");
      //DEBUG
    }
  } else if (time > tleft && time < tright || env_type.compare("gauss") == 0) {
    if (gauge.compare("length1") == 0) {
      get_eder_explicit(time, dfield);
    } else if (gauge.compare("length2") == 0) {
      get_eder_implicit(time, dfield);
    } else if (gauge.compare("velocity1") == 0) {
      get_ader_explicit(time, dfield);
    } else if (gauge.compare("velocity2") == 0) {
      get_ader_implicit(time, dfield);
    }
  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_value(double t1, double t2, double c1, double c2, double* lfield) const
{
  double lfield1[9];
  double lfield2[9];
  for (long i = 0; i < 9; i++) {
    lfield1[i] = ZERO;
    lfield2[i] = ZERO;
  }
  get_value(t1, lfield1);
  get_value(t2, lfield2);
  for (long i = 0; i < 9; i++) {
    lfield[i] += c1 * lfield1[i] + c2 * lfield2[i];
  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_ader_explicit(double time, double* dfield) const
{
  double tfield[9];
  for (long i = 0; i < 9; i++) {
    tfield[i] = ZERO;
  }
  get_evalue_implicit(time, tfield);
  dfield[0] += -tfield[0];
  dfield[1] += -tfield[1];
  dfield[2] += -tfield[2];
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_ader_implicit(double time, double* dfield) const
{
  double tfield[9];
  for (long i = 0; i < 9; i++) {
    tfield[i] = ZERO;
  }
  get_evalue_explicit(time, tfield);
  dfield[0] += -tfield[0];
  dfield[1] += -tfield[1];
  dfield[2] += -tfield[2];
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_evalue_explicit(double time, double* lfield) const
{
  double env;
  double treltl = time - tleft;
  double cyc = treltl / period;

  if (env_type.compare("pulse") == 0) {
    std::cout << "clfield1::get_evalue_explicit: DON't call me with env_type = pulse." << std::endl;
    abort();
  } else if (env_type.compare("flat") == 0) {
    lfield[2] += famp * sin(freq * treltl + cep);
  } else if (env_type.compare("trape") == 0) {
    if (cyc < cyc1) {
  	env = cyc / cyc1;
    } else if (cyc < cyc2) {
  	env = ONE;
    } else if (cyc < cyc3) {
  	env = (cyc - cyc3) / (cyc2 - cyc3);
    }
    lfield[2] += env * famp * sin(freq * treltl + cep);
  } else if (env_type.compare("sin2") == 0) {
    env = sin(PI * treltl / (numcyc * period));
    env = env * env;
    lfield[2] += env * famp * sin(freq * treltl + cep);
  } else if (env_type.compare("sin2flat") == 0) {
    double ttime;
    double sin2cyc;
    if (cyc < cyc1) {
  	ttime = treltl;
  	sin2cyc = cyc1 * TWO;
  	env = sin(PI * ttime / (sin2cyc * period));
  	env = env * env;
    } else if (cyc < cyc2) {
  	env = ONE;
    } else if (cyc < cyc3) {
  	sin2cyc = (cyc3 - cyc2) * TWO;
  	ttime = treltl - cyc2 * period + (cyc3 - cyc2) * period;
  	env = sin(PI * ttime / (sin2cyc * period));
  	env = env * env;
    }
    lfield[2] += env * famp * sin(freq * treltl + cep);
  } else if (env_type.compare("sinc") == 0) {
    double trelt0 = time - tcenter;
    if (std::abs(trelt0) > 1.E-10) {
      env = cos(dfreq * trelt0) / (freq * trelt0);
      lfield[2] += env * famp * sin(freq * trelt0);
    } else {
      lfield[2] += famp * cos(dfreq * trelt0);
    }
  } else if (env_type.compare("gauss") == 0) {
    double trelt0 = time - tcenter;
    env = exp(-trelt0*trelt0/(2*sigma*sigma));
    lfield[2] += env * famp * sin(freq * trelt0 + cep);
  } else {
    std::cout << "clfield1::get_evalue_explicit: bad env_type." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_eder_explicit(double time, double* dfield) const
{
  double env, denv, fac;
  double treltl = time - tleft;
  double cyc = treltl / period;

  if (env_type.compare("pulse") == 0) {
    std::cout << "clfield1::get_eder_explicit: DON't call me with env_type = pulse." << std::endl;
    abort();
    dfield[2] += ZERO;
  } else if (env_type.compare("flat") == 0) {
    dfield[2] += famp * freq * cos(freq * treltl + cep);
  } else if (env_type.compare("sin2") == 0) {
    fac = PI / (numcyc * period);
    env = sin(fac*treltl);
    env = env * env;
    denv = TWO * fac * sin(fac*treltl) * cos(fac*treltl);
    dfield[2] += famp * (freq * env * cos(freq * treltl + cep) + denv * sin(freq * treltl + cep));
  } else {
    double field_now[3];
    double k1, k2, k4;
    double etmp;
    double dstep = TWO*PI/(1239.84190/(800.0*27.2113845))/100000.0; // 100000 cycles for 800 nm pulse
  
    etmp = ZERO;
    field_now[2] = ZERO;
    get_evalue_explicit(time + 3*dstep, field_now);
    etmp += field_now[2] / 60.0;
  
    field_now[2] = ZERO;
    get_evalue_explicit(time + 2*dstep, field_now);
    etmp -= field_now[2] / 20.0 * 3.0;
  
    field_now[2] = ZERO;
    get_evalue_explicit(time + dstep, field_now);
    etmp += field_now[2] / 4.0 * 3.0;
  
    field_now[2] = ZERO;
    get_evalue_explicit(time - dstep, field_now);
    etmp -= field_now[2] / 4.0 * 3.0;
  
    field_now[2] = ZERO;
    get_evalue_explicit(time - 2*dstep, field_now);
    etmp += field_now[2] / 20.0 * 3.0;
  
    field_now[2] = ZERO;
    get_evalue_explicit(time - 3*dstep, field_now);
    etmp -= field_now[2] / 60.0;
    dfield[2] += etmp / dstep;
  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_eder_implicit(double time, double* dfield) const
{
  double env, denv, fac;
  double treltl = time - tleft;
  double cyc = treltl / period;

  if (env_type.compare("flat") == 0) {
    dfield[2] += famp * cos(freq*treltl + cep);
  } else if (env_type.compare("gauss") == 0) {
    double trelt0 = time - tcenter;
    env = exp(-trelt0*trelt0/(2*sigma*sigma));
    double sig2 = sigma * sigma;
    double sig4 = sig2 * sig2;
    dfield[2] += env * famp * (cos(freq*treltl + cep) * TWO * treltl / sig2
			     + sin(freq*treltl + cep) * (sig4*freq*freq + sig2 - treltl*treltl) / (sig4 * freq));
  } else {
    double field_now[3];
    double k1, k2, k4;
    double etmp;
    double dstep = TWO*PI/(1239.84190/(800.0*27.2113845))/100000.0; // 100000 cycles for 800 nm pulse
  
    etmp = ZERO;
    field_now[2] = ZERO;
    get_evalue_implicit(time + 3*dstep, field_now);
    etmp += field_now[2] / 60.0;
  
    field_now[2] = ZERO;
    get_evalue_implicit(time + 2*dstep, field_now);
    etmp -= field_now[2] / 20.0 * 3.0;
  
    field_now[2] = ZERO;
    get_evalue_implicit(time + dstep, field_now);
    etmp += field_now[2] / 4.0 * 3.0;
  
    field_now[2] = ZERO;
    get_evalue_implicit(time - dstep, field_now);
    etmp -= field_now[2] / 4.0 * 3.0;
  
    field_now[2] = ZERO;
    get_evalue_implicit(time - 2*dstep, field_now);
    etmp += field_now[2] / 20.0 * 3.0;
  
    field_now[2] = ZERO;
    get_evalue_implicit(time - 3*dstep, field_now);
    etmp -= field_now[2] / 60.0;
    dfield[2] += etmp / dstep;
  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_evalue_implicit(double time, double* lfield) const
{
  if (env_type.compare("zero") != 0) {
    get_evalue_numder(time, lfield);
  }
//  double env;
//  double treltl = time - tleft;
//  double cyc = treltl / period;
//
//  if (env_type.compare("flat") == 0) {
//    lfield[2] += famp * cos(freq*treltl + cep);
//  } else if (env_type.compare("trape") == 0) {
//    double phase;
//    double slope;
//    phase = freq * treltl + cep;
//    if (cyc < cyc1) {
//      slope = ONE / (cyc1 * period);
//      lfield[2] += famp * slope * (cos(phase) * treltl + sin(phase) / freq);
//    } else if (cyc < cyc2) {
//      lfield[2] += famp * cos(phase);
//    } else if (cyc < cyc3) {
//      slope = ONE / ((cyc2 - cyc3) * period);
//      lfield[2] += famp * slope * (cos(phase) * (treltl - cyc3 * period) + sin(phase) / freq);
//    }
//  } else if (env_type.compare("sin2") == 0) {
//    double fac0, facp, facm;
//    double wpls = (ONE + ONE / numcyc) * freq;
//    double wmns = (ONE - ONE / numcyc) * freq;
//    fac0 = cos(freq * treltl + cep);
//    facp = cos(wpls * treltl + cep) * wpls / freq * HALF;
//    facm = cos(wmns * treltl + cep) * wmns / freq * HALF;
//    lfield[2] += - HALF * famp * (fac0 - facp - facm);
//  } else if (env_type.compare("gauss") == 0) {
//    double trelt0 = time - tcenter;
//    env = exp(-trelt0*trelt0/(2*sigma*sigma));
//    lfield[2] += env * famp * (sin(freq*trelt0 + cep) * trelt0/(freq*sigma*sigma)
//  			       - cos(freq*trelt0 + cep));
//  } else {
//    std::cout << "clfield1::get_evalue_implicit: bad env_type." << std::endl;
//    abort();
//  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_avalue_explicit(double time, double* lfield) const
{
  double env;
  double treltl = time - tleft;
  double cyc = treltl / period;

  if (env_type.compare("flat") == 0) {
    lfield[2] += famp / freq * sin(freq * treltl + cep);
  } else if (env_type.compare("trape") == 0) {
    if (cyc < cyc1) {
      env = cyc / cyc1;
    } else if (cyc < cyc2) {
      env = ONE;
    } else if (cyc < cyc3) {
      env = (cyc - cyc3) / (cyc2 - cyc3);
    }
    lfield[2] += env * famp / freq * sin(freq * treltl + cep);
  } else if (env_type.compare("sin2") == 0) {
    env = sin(PI * treltl / (numcyc * period));
    env = env * env;
    lfield[2] += env * famp / freq * sin(freq * treltl + cep);
  } else if (env_type.compare("sin2flat") == 0) {
    double ttime;
    double sin2cyc;
    if (cyc < cyc1) {
      ttime = treltl;
      sin2cyc = cyc1 * TWO;
      env = sin(PI * ttime / (sin2cyc * period));
      env = env * env;
    } else if (cyc < cyc2) {
      env = ONE;
    } else if (cyc < cyc3) {
      sin2cyc = (cyc3 - cyc2) * TWO;
      ttime = treltl - cyc2 * period + (cyc3 - cyc2) * period;
      env = sin(PI * ttime / (sin2cyc * period));
      env = env * env;
    }
    lfield[2] += env * famp / freq * sin(freq * treltl + cep);
  } else if (env_type.compare("gauss") == 0) {
    double trelt0 = time - tcenter;
    env = exp(-trelt0*trelt0/(2*sigma*sigma));
    lfield[2] += env * famp / freq * sin(freq * trelt0 + cep);
//iwatsu version if (field_num == 1) {
//iwatsu version   lfield[2] += env * famp / freq * sin(freq * treltl + cep);
//iwatsu version } else {
//iwatsu version   //    lfield[2] += env * famp / freq * sin(freq * time + cep);
//iwatsu version   //    lfield[2] += env * famp / freq * sin(freq * (time+1000.0/24.188843265) + cep);
//iwatsu version   //    lfield[2] += env * famp / freq * sin(freq * (time+1000.0/24.188843265+231.74006) + cep);
//iwatsu version   lfield[2] += env * famp / freq * sin(freq * (time-1000.0/24.188843265) + cep);
//iwatsu version }
  } else {
    std::cout << "clfield1::get_avalue_explicit: bad env_type." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clfield1::get_avalue_implicit(double time, double* lfield) const
{
  if (env_type.compare("zero") != 0) {
    get_avalue_numint(time, lfield);
  }
//  double env;
//  double treltl = time - tleft;
//  double cyc = treltl / period;
//  if (env_type.compare("zero") == 0) {
//  } else if (env_type.compare("trape") == 0) {
//    double slope, phase, p1, p2;
//    phase = freq * treltl + cep;
//    p1 = freq * cyc1 * period + cep;
//    p2 = freq * cyc2 * period + cep;
//    if (cyc < cyc1) {
//      slope = ONE / (cyc1 * period);
//      lfield[2] += famp / freq * slope * (cos(phase) * treltl - (sin(phase) - sin(cep)) / freq);
//    } else if (cyc < cyc2) {
//      lfield[2] += famp / freq * (cos(phase) - (sin(p1) - sin(cep)) / freq);
//    } else if (cyc < cyc3) {
//      slope = ONE / ((cyc2 - cyc3) * period);
//      lfield[2] += famp / freq * (
//				 slope * (cos(phase) * (treltl - cyc3 * period) - (sin(phase) - sin(p2)) / freq)
//				 - (sin(p1) - sin(cep)) / freq
//				);
//    }
//  } else if (env_type.compare("sin2") == 0) {
//    double fac0, facp, facm;
//    double wpls = (ONE + ONE / numcyc) * freq;
//    double wmns = (ONE - ONE / numcyc) * freq; 
//    fac0 = (cos(freq * treltl + cep) - cos(cep)) / freq;
//    facp = (cos(wpls * treltl + cep) - cos(cep)) / wpls * HALF;
//    facm = (cos(wmns * treltl + cep) - cos(cep)) / wmns * HALF;
//    lfield[2] += HALF * famp * (fac0 - facp - facm);
//  } else if (env_type.compare("gauss") == 0) {
//  } else {
//    std::cout << "clfield1::get_avalue_implicit: bad env_type." << std::endl;
//    abort();
//  }
}
//////////////////////////////////////////////////////////////////////////
void clfield1::get_evalue_numder(double time, double* lfield) const
{
  double field_now[3];
  double k1, k2, k4;
  double etmp;
  double dstep = TWO*PI/(1239.84190/(800.0*27.2113845))/100000.0; // 100000 cycles for 800 nm pulse

  etmp = ZERO;
  field_now[2] = ZERO;
  get_avalue_explicit(time + 3*dstep, field_now);
  etmp += field_now[2] / 60.0;

  field_now[2] = ZERO;
  get_avalue_explicit(time + 2*dstep, field_now);
  etmp -= field_now[2] / 20.0 * 3.0;

  field_now[2] = ZERO;
  get_avalue_explicit(time + dstep, field_now);
  etmp += field_now[2] / 4.0 * 3.0;

  field_now[2] = ZERO;
  get_avalue_explicit(time - dstep, field_now);
  etmp -= field_now[2] / 4.0 * 3.0;

  field_now[2] = ZERO;
  get_avalue_explicit(time - 2*dstep, field_now);
  etmp += field_now[2] / 20.0 * 3.0;

  field_now[2] = ZERO;
  get_avalue_explicit(time - 3*dstep, field_now);
  etmp -= field_now[2] / 60.0;

  lfield[2] += -etmp / dstep;
}
//////////////////////////////////////////////////////////////////////////
void clfield1::get_avalue_numint(double time, double* lfield) const
{
  double field_now[3];
  double k1, k2, k4;
  double time0, atmp;
  double dstep = time - time_prev[field_num];

  double dsmall = TWO*PI/(1239.84190/(800.0*27.2113845))/100000.0; // 100000 cycles for 800 nm pulse
  long numddt = 1;
  while (dstep / numddt > dsmall) {
    numddt *= 10;
  }
  double ddstep = dstep / numddt;

  atmp = aval_prev[field_num];
  for (long istep = 0; istep < numddt; istep ++) {
    time0 = time_prev[field_num] + istep * ddstep;
    field_now[2] = ZERO;
    get_evalue_explicit(time0, field_now);
    k1 = - ddstep * field_now[2];
    field_now[2] = ZERO;
    get_evalue_explicit(time0 + ddstep * HALF, field_now);
    k2 = - ddstep * field_now[2];
    field_now[2] = ZERO;
    get_evalue_explicit(time0 + ddstep, field_now);
    k4 = - ddstep * field_now[2];
    atmp += (k1 + FOUR * k2 + k4) / SIX;
  }
  lfield[2] += atmp;
  //debug
  //  printf("clfield1: %10ld%20.10e\n", numddt, ddstep);
  //  printf("time: %20.10e%20.10e\n", time_prev[field_num], time);
  //  printf("aval: %20.10e%20.10e\n", aval_prev[field_num], atmp);
  //debug

  aval_prev[field_num] = atmp;
  time_prev[field_num] = time;
}
//////////////////////////////////////////////////////////////////////////
