////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
// Sato_tSURFF
//// Orimo_ECS
//#include "surff.hpp"
//// Orimo_ECS
// Sato_tSURFF
////////////////////////////////////////////////////////////////////////
clcrk::clcrk()
{
  std::cout << "clcrk" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clcrk::~clcrk()
{
  std::cout << "~clcrk" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clcrk::clcrk(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
	     const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void clcrk::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		const clfield& Field, const clhprod& HPW)
{
  std::cout << "clcrk::gen" << std::endl;
  clh2prop::gen_basic(MPIP, IO, Bas, Field, HPW);

  IO.read_info("rk_order", LFOUR, rk_order);
  IO.read_info("rk_formula", "standard", rk_formula);

  rk_max_stage = 10;
  rk_nodes.resize(rk_max_stage);
  rk_wghts.resize(2);
  rk_wghts[0].resize(rk_max_stage);
  rk_wghts[1].resize(rk_max_stage);
  rk_coeff.resize(rk_max_stage);
  for (int is = 0; is < rk_max_stage; is ++) {
    rk_coeff[is].resize(rk_max_stage);
  }
  set_coeff();

  Wfn0.gen(MPIP, IO, Bas);
  tWfn.gen(MPIP, IO, Bas);
  kWfn.resize(rk_nstage);
  for (int is = 0; is < rk_nstage; is ++) {
    kWfn[is].gen(MPIP, IO, Bas);
  }
}
////////////////////////////////////////////////////////////////////////
void clcrk::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field, 
		 clhprod& HProd, clwfn& Wfn)
{
  double ttime;
  double lfield[9];
  HProd.copy(Proc, Bas, Wfn, Wfn0);

  dcomplex fac, dtwgt;
  for (int is = 0; is < rk_nstage; is ++) {
    dtwgt = Field.dtime * rk_wghts[rk_level][is];
    ttime = Field.time + Field.dtime * rk_nodes[is];

    HProd.copy(Proc, Bas, Wfn0, tWfn);
    for (int js = 0; js < is; js ++) {
      fac = Field.dtime * rk_coeff[is][js];
      HProd.axpy(Proc, Bas, fac, kWfn[js], tWfn);
    }

    Field.get_value(ttime, lfield);
    HProd.htot(Proc, Bas, ONE, lfield, tWfn, kWfn[is]);
// Sato_tSURFF
//// Orimo_ECS
//    //=============================================
//    if( (Wfn.surff->nprint > 0) && (is == 0) ){
//      //Wfn.surff->rec_v2xmat(ONE, HProd.xmat);
//      Wfn.surff->rec_v2xmat(ONE);
//    }
//    //=============================================
//// Orimo_ECS
// Sato_tSURFF
    HProd.axpy(Proc, Bas, dtwgt, kWfn[is], Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void clcrk::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field, 
		 double time, double dtime, clhprod& HProd, clwfn& Wfn)
{
  double ttime;
  double lfield[9];
  HProd.copy(Proc, Bas, Wfn, Wfn0);

  dcomplex fac, dtwgt;
  for (int is = 0; is < rk_nstage; is ++) {
    dtwgt = dtime * rk_wghts[rk_level][is];
    ttime = time + dtime * rk_nodes[is];

    HProd.copy(Proc, Bas, Wfn0, tWfn);
    for (int js = 0; js < is; js ++) {
      fac = dtime * rk_coeff[is][js];
      HProd.axpy(Proc, Bas, fac, kWfn[js], tWfn);
    }

    Field.get_value(ttime, lfield);
    HProd.htot(Proc, Bas, ONE, lfield, tWfn, kWfn[is]);
// Sato_tSURFF
//// Orimo_ECS
//    //=============================================
//    if( (Wfn.surff->nprint > 0) && (is == 0) ){
//      //Wfn.surff->rec_v2xmat(ONE, HProd.xmat);
//      Wfn.surff->rec_v2xmat(ONE);
//    }
//    //=============================================
//// Orimo_ECS
// Sato_tSURFF
    HProd.axpy(Proc, Bas, dtwgt, kWfn[is], Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff()
{
  if (rk_order == 1) {
    // 1st order: Euler
    rk_level = 0;
    rk_nstage = 1; 
    set_coeff_01();
  } else if (rk_order == 2) {
    // 2nd order: midpoint (better than Henu's method)
    rk_level = 0;
    rk_nstage = 2;
    set_coeff_02();
  } else if (rk_order == 3) {
    // 3rd order: Kutta
    rk_level = 0;
    rk_nstage = 3;
    set_coeff_03();
  } else if (rk_order == 4) {
    if (rk_formula.compare("standard") == 0) {
      // 4th order: standard
      rk_level = 0;
      rk_nstage = 4;
      set_coeff_04();
    } else if (rk_formula.compare("kutta") == 0) {
      // 4th order: Kutta's 3/8 rule
      rk_level = 0;
      rk_nstage = 4;
      set_coeff_04_kutta();
    } else {
      std::cout << "clcrk::gen. bad rk_formula." << std::endl;
    }
  } else if (rk_order == 5) {
    if (rk_formula.compare("fehlberg") == 0) {
      // 5th order: taken from Fehlberg's RK5(4)
      rk_level = 0;
      rk_nstage = 6;
      set_coeff_54_fehlberg();
    } else if (rk_formula.compare("dormand") == 0) {
      // 5th order: taken from Dormand-Prince RK5(4)
      rk_level = 0;
      rk_nstage = 6;
      set_coeff_54_dormand();
    } else {
      std::cout << "clcrk::gen. bad rk_formula." << std::endl;
    }
  } else if (rk_order == 6) {
    if (rk_formula.compare("prince") == 0) {
      // 6th order: taken from Prince's RK6(5)
      rk_level = 0;
      rk_nstage = 8;
      set_coeff_65_prince();
    } else if (rk_formula.compare("calvo") == 0) {
      // 6th order: taken from Calvo's RK6(5)
      rk_level = 0;
      rk_nstage = 8;
      set_coeff_65_calvo();
    } else {
      std::cout << "clcrk::gen. bad rk_formula." << std::endl;
    }
  } else {
    std::cout << "clcrk::gen. bad rk_order." << std::endl;
    abort();
  }

  std::cout << "# rk_level = " << rk_level << std::endl;
  std::cout << "# rk_nstage = " << rk_nstage << std::endl;
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff_01()
{
  // First-order Euler method
  rk_nodes[0] = ZERO;
  rk_wghts[0][0] = ONE;
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff_02()
{
  // Second-order midpoint method
  rk_nodes[0] = ZERO;
  rk_nodes[1] = HALF;
  rk_wghts[0][0] = ZERO;
  rk_wghts[0][1] = ONE;
  rk_coeff[1][0] = HALF;
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff_03()
{
  // Third-order method by Kutta
  rk_nodes[0] = ZERO;
  rk_nodes[1] = HALF;
  rk_nodes[2] = ONE;
  rk_wghts[0][0] = 1.0 / 6.0;
  rk_wghts[0][1] = 2.0 / 3.0;
  rk_wghts[0][2] = 1.0 / 6.0;
  rk_coeff[1][0] = HALF;
  rk_coeff[2][0] = -ONE;  rk_coeff[2][1] = TWO;
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff_04()
{
  // Forth-order standard RK4
  rk_nodes[0] = ZERO;
  rk_nodes[1] = HALF;
  rk_nodes[2] = HALF;
  rk_nodes[3] = ONE;
  rk_wghts[0][0] = ONE / SIX;
  rk_wghts[0][1] = ONE / THREE;
  rk_wghts[0][2] = ONE / THREE;
  rk_wghts[0][3] = ONE / SIX;
  rk_coeff[1][0] = HALF;
  rk_coeff[2][0] = ZERO; rk_coeff[2][1] = HALF;
  rk_coeff[3][0] = ZERO; rk_coeff[3][1] = ZERO; rk_coeff[3][2] = ONE;
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff_04_kutta()
{
  // Forth-order 3/8-rule by Kutta
  rk_nodes[0] = ZERO;
  rk_nodes[1] = 1.0 / 3.0;
  rk_nodes[2] = 2.0 / 3.0;
  rk_nodes[3] = ONE;
  rk_wghts[0][0] = 1.0 / 8.0;
  rk_wghts[0][1] = 3.0 / 8.0;
  rk_wghts[0][2] = 3.0 / 8.0;
  rk_wghts[0][3] = 1.0 / 8.0;
  rk_coeff[1][0] = 1.0 / 3.0;
  rk_coeff[2][0] =-1.0 / 3.0; rk_coeff[2][1] = ONE;
  rk_coeff[3][0] = ONE;       rk_coeff[3][1] =-ONE; rk_coeff[3][2] = ONE;
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff_54_fehlberg()
{
  rk_nodes[0] = ZERO;
  rk_nodes[1] = 1.0 / 4.0;
  rk_nodes[2] = 3.0 / 8.0;
  rk_nodes[3] = 12.0 / 13.0;
  rk_nodes[4] = ONE;
  rk_nodes[5] = 1.0 / 2.0;

  // 5th order weights
  rk_wghts[0][0] = 16.0 / 135.0;
  rk_wghts[0][1] = ZERO;
  rk_wghts[0][2] = 6656.0 / 12825.0;
  rk_wghts[0][3] = 28561.0 / 56430.0;
  rk_wghts[0][4] =-9.0 / 50.0;
  rk_wghts[0][5] = 2.0 / 55.0;

  // 4th order weights
  rk_wghts[1][0] = 25.0 / 216.0;
  rk_wghts[1][1] = ZERO;
  rk_wghts[1][2] = 1408.0 / 2565.0;
  rk_wghts[1][3] = 2197.0 / 4104.0;
  rk_wghts[1][4] =-1.0 / 5.0;
  rk_wghts[1][5] = ZERO;

  rk_coeff[1][0] = 1.0 / 4.0; // 1st
  rk_coeff[2][0] = 3.0 / 32.0; // 2nd
  rk_coeff[2][1] = 9.0 / 32.0;
  rk_coeff[3][0] = 1932.0 / 2197.0; // 3rd
  rk_coeff[3][1] =-7200.0 / 2197.0;
  rk_coeff[3][2] = 7296.0 / 2197.0;
  rk_coeff[4][0] = 439.0 / 216.0; // 4th
  rk_coeff[4][1] =-8.0;
  rk_coeff[4][2] = 3680.0 / 513.0;
  rk_coeff[4][3] =-845.0 / 4104.0;
  rk_coeff[5][0] =-8.0 / 27.0; // 5th
  rk_coeff[5][1] = 2.0;
  rk_coeff[5][2] =-3544.0 / 2565.0;
  rk_coeff[5][3] = 1859.0 / 4104.0;
  rk_coeff[5][4] =-11.0 / 40.0;
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff_54_dormand()
{
  rk_nodes[0] = ZERO;
  rk_nodes[1] = 1.0 / 5.0;
  rk_nodes[2] = 3.0 / 10.0;
  rk_nodes[3] = 4.0 / 5.0;
  rk_nodes[4] = 8.0 / 9.0;
  rk_nodes[5] = ONE;
  rk_nodes[6] = ONE;

  // 5th order weights
  rk_wghts[0][0] = 35.0 / 384.0;   
  rk_wghts[0][1] = ZERO;
  rk_wghts[0][2] = 500.0 / 1113.0;
  rk_wghts[0][3] = 125.0 / 192.0;
  rk_wghts[0][4] =-2187.0 / 6784.0;
  rk_wghts[0][5] = 11.0 / 84.0;
  rk_wghts[0][6] = ZERO;

  // 4th order weights
  rk_wghts[1][0] = 5179.0 / 57600.0;
  rk_wghts[1][1] = ZERO;
  rk_wghts[1][2] = 7571.0 / 16695.0;
  rk_wghts[1][3] = 393.0 / 640.0;
  rk_wghts[1][4] =-92097.0 / 339200.0;
  rk_wghts[1][5] = 187.0 / 2100.0;
  rk_wghts[1][6] = 1.0 / 40.0;

  // taken from wiki (which is different from ppt of Douglas Wilhelm Harder...)
  rk_coeff[1][0] = 1.0 / 5.0; // 1st
  rk_coeff[2][0] = 3.0 / 40.0; // 2nd
  rk_coeff[2][1] = 9.0 / 40.0;
  rk_coeff[3][0] = 44.0 / 45.0; // 3rd
  rk_coeff[3][1] =-56.0 / 15.0;
  rk_coeff[3][2] = 32.0 / 9.0;
  rk_coeff[4][0] = 19372.0 / 6561.0; // 4th
  rk_coeff[4][1] =-25360.0 / 2187.0;
  rk_coeff[4][2] = 64448.0 / 6561.0;
  rk_coeff[4][3] =-212.0 / 729.0;
  rk_coeff[5][0] = 9017.0 / 3168.0; // 5th
  rk_coeff[5][1] =-355.0 / 33.0;
  rk_coeff[5][2] = 46732.0 / 5247.0;
  rk_coeff[5][3] = 49.0 / 176.0;
  rk_coeff[5][4] =-5103.0 / 18656.0;
  rk_coeff[6][0] = 35.0 / 384.0;    // 6th
  rk_coeff[6][1] = ZERO;
  rk_coeff[6][2] = 500.0 / 1113.0;
  rk_coeff[6][3] = 125.0 / 192.0;
  rk_coeff[6][4] =-2187.0 / 6784.0;
  rk_coeff[6][5] = 11.0 / 84.0;
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff_65_prince()
// P. J. Prince and J. R. Dormand,
// "High order embedded Runge-Kutta formulae",
// J. Comput. Appl. Math, 7, 67 (1981).
{
  rk_nodes[0] = ZERO;
  rk_nodes[1] = 1.0 / 10.0;
  rk_nodes[2] = 2.0 / 9.0;
  rk_nodes[3] = 3.0 / 7.0;
  rk_nodes[4] = 3.0 / 5.0;
  rk_nodes[5] = 4.0 / 5.0;
  rk_nodes[6] = ONE;
  rk_nodes[7] = ONE;

  // 6th order weights
  rk_wghts[0][0] = 61.0 / 864.0;
  rk_wghts[0][1] = ZERO;
  rk_wghts[0][2] = 98415.0 / 321776.0;
  rk_wghts[0][3] = 16807.0 / 146016.0;
  rk_wghts[0][4] = 1375.0 / 7344.0;
  rk_wghts[0][5] = 1375.0 / 5408.0;
  rk_wghts[0][6] =-37.0 / 1120.0;
  rk_wghts[0][7] = 1.0 / 10.0;

  // 5th order weights
  rk_wghts[0][0] = 821.0 / 10800.0;
  rk_wghts[0][1] = ZERO;
  rk_wghts[0][2] = 19683.0 / 71825.0;
  rk_wghts[0][3] = 175273.0 / 912600.0;
  rk_wghts[0][4] = 395.0 / 3672.0;
  rk_wghts[0][5] = 785.0 / 2704.0;
  rk_wghts[0][6] = 3.0 / 50.0;
  rk_wghts[0][7] = ZERO;

  rk_coeff[1][0] = 1.0 / 10.0; // row 1
  rk_coeff[2][0] =-2.0 / 81.0; // row 2
  rk_coeff[2][1] = 20.0 / 81.0;
  rk_coeff[3][0] = 615.0 / 1372.0; // row 3
  rk_coeff[3][1] =-270.0 / 343.0;
  rk_coeff[3][2] = 1053.0 / 1372.0;
  rk_coeff[4][0] = 3243.0 / 5500.0;   // row 4
  rk_coeff[4][1] =-54.0 / 55.0;
  rk_coeff[4][2] = 50949.0 / 71500.0;
  rk_coeff[4][3] = 4998.0 / 17875.0;
  rk_coeff[5][0] =-26492.0 / 37125.0;    // row 5
  rk_coeff[5][1] = 72.0 / 55.0;
  rk_coeff[5][2] = 2808.0 / 23375.0;
  rk_coeff[5][3] =-24206.0 / 37125.0;
  rk_coeff[5][4] = 338.0 / 459.0;
  rk_coeff[6][0] = 5561.0 / 2376.0; // row 6
  rk_coeff[6][1] =-35.0 / 11.0;
  rk_coeff[6][2] =-24117.0 / 31603.0;
  rk_coeff[6][3] = 899983.0 / 200772.0;
  rk_coeff[6][4] =-5225.0 / 1836.0;
  rk_coeff[6][5] = 3925.0 / 4056.0;
  rk_coeff[7][0] = 465467.0 / 266112.0; // row 7
  rk_coeff[7][1] =-2945.0 / 1232.0;
  rk_coeff[7][2] =-5610201.0 / 14158144.0;
  rk_coeff[7][3] = 10513573.0 / 3212352.0;
  rk_coeff[7][4] =-424325.0 / 205632.0;
  rk_coeff[7][5] = 376225.0 / 454272.0;
  rk_coeff[7][6] = ZERO;
}
////////////////////////////////////////////////////////////////////////
void clcrk::set_coeff_65_calvo()
// M. Calvo, J. I. Montijano, and L. Randez, 
// "A NEW EMBEDDED PAIR OF RUNGE-KUTTA FORMULAS OF ORDERS 5 AND 6",
// Comput. Math. Applic. 20, 15 (1990) 
// 06_calvo has FSAL (First Same As Last) property.
{
  rk_nodes[0] = ZERO; //nodes
  rk_nodes[1] = 2.0 / 15.0;
  rk_nodes[2] = 1.0 / 5.0;
  rk_nodes[3] = 3.0 / 10.0;
  rk_nodes[4] = 14.0 / 25.0;
  rk_nodes[5] = 19.0 / 25.0;
  rk_nodes[6] = 35226607.0 / 35688279.0;
  rk_nodes[7] = ONE;
  rk_nodes[8] = ONE;

  // 6th order weights
  rk_wghts[0][0] = 17572349.0 / 289262523.0;
  rk_wghts[0][1] = ZERO;
  rk_wghts[0][2] = 57513011.0 / 201864250.0;
  rk_wghts[0][3] = 15587306.0 / 354501571.0;
  rk_wghts[0][4] = 71783021.0 / 234982865.0;
  rk_wghts[0][5] = 29672000.0 / 180480167.0;
  rk_wghts[0][6] = 65567621.0 / 127060952.0;
  rk_wghts[0][7] =-79074570.0 / 210557597.0;
  rk_wghts[0][8] = ZERO;

  // 5th order weights
  rk_wghts[1][0] = 15231665.0 / 510830334.0;
  rk_wghts[1][1] = ZERO;
  rk_wghts[1][2] = 59452991.0 / 116050448.0;
  rk_wghts[1][3] =-28398517.0 / 122437738.0;
  rk_wghts[1][4] = 56673824.0 / 137010559.0;
  rk_wghts[1][5] = 68003849.0 / 426673583.0;
  rk_wghts[1][6] = 7097631.0 / 37564021.0;
  rk_wghts[1][7] =-71226429.0 / 583093742.0;
  rk_wghts[1][8] = 1.0 / 20.0;

  rk_coeff[1][0] = 2.0 / 15.0; // row 1
  rk_coeff[2][0] = 1.0 / 20.0; // row 2
  rk_coeff[2][1] = 3.0 / 20.0;
  rk_coeff[3][0] = 3.0 / 40.0; // row 3
  rk_coeff[3][1] = ZERO;
  rk_coeff[3][2] = 9.0 / 40.0;
  rk_coeff[4][0] = 86727015.0 / 196851553.0;   // row 4
  rk_coeff[4][1] =-60129073.0 / 52624712.0;
  rk_coeff[4][2] = 957436434.0 / 1378352377.0;
  rk_coeff[4][3] = 83886832.0 / 147842441.0;
  rk_coeff[5][0] =-86860849.0 / 45628967.0;    // row 5
  rk_coeff[5][1] = 111022885.0 / 25716487.0;
  rk_coeff[5][2] = 108046682.0 / 101167669.0;
  rk_coeff[5][3] =-141756746.0 / 36005461.0;
  rk_coeff[5][4] = 73139862.0 / 60170633.0;
  rk_coeff[6][0] = 77759591.0 / 16096467.0; // row 6
  rk_coeff[6][1] =-49252809.0 / 6452555.0;
  rk_coeff[6][2] =-381680111.0 / 51572984.0;
  rk_coeff[6][3] = 879269579.0 / 66788831.0;
  rk_coeff[6][4] =-90453121.0 / 33722162.0;
  rk_coeff[6][5] = 111179552.0 / 157155827.0;
  rk_coeff[7][0] = 237564263.0 / 39280295.0; // row 7
  rk_coeff[7][1] =-100523239.0 / 10677940.0;
  rk_coeff[7][2] =-265574846.0 / 27330247.0;
  rk_coeff[7][3] = 317978411.0 / 18988713.0;
  rk_coeff[7][4] =-124494385.0 / 35453627.0;
  rk_coeff[7][5] = 86822444.0 / 100138635.0;
  rk_coeff[7][6] =-12873523.0 / 724232625.0;
  rk_coeff[8][0] = 17572349.0 / 289262523.0; // row 8
  rk_coeff[8][1] = ZERO;
  rk_coeff[8][2] = 57513011.0 / 201864250.0;
  rk_coeff[8][3] = 15587306.0 / 354501571.0;
  rk_coeff[8][4] = 71783021.0 / 234982865.0;
  rk_coeff[8][5] = 29672000.0 / 180480167.0;
  rk_coeff[8][6] = 65567621.0 / 127060952.0;
  rk_coeff[8][7] =-79074570.0 / 210557597.0;
}
////////////////////////////////////////////////////////////////////////
