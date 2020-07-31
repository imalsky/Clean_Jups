#include <stdio.h>

double Planck(double T, double lambda);

double two_stream(int NLAYER, int kmin, double w0_val, double g0_val, \
                 double *temperature_array, double *tau_array, \
                 double NU, double NU_BIN, double* TMI, double incident_frac)
{

  kmin = 0;
  NLAYER = 102;






  double mu_1 = 1.0;
  double mu_0 = 1.0;  // Needs to be adjusted for the solid angle stuff
  double mu = 1.0;

  // These are indexing values
  int J, L, KINDEX, Z, M;
  int NEW_NLAYER;

  // These are just constants
  double bolz_constant = 1.380649e-23;
  double h_constant    = 6.62607015e-34;

  // These are boundary conditions are values for the flux stuff
  double EMIS = 0.;
  double RSFX = 0.;
  double redistribution_param = 1.0;

  double SFCS;
  double DIRECT;
  double STELLAR_BB;
  double BB_TOP_OF_ATM;
  double BB_BOTTOM_OF_ATM;

  // Frequency Stuff
  double B_SPECTRAL_DENSITY_VAL;

  // How to generate the matrix from the toon paper
  double y1[NLAYER - kmin];
  double y2[NLAYER - kmin];
  double LAMBDAS[NLAYER - kmin];
  double GAMMA[NLAYER - kmin];
  double temp_e_val[NLAYER - kmin];
  double e1[NLAYER - kmin];
  double e2[NLAYER - kmin];
  double e3[NLAYER - kmin];
  double e4[NLAYER - kmin];

  // Matrix Coefficients
  double A[2 * (NLAYER - kmin)];
  double B[2 * (NLAYER - kmin)];
  double D[2 * (NLAYER - kmin)];
  double E[2 * (NLAYER - kmin)];
  
  // Solution Coefficients
  double AS[2 * (NLAYER - kmin)];
  double DS[2 * (NLAYER - kmin)];
  double X[2 * (NLAYER - kmin)];
  double Y[2 * (NLAYER - kmin)];

  // More temperatoray matrix soluton stuff
  double temp_gamma_val[NLAYER - kmin];
  double CP[NLAYER - kmin];
  double CPB[NLAYER - kmin];
  double CM[NLAYER - kmin];
  double CMB[NLAYER - kmin];

  // Planck Function Stuff, and the slope
  double B0[NLAYER - kmin];
  double B1[NLAYER - kmin];
  double FNET[NLAYER - kmin];
  double Bnu, twohnu3_c2, hc_Tkla;
  double temp_val_1, temp_val_2;

  // The Source Function Variables
  double SOURCE_G[NLAYER - kmin];
  double SOURCE_H[NLAYER - kmin];
  double SOURCE_J[NLAYER - kmin];
  double SOURCE_K[NLAYER - kmin];
  double ALPHA_1[NLAYER - kmin];
  double ALPHA_2[NLAYER - kmin];
  double SIGMA_1[NLAYER - kmin];
  double SIGMA_2[NLAYER - kmin];
  double SOURCE_Y1[NLAYER - kmin];
  double SOURCE_Y2[NLAYER - kmin];
  double source_temp[NLAYER - kmin];

  // Upward and downwards intensities
  double INTENSITY_DOWN[NLAYER - kmin];
  double INTENSITY_UP[NLAYER - kmin];

  // Top layer values
  double TWO_STREAM_INTENSITY;
  double SOURCE_INTENSITY;
  double RUNNING_SOURCE;

  // The number of layers
  // Sometimes the inputs are bad and the top
  // Of the layers need to be cut off
  NEW_NLAYER = NLAYER - kmin;

  // Scattering and atmosphere parameters


  double TAUCS[NLAYER - kmin];


  double W0[102] = {7.434569304525805e-002,7.434569304525805e-002,8.489218091995920e-002,\
  8.489218091995920e-002,7.999692917570497e-002,7.999692917570500e-002,\
  7.999700948909487e-002,7.999700948909487e-002,7.999707041313750e-002,\
  7.999707041313750e-002,7.999711662881116e-002,7.999711662881115e-002,\
  7.999715168701992e-002,7.999715168701992e-002,7.999717828141110e-002,\
  7.999717828141109e-002,7.999719845533153e-002,7.999719845533154e-002,\
  7.999721375882182e-002,7.999721375882182e-002,7.999722536771026e-002,\
  7.999722536771026e-002,7.999723417395492e-002,7.999723417395491e-002,\
  7.999724085417591e-002,7.999724085417591e-002,7.999724592164248e-002,\
  7.999724592164248e-002,7.999724976570943e-002,7.999724976570943e-002,\
  7.999725268173263e-002,7.999725268173262e-002,7.999725489376257e-002,\
  7.999725489376255e-002,8.472911245051135e-002,8.472911245051135e-002,\
  9.297652643092529e-002,9.297652643092529e-002,0.104106028376184,\
  0.104106028376184,0.119699521734725,0.119699521734725,0.141720718398443,\
  0.141720718398443,0.175227814761669,0.175227814761669,0.225587543769275,\
  0.225587543769275,0.294888306376572,0.294888306376572,0.388674297591101,\
  0.388674297591101,0.505744799940412,0.505744799940412,0.625316197772204,\
  0.625316197772204,0.724767146028241,0.724767146028240,0.794301005597684,\
  0.794301005597684,0.838774048523858,0.838774048523858,0.867759950554337,\
  0.867759950554337,0.875664898033914,0.875664898033914,0.869361974930026,\
  0.869361974930026,0.847762272151204,0.847762272151204,0.812950359496062,\
  0.812950359496062,0.768185696870480,0.768185696870480,0.719200708350483,\
  0.719200708350483,0.668192427439879,0.668192427439879,0.607363987381007,\
  0.607363987381007,0.532525735908775,0.532525735908775,0.447442135979522,\
  0.447442135979521,0.202768635530217,0.202768635530217,0.105723381685748,\
  0.105723381685748,7.817670187228151e-002,7.817670187228165e-002,5.979540194933848e-002,\
  5.979540194933853e-002,3.390544549138463e-002,3.390544549138463e-002,0.000000000000000e000,\
  0.000000000000000e000,0.000000000000000e000,0.000000000000000e000,0.000000000000000e000,\
  0.000000000000000e000,0.000000000000000e000,0.000000000000000e000};

  double G0[102] = {0.110326882158668,0.110326882158668,0.110326956918972,0.110326956918972,0.110326772109971,\
  0.110326772109971,0.110326873438366,0.110326873438366,0.110326950303938,0.110326950303938,\
  0.110327008612508,0.110327008612508,0.110327052844123,0.110327052844123,0.110327086397258,\
  0.110327086397258,0.110327111849928,0.110327111849928,0.110327131157760,0.110327131157760,\
  0.110327145804253,0.110327145804253,0.110327156914757,0.110327156914757,0.110327165342937,\
  0.110327165342937,0.110327171736366,0.110327171736366,0.110327176586279,0.110327176586279,\
  0.110327180265315,0.110327180265315,0.110327183056148,0.110327183056148,0.112254085593167,\
  0.112254085593167,0.114959453958986,0.114959453958986,0.118358383332087,0.118358383332087,\
  0.122183827464402,0.122183827464402,0.126155103807160,0.126155103807160,0.130313245599292,\
  0.130313245599292,0.133748779520810,0.133748779520810,0.135375962548969,0.135375962548969,\
  0.134183699412485,0.134183699412485,0.129238005738520,0.129238005738520,0.121600030011815,\
  0.121600030011815,0.115855879644724,0.115855879644724,0.122830902173550,0.122830902173550,\
  0.154889553165107,0.154889553165107,0.225947922855962,0.225947922855962,0.315123504542173,\
  0.315123504542173,0.358428805263750,0.358428805263750,0.388659181988818,0.388659181988818,\
  0.408227725579857,0.408227725579857,0.418568543084830,0.418568543084830,0.419002832244987,\
  0.419002832244987,0.410708458415535,0.410708458415535,0.402779488738229,0.402779488738229,\
  0.408189397139313,0.408189397139313,0.421549423462755,0.421549423462755,0.421108712660989,\
  0.421108712660989,0.375638623840812,0.375638623840812,0.379455427694435,0.379455427694435,\
  0.381757089883631,0.381757089883631,0.345678203649185,0.345678203649185,0.0e000,0.0e000,\
  0.0e000,0.0e000,0.0e000,0.0e000,0.0e000,0.0e000};

  // Scattering and atmosphere parameters
  double TEMPS[102] = {928.729979382128,929.031000326966,929.095129933947,929.151000326972,\
  929.235436839383,929.309000326988,929.420690237127,929.518000327016,929.664423952239,\
  929.792000327059,929.984911928262,930.153000327118,930.407358453202,930.629000327199,\
  930.962966537640,931.254000327308,931.693211116800,932.076000327453,932.652494407311,\
  933.155000327645,933.911485328685,934.571000327897,935.561375424113,936.425000328230,\
  937.719143920565,938.848000328666,940.534942801061,942.007000329237,944.197545837884,\
  946.110000329982,948.943590249928,951.419000330952,955.062542388501,958.248000332209,\
  962.901590036104,966.974000333832,972.865950130157,978.028000335912,985.410924864053,\
  991.888000338559,1001.02324580884,1009.05000034190,1020.18213963746,1029.98000034606,\
  1043.32954525886,1055.10000035120,1070.78554001748,1084.64000035743,1102.69580586505,\
  1118.67000036488,1138.95188187402,1156.92000037361,1179.12707580107,1198.82000038366,\
  1222.42652631029,1243.37000039496,1267.65644529934,1289.20000040738,1313.24651042763,\
  1334.56000042068,1357.34443422876,1377.51000043451,1398.03788093586,1416.17000044845,\
  1433.71686738093,1449.18000046204,1463.61115832364,1476.30000047489,1488.28231649120,\
  1498.80000048672,1509.75301617203,1519.36000049752,1531.02500863645,1541.26000050752,\
  1555.12756017287,1567.31000051716,1584.39479232520,1599.43000052697,1620.44646615001,\
  1638.98000053749,1664.52346097883,1687.10000054922,1717.71087855874,1744.83000056266,\
  1781.00513455798,1813.13000057824,1855.28158028012,1892.80000059640,1941.26196916547,\
  1984.49000061750,2039.53989352302,2088.74000064188,2150.59432335360,2205.97000066983,\
  2274.83142709262,2336.57000070162,2412.62212979420,2480.89000073746,2564.34915782821,\
  2639.34000077758,2722.32051132831,2722.32051132831};
  
  double TAULS[102] = {3.785206799298454e-003,3.785206799298454e-003,4.369965972797778e-003,\
  4.369965972797778e-003,2.595438703931929e-003,2.595438703931928e-003,\
  3.421454547729317e-003,3.421454547729317e-003,4.510355499969544e-003,\
  4.510355499969544e-003,5.945806517880663e-003,5.945806517880664e-003,\
  7.838099495109633e-003,7.838099495109633e-003,1.033262746366466e-002,\
  1.033262746366466e-002,1.362105576772974e-002,1.362105576772974e-002,\
  1.795604853882313e-002,1.795604853882313e-002,2.367068197091029e-002,\
  2.367068197091029e-002,3.120403600113891e-002,3.120403600113891e-002,\
  4.113493071541886e-002,4.113493071541886e-002,5.422639959241263e-002,\
  5.422639959241263e-002,7.148431665711115e-002,7.148431665711115e-002,\
  9.423468212109742e-002,9.423468212109744e-002,0.124225504698618,\
  0.124225504698618,0.164346850475370,0.164346850475370,0.217802839656924,\
  0.217802839656924,0.288312940957087,0.288312940957087,0.381698208868834,\
  0.381698208868835,0.504313298985020,0.504313298985020,0.664601927142552,\
  0.664601927142552,0.870448149355705,0.870448149355705,1.13363591367261,\
  1.13363591367261,1.46837177270562,1.46837177270562,1.89879903080063,\
  1.89879903080063,2.46800061035889,2.46800061035889,3.26133011659987,\
  3.26133011659987,4.44504276493956,4.44504276493956,6.24956799366401,\
  6.24956799366401,8.99077758138895,8.99077758138895,11.5088199165136,\
  11.5088199165136,13.7880334308601,13.7880334308601,14.9228878340194,\
  14.9228878340194,15.4293229311154,15.4293229311154,16.0014607231307,\
  16.0014607231307,17.2370839545222,17.2370839545222,19.2839577291430,\
  19.2839577291430,21.6324646844230,21.6324646844230,23.9839555160984,\
  23.9839555160984,26.5834619303397,26.5834619303397,23.1690453010524,\
  23.1690453010523,24.7407793725123,24.7407793725123,31.5058061322877,\
  31.5058061322876,40.5853307621816,40.5853307621816,51.4183088548791,\
  51.4183088548791,65.3288657406461,65.3288657406460,86.1202174848536,\
  86.1202174848534,113.528557025352,113.528557025351,139.338065787836,\
  1.000000000e-015};



  //**********************************************
  //* Data Sanitation (Sorry, kinda gross code)  *
  //**********************************************
  if (TAULS[1] < 1e-8)
  {
    TAULS[1] = TAULS[2];
  }

  if (TAULS[0] < 1e-8)
  {
    TAULS[0] = TAULS[1];
  }

  if (TAULS[NEW_NLAYER-1] < 1e-8)
  {
    TAULS[NEW_NLAYER-1] = TAULS[NEW_NLAYER-2] + abs(TAULS[NEW_NLAYER-2] - TAULS[NEW_NLAYER-3]);
  }

  if (TEMPS[2] < 1.0)
  {
    TEMPS[2] = TEMPS[3] - (TEMPS[4] - TEMPS[3]);
  }

  if (TEMPS[1] < 1.0)
  {
    TEMPS[1] = TEMPS[2] - (TEMPS[3] - TEMPS[2]);
  }

  if (TEMPS[0] < 1.0)
  {
    TEMPS[0] = TEMPS[1] - (TEMPS[2] - TEMPS[1]);
  }
  
  // Calculate the intensity at the top of the atmosphere
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[NEW_NLAYER-1])) - 1.0;
  BB_TOP_OF_ATM = temp_val_1 * (1.0 / temp_val_2);

  // Calculate the intensity at the bottom of the atmosphere
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[NEW_NLAYER-1])) - 1.0;
  BB_BOTTOM_OF_ATM = temp_val_1 * (1.0 / temp_val_2);

  // Calculate the flux at the top of the atmosphere from the star
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * STELLAR_TEMP)) - 1.0;
  STELLAR_BB = temp_val_1 * (1.0 / temp_val_2);

  // Define the energy from other sources (eq 37 and 38, Toon)
  if (NU > 430.0e12)
  {
  	//redistribution param is unused
    //DIRECT = incident_frac * PI * STELLAR_BB * pow(R_STAR / ORB_SEP, 2.0);
    //SFCS = RSFX * DIRECT * mu_0 * exp(-(TAUCS[NEW_NLAYER-1] + TAULS[NEW_NLAYER-1]) / mu_0);

    SFCS = 0.0;
    DIRECT = 0.0;
  }
  else
  {
    SFCS = EMIS * PI * BB_BOTTOM_OF_ATM;
    DIRECT = 0.0;

    SFCS = 0.0;
    DIRECT = 0.0;

  }

  // HERE WE FIND LAYER PROPERTIES FOLLOWING GENERAL SCHEME
  // OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
  // NEEDED FOR MATRIX.

  for(J=0; J<NEW_NLAYER; J++)
  {
    y1[J]    =  (2.0 - (W0[J] * (1.0 + G0[J])));
    y2[J]    =  (W0[J] * (1.0 - G0[J]));

    LAMBDAS[J]    =  sqrt(fabs(pow(y1[J], 2.0) - pow(y2[J], 2.0)));
    GAMMA[J]  =  y2[J] / (y1[J] + LAMBDAS[J]);
    temp_e_val[J]   =  exp(-LAMBDAS[J] * TAULS[J]);

    e1[J]   =  1.0 + GAMMA[J] * temp_e_val[J];  //e1                          
    e2[J]   =  1.0 - GAMMA[J] * temp_e_val[J];  //e2                      
    e3[J]   =  GAMMA[J] + temp_e_val[J];        //e3                        
    e4[J]   =  GAMMA[J] - temp_e_val[J];        //e4
  }


  J = 0;
  for(L=1; L<2*NEW_NLAYER -1; L+=2)
  {
    // HERE ARE THE EVEN MATRIX ELEMENTS
    A[L]   =  e2[J+1] * e1[J]   - e4[J+1] * e3[J];
    B[L]   =  e2[J+1] * e2[J]   - e4[J+1] * e4[J];
    D[L]   =  e1[J+1] * e4[J+1] - e3[J+1] * e2[J+1];
  
    // HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP. 
    A[L+1] =  e2[J]   * e3[J]   - e1[J]   * e4[J];
    B[L+1] =  e1[J+1] * e1[J]   - e3[J+1] * e3[J]; 
    D[L+1] =  e3[J]   * e4[J+1] - e1[J]   * e2[J+1];
    J = J + 1;
   }

  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME
  // NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
  A[0] = 0.0;
  B[0] = e1[0];
  D[0] = -e2[0];

  A[2*NEW_NLAYER-1] = e1[NEW_NLAYER-1] - RSFX * e3[NEW_NLAYER-1];
  B[2*NEW_NLAYER-1] = e2[NEW_NLAYER-1] - RSFX * e4[NEW_NLAYER-1];
  D[2*NEW_NLAYER-1] = 0.0;

  // This is the part of the code that solves for the blackbody stuff
  for(J=0; J<NEW_NLAYER; J++)
  {
    if(0 >= J-1)
      KINDEX = 0;
	  else
      KINDEX = J-1;
    
    temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
    temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[J])) - 1.0;
    B_SPECTRAL_DENSITY_VAL = temp_val_1 * (1.0 / temp_val_2);



    // This needs to be changed!
    B_SPECTRAL_DENSITY_VAL = 5.67E-8 * (pow(TEMPS[J], 4.0)) / PI;

    B0[J] = B_SPECTRAL_DENSITY_VAL;
    B1[J] = (B0[J] - B0[KINDEX]) / TAULS[J];
  }

  // The very top of the atmosphere is isothermal
  B1[0] = 0;

  // This solves for the C values in the toon code
  for(J=0; J<NEW_NLAYER; J++)
  {
    if(0 > J-1)
      KINDEX = 0;
    else
      KINDEX = J-1;

    temp_gamma_val[J]   = 1.0 / (y1[J] + y2[J]);

    CP[J]  = (B0[KINDEX] + B1[J] * temp_gamma_val[J]) * 2.0 * PI * mu_1;
    CPB[J] = CP[J] + B1[J] * TAULS[J] * 2.0 * PI * mu_1;

    CM[J]  = (B0[KINDEX] - B1[J] * temp_gamma_val[J]) * 2.0 * PI * mu_1;
    CMB[J] = CM[J] + B1[J] * TAULS[J] * 2.0 * PI * mu_1;
  }

  J = 0;
  for(L=1; L<2*NEW_NLAYER; L+=2)
  {
    // HERE ARE THE EVEN MATRIX ELEMENTS
    E[L]   = (CP[J+1] - CPB[J]) * e2[J+1] - (CM[J+1] - CMB[J]) * e4[J+1];

    // HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
    E[L+1] = e3[J] * (CP[J+1] - CPB[J]) + e1[J] * (CMB[J] - CM[J+1]);
    J = J + 1;
  }

  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME NO
  // DIFFUSE RADIATION IS INCIDENT AT THE TOP.

  E[0] = -CM[0];
  E[2*NEW_NLAYER-1]  = SFCS + RSFX * CMB[NEW_NLAYER-1] - CPB[NEW_NLAYER-1];
  DS[2*NEW_NLAYER-1] = E[2*NEW_NLAYER-1] / B[2*NEW_NLAYER-1];
  AS[2*NEW_NLAYER-1] = A[2*NEW_NLAYER-1] / B[2*NEW_NLAYER-1];

  //********************************************
  //*     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
  //********************************************

  for(L=2; L<2*NEW_NLAYER+1; L++)
  {
    X[2*NEW_NLAYER-L]  = 1.0 / (B[2*NEW_NLAYER-L] - D[2*NEW_NLAYER-L] * AS[2*NEW_NLAYER-L+1]);
    AS[2*NEW_NLAYER-L] = A[2*NEW_NLAYER-L] * X[2*NEW_NLAYER-L];
    DS[2*NEW_NLAYER-L] = (E[2*NEW_NLAYER-L]-D[2*NEW_NLAYER-L] * DS[1+2*NEW_NLAYER-L]) * X[2*NEW_NLAYER-L];
  }
  
  Y[0] = DS[0];
  for(L=1; L<2*NEW_NLAYER; L++)
  {
    Y[L] = DS[L] - AS[L] * Y[L-1];
  }

  //***************************************************************
  //  CALCULATE LAYER CODFICIENTS, NET FLUX AND MEAN INTENSITY
  //***************************************************************

  // I reverse the list here because that's how it wants it in Eliza's code
  for(J=1; J<NEW_NLAYER+1; J++)
  {
    FNET[NEW_NLAYER-J] = Y[2*J-2] * (e1[J-1]-e3[J-1]) + Y[2*J-2] \
                 * (e2[J-1]-e4[J-1]) + CPB[J-1] - CMB[J-1] - DIRECT;

    TMI[NEW_NLAYER-J] = (1.0 / mu_1) * (Y[2*J-2]*(e1[J-1] + e3[J-1]) + Y[2*J-1] * (e2[J-1]+e4[J-1]) \
             + CPB[J-1] + CMB[J-1]) +  (DIRECT / mu_0);

    TMI[NEW_NLAYER-J] = TMI[NEW_NLAYER-J];

  }

  //********************************************
  //*    Source Function Technique Solution    *
  //********************************************
  for(J=1; J<NEW_NLAYER+1; J++)
  {
    SOURCE_Y1[J-1] = Y[2*J-2];
    SOURCE_Y2[J-1] = Y[2*J-1];
  }
  
  for(J=0; J<NEW_NLAYER; J++)
  {
    SOURCE_G[J] = (SOURCE_Y1[J] + SOURCE_Y2[J]) * (1.0/mu_1 - LAMBDAS[J]);
    SOURCE_H[J] = (SOURCE_Y1[J] - SOURCE_Y2[J]) * GAMMA[J] * (1.0/mu_1 + LAMBDAS[J]);
    SOURCE_J[J] = (SOURCE_Y1[J] + SOURCE_Y2[J]) * GAMMA[J] * (1.0/mu_1 + LAMBDAS[J]);
    SOURCE_K[J] = (SOURCE_Y1[J] - SOURCE_Y2[J]) * GAMMA[J] * (1.0/mu_1 - LAMBDAS[J]);

    source_temp[J] = (1.0 / (y1[J] + y2[J])) - mu_1;
    ALPHA_1[J]     = 2.0 * PI * (B0[J] + (B1[J] * source_temp[J]));
    ALPHA_2[J]     = 2.0 * PI * B1[J];

    SIGMA_1[J] = 2.0 * PI * (B0[J] - (B1[J] * source_temp[J]));
    SIGMA_2[J] = 2.0 * PI * B1[J];
  }

  
  INTENSITY_DOWN[0] = BB_TOP_OF_ATM * exp(-TAULS[0]) + \
                      SOURCE_J[0]/(LAMBDAS[0] + 1.0) * (1.0 - exp(-TAULS[0]*LAMBDAS[0]+1.0)) + \
                      SOURCE_K[0]/(LAMBDAS[0] - 1.0) * (exp(-TAULS[0]) - exp(-TAULS[0]*LAMBDAS[0])) + \
                      SIGMA_1[0] * (1.0 - exp(-TAULS[0])) + \
                      SIGMA_2[0] * (exp(-TAULS[0]) + TAULS[0] + 1.0);


  INTENSITY_UP[NEW_NLAYER-1] = 2.0 * BB_BOTTOM_OF_ATM * EMIS * PI;

  // Do the downward intensity first
  for(J=1; J<NEW_NLAYER; J++)
  {
    INTENSITY_DOWN[J] = INTENSITY_DOWN[J-1] * exp(-TAULS[J]) + \
                    SOURCE_J[J]/(LAMBDAS[J] + 1.0) * (1.0 - exp(-TAULS[J]*LAMBDAS[J]+1.0)) + \
                   SOURCE_K[J]/(LAMBDAS[J] - 1.0) * (exp(-TAULS[J]) - exp(-TAULS[J]*LAMBDAS[J])) + \
                    SIGMA_1[J] * (1.0 - exp(-TAULS[J])) + \
                    SIGMA_2[J] * (exp(-TAULS[J]) + TAULS[J] - 1.0);
  }

  // Calculate the upward intensity next
  for(Z=1; Z<NEW_NLAYER; Z++)
  {
    J = NEW_NLAYER - Z - 1;
    INTENSITY_UP[J] = INTENSITY_UP[J+1] * exp(-TAULS[J+1]) + \
                      SOURCE_G[J+1]/(LAMBDAS[J+1]-1.0)*(exp(-TAULS[J+1])-exp(-TAULS[J+1]*LAMBDAS[J+1])) + \
                      SOURCE_H[J+1]/(LAMBDAS[J+1]+1.0) * (1.0 - exp(-TAULS[J+1] * (LAMBDAS[J+1] + 1.0))) + \
                      ALPHA_1[J+1] * (1.0 - exp(-TAULS[J+1])) + \
                      ALPHA_2[J+1] * (1.0 - ((TAULS[J+1] + 1.0) * (exp(-TAULS[J+1]))));
  }
  
  printf("\n\n\n\n\n\n\n");
  for(J=1; J<NEW_NLAYER; J++)
  {
    printf("%.8e, %.8e,\n", TMI[NEW_NLAYER - J - 1] / (4.0 * PI), INTENSITY_UP[J] / (2.0 * PI));
  }


  TWO_STREAM_INTENSITY = TMI[0] / (4.0 * PI);
  SOURCE_INTENSITY     = INTENSITY_UP[0] / (2.0 * PI);
  // Define the energy from other sources (eq 37 and 38, Toon)
  if (NU > 430.0e12)
  {
    return TWO_STREAM_INTENSITY;
  }
  else
  {
    return SOURCE_INTENSITY;
  } 
}