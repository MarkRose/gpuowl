// Copyright (C) Mihai Preda and George Woltman

#if FFT_FP64

// 902 FMA + 90 ADD.  The method of fft11.cl (prime95's gwnum/zr11.mac): pair inputs k and 31-k, then each
// output pair j, 31-j is (u0 + sum_k cos(j*k*tau/31) sum_k) +- sum_k sin(j*k*tau/31) dif_k, with w = e^(2*pi*i/31).  As in fft11.cl, cbase = u0 +
// cos(8*tau/31) * sum zeroes one term per cos row (cos(8*tau/31) = -0.051 is the cosine nearest the middle of the range: every coefficient
// stays below 1.03 in magnitude), and the sin rows are scaled by sin(9*tau/31), the most accurately represented of the 15 sines
// (relative rounding error 1.6e-17).
// Accumulation order as fft13.cl: cos terms by |coef| rank 1,0,3,2,5,4,7,6,9,8,11,10,13,12 (ranks from the smallest), sin terms smallest first.

void OVERLOAD fft31(T2 *u) {
  const double COS8 = -0.05064916883871271227875; // cos(8*tau/31): subtracted from every cos coefficient, it zeroes one term per cos row
  const double C1_C8 = 1.030179110091207206217;   // cos(1*tau/31) - cos(8*tau/31)
  const double C2_C8 = 0.9696069804589433414059;  // cos(2*tau/31) - cos(8*tau/31)
  const double C3_C8 = 0.8714126100459890386423;  // cos(3*tau/31) - cos(8*tau/31)
  const double C4_C8 = 0.7396160879143992800796;  // cos(4*tau/31) - cos(8*tau/31)
  const double C5_C8 = 0.5796131791656751696442;  // cos(5*tau/31) - cos(8*tau/31)
  const double C6_C8 = 0.3979544216835329978206;  // cos(6*tau/31) - cos(8*tau/31)
  const double C7_C8 = 0.2020769463432893759362;  // cos(7*tau/31) - cos(8*tau/31)
  const double C9_C8 = -0.2000033634200078270361; // cos(9*tau/31) - cos(8*tau/31)
  const double C10_C8 = -0.3897449827189215972374;// cos(10*tau/31) - cos(8*tau/31)
  const double C11_C8 = -0.5614568137089501318680;// cos(11*tau/31) - cos(8*tau/31)
  const double C12_C8 = -0.7081089538540781896345;// cos(12*tau/31) - cos(8*tau/31)
  const double C13_C8 = -0.8236974473058694065487;// cos(13*tau/31) - cos(8*tau/31)
  const double C14_C8 = -0.9034900875613361391971;// cos(14*tau/31) - cos(8*tau/31)
  const double C15_C8 = -0.9442201545531824340426;// cos(15*tau/31) - cos(8*tau/31)
  const double SIN9 = 0.9680771188662043051530;   // sin(9*tau/31): factored out of the sin rows, it makes one term per row a plain add
  const double S1_S9 = 0.2079364506873352525185;  // sin(1*tau/31) / sin(9*tau/31)
  const double S2_S9 = 0.4073599586520354364685;  // sin(2*tau/31) / sin(9*tau/31)
  const double S3_S9 = 0.5901061020467584619302;  // sin(3*tau/31) / sin(9*tau/31)
  const double S4_S9 = 0.7486932322891642376550;  // sin(4*tau/31) / sin(9*tau/31)
  const double S5_S9 = 0.8766287736339320569721;  // sin(5*tau/31) / sin(9*tau/31)
  const double S6_S9 = 0.9686750299866192862133;  // sin(6*tau/31) / sin(9*tau/31)
  const double S7_S9 = 1.021063616797171006429;   // sin(7*tau/31) / sin(9*tau/31)
  const double S8_S9 = 1.031649739166165645980;   // sin(8*tau/31) / sin(9*tau/31)
  const double S10_S9 = 0.9274101433388233418963; // sin(10*tau/31) / sin(9*tau/31)
  const double S11_S9 = 0.8168520064432902518249; // sin(11*tau/31) / sin(9*tau/31)
  const double S12_S9 = 0.6728518524279333665138; // sin(12*tau/31) / sin(9*tau/31)
  const double S13_S9 = 0.5013050645174410786296; // sin(13*tau/31) / sin(9*tau/31)
  const double S14_S9 = 0.3092347883647606770760; // sin(14*tau/31) / sin(9*tau/31)
  const double S15_S9 = 0.1045044036428821091535; // sin(15*tau/31) / sin(9*tau/31)

  X2_mul_t4(u[1], u[30]);
  X2_mul_t4(u[2], u[29]);
  X2_mul_t4(u[3], u[28]);
  X2_mul_t4(u[4], u[27]);
  X2_mul_t4(u[5], u[26]);
  X2_mul_t4(u[6], u[25]);
  X2_mul_t4(u[7], u[24]);
  X2_mul_t4(u[8], u[23]);
  X2_mul_t4(u[9], u[22]);
  X2_mul_t4(u[10], u[21]);
  X2_mul_t4(u[11], u[20]);
  X2_mul_t4(u[12], u[19]);
  X2_mul_t4(u[13], u[18]);
  X2_mul_t4(u[14], u[17]);
  X2_mul_t4(u[15], u[16]);

  T2 sum = u[1] + u[2] + u[3] + u[4] + u[5] + u[6] + u[7] + u[8] + u[9] + u[10] + u[11] + u[12] + u[13] + u[14] + u[15];
  T2 cbase = fmaT2(COS8, sum, u[0]);                    // u[0] + COS8 * sum: each cos row below then has one term fewer
  u[0] = u[0] + sum;

  T2 tmp130a = fmaT2(C7_C8, u[7], cbase);
  T2 tmp229a = fmaT2(C7_C8, u[12], cbase);
  T2 tmp328a = fmaT2(C7_C8, u[8], cbase);
  T2 tmp427a = fmaT2(C7_C8, u[6], cbase);
  T2 tmp526a = fmaT2(C7_C8, u[11], cbase);
  T2 tmp625a = fmaT2(C7_C8, u[4], cbase);
  T2 tmp724a = fmaT2(C7_C8, u[1], cbase);
  T2 tmp823a = fmaT2(C7_C8, u[3], cbase);
  T2 tmp922a = fmaT2(C7_C8, u[13], cbase);
  T2 tmp1021a = fmaT2(C7_C8, u[10], cbase);
  T2 tmp1120a = fmaT2(C7_C8, u[5], cbase);
  T2 tmp1219a = fmaT2(C7_C8, u[2], cbase);
  T2 tmp1318a = fmaT2(C7_C8, u[9], cbase);
  T2 tmp1417a = fmaT2(C7_C8, u[15], cbase);
  T2 tmp1516a = fmaT2(C7_C8, u[14], cbase);

  tmp130a = fmaT2(C9_C8, u[9], tmp130a);
  tmp229a = fmaT2(C9_C8, u[11], tmp229a);
  tmp328a = fmaT2(C9_C8, u[3], tmp328a);
  tmp427a = fmaT2(C9_C8, u[10], tmp427a);
  tmp526a = fmaT2(C9_C8, u[8], tmp526a);
  tmp625a = fmaT2(C9_C8, u[14], tmp625a);
  tmp724a = fmaT2(C9_C8, u[12], tmp724a);
  tmp823a = fmaT2(C9_C8, u[5], tmp823a);
  tmp922a = fmaT2(C9_C8, u[1], tmp922a);
  tmp1021a = fmaT2(C9_C8, u[4], tmp1021a);
  tmp1120a = fmaT2(C9_C8, u[2], tmp1120a);
  tmp1219a = fmaT2(C9_C8, u[7], tmp1219a);
  tmp1318a = fmaT2(C9_C8, u[15], tmp1318a);
  tmp1417a = fmaT2(C9_C8, u[6], tmp1417a);
  tmp1516a = fmaT2(C9_C8, u[13], tmp1516a);

  tmp130a = fmaT2(C6_C8, u[6], tmp130a);
  tmp229a = fmaT2(C6_C8, u[3], tmp229a);
  tmp328a = fmaT2(C6_C8, u[2], tmp328a);
  tmp427a = fmaT2(C6_C8, u[14], tmp427a);
  tmp526a = fmaT2(C6_C8, u[5], tmp526a);
  tmp625a = fmaT2(C6_C8, u[1], tmp625a);
  tmp724a = fmaT2(C6_C8, u[8], tmp724a);
  tmp823a = fmaT2(C6_C8, u[7], tmp823a);
  tmp922a = fmaT2(C6_C8, u[11], tmp922a);
  tmp1021a = fmaT2(C6_C8, u[13], tmp1021a);
  tmp1120a = fmaT2(C6_C8, u[9], tmp1120a);
  tmp1219a = fmaT2(C6_C8, u[15], tmp1219a);
  tmp1318a = fmaT2(C6_C8, u[10], tmp1318a);
  tmp1417a = fmaT2(C6_C8, u[4], tmp1417a);
  tmp1516a = fmaT2(C6_C8, u[12], tmp1516a);

  tmp130a = fmaT2(C10_C8, u[10], tmp130a);
  tmp229a = fmaT2(C10_C8, u[5], tmp229a);
  tmp328a = fmaT2(C10_C8, u[7], tmp328a);
  tmp427a = fmaT2(C10_C8, u[13], tmp427a);
  tmp526a = fmaT2(C10_C8, u[2], tmp526a);
  tmp625a = fmaT2(C10_C8, u[12], tmp625a);
  tmp724a = fmaT2(C10_C8, u[3], tmp724a);
  tmp823a = fmaT2(C10_C8, u[9], tmp823a);
  tmp922a = fmaT2(C10_C8, u[8], tmp922a);
  tmp1021a = fmaT2(C10_C8, u[1], tmp1021a);
  tmp1120a = fmaT2(C10_C8, u[15], tmp1120a);
  tmp1219a = fmaT2(C10_C8, u[6], tmp1219a);
  tmp1318a = fmaT2(C10_C8, u[4], tmp1318a);
  tmp1417a = fmaT2(C10_C8, u[14], tmp1417a);
  tmp1516a = fmaT2(C10_C8, u[11], tmp1516a);

  tmp130a = fmaT2(C5_C8, u[5], tmp130a);
  tmp229a = fmaT2(C5_C8, u[13], tmp229a);
  tmp328a = fmaT2(C5_C8, u[12], tmp328a);
  tmp427a = fmaT2(C5_C8, u[9], tmp427a);
  tmp526a = fmaT2(C5_C8, u[1], tmp526a);
  tmp625a = fmaT2(C5_C8, u[6], tmp625a);
  tmp724a = fmaT2(C5_C8, u[14], tmp724a);
  tmp823a = fmaT2(C5_C8, u[11], tmp823a);
  tmp922a = fmaT2(C5_C8, u[4], tmp922a);
  tmp1021a = fmaT2(C5_C8, u[15], tmp1021a);
  tmp1120a = fmaT2(C5_C8, u[8], tmp1120a);
  tmp1219a = fmaT2(C5_C8, u[3], tmp1219a);
  tmp1318a = fmaT2(C5_C8, u[2], tmp1318a);
  tmp1417a = fmaT2(C5_C8, u[7], tmp1417a);
  tmp1516a = fmaT2(C5_C8, u[10], tmp1516a);

  tmp130a = fmaT2(C11_C8, u[11], tmp130a);
  tmp229a = fmaT2(C11_C8, u[10], tmp229a);
  tmp328a = fmaT2(C11_C8, u[14], tmp328a);
  tmp427a = fmaT2(C11_C8, u[5], tmp427a);
  tmp526a = fmaT2(C11_C8, u[4], tmp526a);
  tmp625a = fmaT2(C11_C8, u[7], tmp625a);
  tmp724a = fmaT2(C11_C8, u[6], tmp724a);
  tmp823a = fmaT2(C11_C8, u[13], tmp823a);
  tmp922a = fmaT2(C11_C8, u[15], tmp922a);
  tmp1021a = fmaT2(C11_C8, u[2], tmp1021a);
  tmp1120a = fmaT2(C11_C8, u[1], tmp1120a);
  tmp1219a = fmaT2(C11_C8, u[12], tmp1219a);
  tmp1318a = fmaT2(C11_C8, u[8], tmp1318a);
  tmp1417a = fmaT2(C11_C8, u[3], tmp1417a);
  tmp1516a = fmaT2(C11_C8, u[9], tmp1516a);

  tmp130a = fmaT2(C4_C8, u[4], tmp130a);
  tmp229a = fmaT2(C4_C8, u[2], tmp229a);
  tmp328a = fmaT2(C4_C8, u[9], tmp328a);
  tmp427a = fmaT2(C4_C8, u[1], tmp427a);
  tmp526a = fmaT2(C4_C8, u[7], tmp526a);
  tmp625a = fmaT2(C4_C8, u[11], tmp625a);
  tmp724a = fmaT2(C4_C8, u[5], tmp724a);
  tmp823a = fmaT2(C4_C8, u[15], tmp823a);
  tmp922a = fmaT2(C4_C8, u[3], tmp922a);
  tmp1021a = fmaT2(C4_C8, u[12], tmp1021a);
  tmp1120a = fmaT2(C4_C8, u[6], tmp1120a);
  tmp1219a = fmaT2(C4_C8, u[10], tmp1219a);
  tmp1318a = fmaT2(C4_C8, u[14], tmp1318a);
  tmp1417a = fmaT2(C4_C8, u[13], tmp1417a);
  tmp1516a = fmaT2(C4_C8, u[8], tmp1516a);

  tmp130a = fmaT2(C12_C8, u[12], tmp130a);
  tmp229a = fmaT2(C12_C8, u[6], tmp229a);
  tmp328a = fmaT2(C12_C8, u[4], tmp328a);
  tmp427a = fmaT2(C12_C8, u[3], tmp427a);
  tmp526a = fmaT2(C12_C8, u[10], tmp526a);
  tmp625a = fmaT2(C12_C8, u[2], tmp625a);
  tmp724a = fmaT2(C12_C8, u[15], tmp724a);
  tmp823a = fmaT2(C12_C8, u[14], tmp823a);
  tmp922a = fmaT2(C12_C8, u[9], tmp922a);
  tmp1021a = fmaT2(C12_C8, u[5], tmp1021a);
  tmp1120a = fmaT2(C12_C8, u[13], tmp1120a);
  tmp1219a = fmaT2(C12_C8, u[1], tmp1219a);
  tmp1318a = fmaT2(C12_C8, u[11], tmp1318a);
  tmp1417a = fmaT2(C12_C8, u[8], tmp1417a);
  tmp1516a = fmaT2(C12_C8, u[7], tmp1516a);

  tmp130a = fmaT2(C3_C8, u[3], tmp130a);
  tmp229a = fmaT2(C3_C8, u[14], tmp229a);
  tmp328a = fmaT2(C3_C8, u[1], tmp328a);
  tmp427a = fmaT2(C3_C8, u[7], tmp427a);
  tmp526a = fmaT2(C3_C8, u[13], tmp526a);
  tmp625a = fmaT2(C3_C8, u[15], tmp625a);
  tmp724a = fmaT2(C3_C8, u[4], tmp724a);
  tmp823a = fmaT2(C3_C8, u[12], tmp823a);
  tmp922a = fmaT2(C3_C8, u[10], tmp922a);
  tmp1021a = fmaT2(C3_C8, u[9], tmp1021a);
  tmp1120a = fmaT2(C3_C8, u[11], tmp1120a);
  tmp1219a = fmaT2(C3_C8, u[8], tmp1219a);
  tmp1318a = fmaT2(C3_C8, u[5], tmp1318a);
  tmp1417a = fmaT2(C3_C8, u[2], tmp1417a);
  tmp1516a = fmaT2(C3_C8, u[6], tmp1516a);

  tmp130a = fmaT2(C13_C8, u[13], tmp130a);
  tmp229a = fmaT2(C13_C8, u[9], tmp229a);
  tmp328a = fmaT2(C13_C8, u[6], tmp328a);
  tmp427a = fmaT2(C13_C8, u[11], tmp427a);
  tmp526a = fmaT2(C13_C8, u[15], tmp526a);
  tmp625a = fmaT2(C13_C8, u[3], tmp625a);
  tmp724a = fmaT2(C13_C8, u[7], tmp724a);
  tmp823a = fmaT2(C13_C8, u[10], tmp823a);
  tmp922a = fmaT2(C13_C8, u[2], tmp922a);
  tmp1021a = fmaT2(C13_C8, u[8], tmp1021a);
  tmp1120a = fmaT2(C13_C8, u[4], tmp1120a);
  tmp1219a = fmaT2(C13_C8, u[14], tmp1219a);
  tmp1318a = fmaT2(C13_C8, u[1], tmp1318a);
  tmp1417a = fmaT2(C13_C8, u[12], tmp1417a);
  tmp1516a = fmaT2(C13_C8, u[5], tmp1516a);

  tmp130a = fmaT2(C15_C8, u[15], tmp130a);
  tmp229a = fmaT2(C15_C8, u[8], tmp229a);
  tmp328a = fmaT2(C15_C8, u[5], tmp328a);
  tmp427a = fmaT2(C15_C8, u[4], tmp427a);
  tmp526a = fmaT2(C15_C8, u[3], tmp526a);
  tmp625a = fmaT2(C15_C8, u[13], tmp625a);
  tmp724a = fmaT2(C15_C8, u[11], tmp724a);
  tmp823a = fmaT2(C15_C8, u[2], tmp823a);
  tmp922a = fmaT2(C15_C8, u[12], tmp922a);
  tmp1021a = fmaT2(C15_C8, u[14], tmp1021a);
  tmp1120a = fmaT2(C15_C8, u[7], tmp1120a);
  tmp1219a = fmaT2(C15_C8, u[9], tmp1219a);
  tmp1318a = fmaT2(C15_C8, u[6], tmp1318a);
  tmp1417a = fmaT2(C15_C8, u[10], tmp1417a);
  tmp1516a = fmaT2(C15_C8, u[1], tmp1516a);

  tmp130a = fmaT2(C14_C8, u[14], tmp130a);
  tmp229a = fmaT2(C14_C8, u[7], tmp229a);
  tmp328a = fmaT2(C14_C8, u[15], tmp328a);
  tmp427a = fmaT2(C14_C8, u[12], tmp427a);
  tmp526a = fmaT2(C14_C8, u[9], tmp526a);
  tmp625a = fmaT2(C14_C8, u[8], tmp625a);
  tmp724a = fmaT2(C14_C8, u[2], tmp724a);
  tmp823a = fmaT2(C14_C8, u[6], tmp823a);
  tmp922a = fmaT2(C14_C8, u[5], tmp922a);
  tmp1021a = fmaT2(C14_C8, u[11], tmp1021a);
  tmp1120a = fmaT2(C14_C8, u[10], tmp1120a);
  tmp1219a = fmaT2(C14_C8, u[4], tmp1219a);
  tmp1318a = fmaT2(C14_C8, u[13], tmp1318a);
  tmp1417a = fmaT2(C14_C8, u[1], tmp1417a);
  tmp1516a = fmaT2(C14_C8, u[3], tmp1516a);

  tmp130a = fmaT2(C1_C8, u[1], tmp130a);
  tmp229a = fmaT2(C1_C8, u[15], tmp229a);
  tmp328a = fmaT2(C1_C8, u[10], tmp328a);
  tmp427a = fmaT2(C1_C8, u[8], tmp427a);
  tmp526a = fmaT2(C1_C8, u[6], tmp526a);
  tmp625a = fmaT2(C1_C8, u[5], tmp625a);
  tmp724a = fmaT2(C1_C8, u[9], tmp724a);
  tmp823a = fmaT2(C1_C8, u[4], tmp823a);
  tmp922a = fmaT2(C1_C8, u[7], tmp922a);
  tmp1021a = fmaT2(C1_C8, u[3], tmp1021a);
  tmp1120a = fmaT2(C1_C8, u[14], tmp1120a);
  tmp1219a = fmaT2(C1_C8, u[13], tmp1219a);
  tmp1318a = fmaT2(C1_C8, u[12], tmp1318a);
  tmp1417a = fmaT2(C1_C8, u[11], tmp1417a);
  tmp1516a = fmaT2(C1_C8, u[2], tmp1516a);

  tmp130a = fmaT2(C2_C8, u[2], tmp130a);
  tmp229a = fmaT2(C2_C8, u[1], tmp229a);
  tmp328a = fmaT2(C2_C8, u[11], tmp328a);
  tmp427a = fmaT2(C2_C8, u[15], tmp427a);
  tmp526a = fmaT2(C2_C8, u[12], tmp526a);
  tmp625a = fmaT2(C2_C8, u[10], tmp625a);
  tmp724a = fmaT2(C2_C8, u[13], tmp724a);
  tmp823a = fmaT2(C2_C8, u[8], tmp823a);
  tmp922a = fmaT2(C2_C8, u[14], tmp922a);
  tmp1021a = fmaT2(C2_C8, u[6], tmp1021a);
  tmp1120a = fmaT2(C2_C8, u[3], tmp1120a);
  tmp1219a = fmaT2(C2_C8, u[5], tmp1219a);
  tmp1318a = fmaT2(C2_C8, u[7], tmp1318a);
  tmp1417a = fmaT2(C2_C8, u[9], tmp1417a);
  tmp1516a = fmaT2(C2_C8, u[4], tmp1516a);

  T2 tmp130b = fmaT2(S15_S9, u[16], u[22]);
  T2 tmp229b = fmaT2(S15_S9, -u[23], -u[20]);
  T2 tmp328b = fmaT2(S15_S9, u[26], u[28]);
  T2 tmp427b = fmaT2(S15_S9, -u[27], u[21]);
  T2 tmp526b = fmaT2(S15_S9, u[28], u[23]);
  T2 tmp625b = fmaT2(S15_S9, -u[18], -u[17]);
  T2 tmp724b = fmaT2(S15_S9, u[20], -u[19]);
  T2 tmp823b = fmaT2(S15_S9, -u[29], u[26]);
  T2 tmp922b = fmaT2(S15_S9, u[19], u[30]);
  T2 tmp1021b = fmaT2(S15_S9, -u[17], u[27]);
  T2 tmp1120b = fmaT2(S15_S9, u[24], -u[29]);
  T2 tmp1219b = fmaT2(S15_S9, u[22], -u[24]);
  T2 tmp1318b = fmaT2(S15_S9, -u[25], u[16]);
  T2 tmp1417b = fmaT2(S15_S9, -u[21], -u[25]);
  T2 tmp1516b = fmaT2(S15_S9, u[30], u[18]);

  tmp130b = fmaT2(S1_S9, u[30], tmp130b);
  tmp229b = fmaT2(S1_S9, -u[16], tmp229b);
  tmp328b = fmaT2(S1_S9, -u[21], tmp328b);
  tmp427b = fmaT2(S1_S9, u[23], tmp427b);
  tmp526b = fmaT2(S1_S9, -u[25], tmp526b);
  tmp625b = fmaT2(S1_S9, -u[26], tmp625b);
  tmp724b = fmaT2(S1_S9, u[22], tmp724b);
  tmp823b = fmaT2(S1_S9, u[27], tmp823b);
  tmp922b = fmaT2(S1_S9, u[24], tmp922b);
  tmp1021b = fmaT2(S1_S9, -u[28], tmp1021b);
  tmp1120b = fmaT2(S1_S9, -u[17], tmp1120b);
  tmp1219b = fmaT2(S1_S9, u[18], tmp1219b);
  tmp1318b = fmaT2(S1_S9, u[19], tmp1318b);
  tmp1417b = fmaT2(S1_S9, -u[20], tmp1417b);
  tmp1516b = fmaT2(S1_S9, -u[29], tmp1516b);

  tmp130b = fmaT2(S14_S9, u[17], tmp130b);
  tmp229b = fmaT2(S14_S9, u[24], tmp229b);
  tmp328b = fmaT2(S14_S9, u[16], tmp328b);
  tmp427b = fmaT2(S14_S9, -u[19], tmp427b);
  tmp526b = fmaT2(S14_S9, u[22], tmp526b);
  tmp625b = fmaT2(S14_S9, -u[23], tmp625b);
  tmp724b = fmaT2(S14_S9, u[29], tmp724b);
  tmp823b = fmaT2(S14_S9, -u[25], tmp823b);
  tmp922b = fmaT2(S14_S9, u[26], tmp922b);
  tmp1021b = fmaT2(S14_S9, -u[20], tmp1021b);
  tmp1120b = fmaT2(S14_S9, -u[21], tmp1120b);
  tmp1219b = fmaT2(S14_S9, -u[27], tmp1219b);
  tmp1318b = fmaT2(S14_S9, u[18], tmp1318b);
  tmp1417b = fmaT2(S14_S9, u[30], tmp1417b);
  tmp1516b = fmaT2(S14_S9, u[28], tmp1516b);

  tmp130b = fmaT2(S2_S9, u[29], tmp130b);
  tmp229b = fmaT2(S2_S9, u[30], tmp229b);
  tmp328b = fmaT2(S2_S9, u[20], tmp328b);
  tmp427b = fmaT2(S2_S9, -u[16], tmp427b);
  tmp526b = fmaT2(S2_S9, -u[19], tmp526b);
  tmp625b = fmaT2(S2_S9, -u[21], tmp625b);
  tmp724b = fmaT2(S2_S9, -u[18], tmp724b);
  tmp823b = fmaT2(S2_S9, u[23], tmp823b);
  tmp922b = fmaT2(S2_S9, u[17], tmp922b);
  tmp1021b = fmaT2(S2_S9, -u[25], tmp1021b);
  tmp1120b = fmaT2(S2_S9, u[28], tmp1120b);
  tmp1219b = fmaT2(S2_S9, -u[26], tmp1219b);
  tmp1318b = fmaT2(S2_S9, -u[24], tmp1318b);
  tmp1417b = fmaT2(S2_S9, u[22], tmp1417b);
  tmp1516b = fmaT2(S2_S9, -u[27], tmp1516b);

  tmp130b = fmaT2(S13_S9, u[18], tmp130b);
  tmp229b = fmaT2(S13_S9, -u[22], tmp229b);
  tmp328b = fmaT2(S13_S9, -u[25], tmp328b);
  tmp427b = fmaT2(S13_S9, u[20], tmp427b);
  tmp526b = fmaT2(S13_S9, u[16], tmp526b);
  tmp625b = fmaT2(S13_S9, -u[28], tmp625b);
  tmp724b = fmaT2(S13_S9, -u[24], tmp724b);
  tmp823b = fmaT2(S13_S9, -u[21], tmp823b);
  tmp922b = fmaT2(S13_S9, -u[29], tmp922b);
  tmp1021b = fmaT2(S13_S9, -u[23], tmp1021b);
  tmp1120b = fmaT2(S13_S9, u[27], tmp1120b);
  tmp1219b = fmaT2(S13_S9, u[17], tmp1219b);
  tmp1318b = fmaT2(S13_S9, u[30], tmp1318b);
  tmp1417b = fmaT2(S13_S9, u[19], tmp1417b);
  tmp1516b = fmaT2(S13_S9, u[26], tmp1516b);

  tmp130b = fmaT2(S3_S9, u[28], tmp130b);
  tmp229b = fmaT2(S3_S9, -u[17], tmp229b);
  tmp328b = fmaT2(S3_S9, u[30], tmp328b);
  tmp427b = fmaT2(S3_S9, -u[24], tmp427b);
  tmp526b = fmaT2(S3_S9, u[18], tmp526b);
  tmp625b = fmaT2(S3_S9, -u[16], tmp625b);
  tmp724b = fmaT2(S3_S9, -u[27], tmp724b);
  tmp823b = fmaT2(S3_S9, u[19], tmp823b);
  tmp922b = fmaT2(S3_S9, -u[21], tmp922b);
  tmp1021b = fmaT2(S3_S9, -u[22], tmp1021b);
  tmp1120b = fmaT2(S3_S9, -u[20], tmp1120b);
  tmp1219b = fmaT2(S3_S9, u[23], tmp1219b);
  tmp1318b = fmaT2(S3_S9, u[26], tmp1318b);
  tmp1417b = fmaT2(S3_S9, -u[29], tmp1417b);
  tmp1516b = fmaT2(S3_S9, -u[25], tmp1516b);

  tmp130b = fmaT2(S12_S9, u[19], tmp130b);
  tmp229b = fmaT2(S12_S9, u[25], tmp229b);
  tmp328b = fmaT2(S12_S9, u[27], tmp328b);
  tmp427b = fmaT2(S12_S9, u[28], tmp427b);
  tmp526b = fmaT2(S12_S9, -u[21], tmp526b);
  tmp625b = fmaT2(S12_S9, u[29], tmp625b);
  tmp724b = fmaT2(S12_S9, u[16], tmp724b);
  tmp823b = fmaT2(S12_S9, -u[17], tmp823b);
  tmp922b = fmaT2(S12_S9, -u[22], tmp922b);
  tmp1021b = fmaT2(S12_S9, -u[26], tmp1021b);
  tmp1120b = fmaT2(S12_S9, -u[18], tmp1120b);
  tmp1219b = fmaT2(S12_S9, u[30], tmp1219b);
  tmp1318b = fmaT2(S12_S9, -u[20], tmp1318b);
  tmp1417b = fmaT2(S12_S9, -u[23], tmp1417b);
  tmp1516b = fmaT2(S12_S9, u[24], tmp1516b);

  tmp130b = fmaT2(S4_S9, u[27], tmp130b);
  tmp229b = fmaT2(S4_S9, u[29], tmp229b);
  tmp328b = fmaT2(S4_S9, -u[22], tmp328b);
  tmp427b = fmaT2(S4_S9, u[30], tmp427b);
  tmp526b = fmaT2(S4_S9, u[24], tmp526b);
  tmp625b = fmaT2(S4_S9, u[20], tmp625b);
  tmp724b = fmaT2(S4_S9, u[26], tmp724b);
  tmp823b = fmaT2(S4_S9, -u[16], tmp823b);
  tmp922b = fmaT2(S4_S9, -u[28], tmp922b);
  tmp1021b = fmaT2(S4_S9, -u[19], tmp1021b);
  tmp1120b = fmaT2(S4_S9, u[25], tmp1120b);
  tmp1219b = fmaT2(S4_S9, -u[21], tmp1219b);
  tmp1318b = fmaT2(S4_S9, -u[17], tmp1318b);
  tmp1417b = fmaT2(S4_S9, -u[18], tmp1417b);
  tmp1516b = fmaT2(S4_S9, -u[23], tmp1516b);

  tmp130b = fmaT2(S11_S9, u[20], tmp130b);
  tmp229b = fmaT2(S11_S9, -u[21], tmp229b);
  tmp328b = fmaT2(S11_S9, u[17], tmp328b);
  tmp427b = fmaT2(S11_S9, -u[26], tmp427b);
  tmp526b = fmaT2(S11_S9, -u[27], tmp526b);
  tmp625b = fmaT2(S11_S9, u[24], tmp625b);
  tmp724b = fmaT2(S11_S9, u[25], tmp724b);
  tmp823b = fmaT2(S11_S9, u[18], tmp823b);
  tmp922b = fmaT2(S11_S9, u[16], tmp922b);
  tmp1021b = fmaT2(S11_S9, -u[29], tmp1021b);
  tmp1120b = fmaT2(S11_S9, u[30], tmp1120b);
  tmp1219b = fmaT2(S11_S9, -u[19], tmp1219b);
  tmp1318b = fmaT2(S11_S9, u[23], tmp1318b);
  tmp1417b = fmaT2(S11_S9, u[28], tmp1417b);
  tmp1516b = fmaT2(S11_S9, u[22], tmp1516b);

  tmp130b = fmaT2(S5_S9, u[26], tmp130b);
  tmp229b = fmaT2(S5_S9, -u[18], tmp229b);
  tmp328b = fmaT2(S5_S9, u[19], tmp328b);
  tmp427b = fmaT2(S5_S9, u[22], tmp427b);
  tmp526b = fmaT2(S5_S9, u[30], tmp526b);
  tmp625b = fmaT2(S5_S9, u[25], tmp625b);
  tmp724b = fmaT2(S5_S9, u[17], tmp724b);
  tmp823b = fmaT2(S5_S9, -u[20], tmp823b);
  tmp922b = fmaT2(S5_S9, u[27], tmp922b);
  tmp1021b = fmaT2(S5_S9, -u[16], tmp1021b);
  tmp1120b = fmaT2(S5_S9, -u[23], tmp1120b);
  tmp1219b = fmaT2(S5_S9, u[28], tmp1219b);
  tmp1318b = fmaT2(S5_S9, -u[29], tmp1318b);
  tmp1417b = fmaT2(S5_S9, u[24], tmp1417b);
  tmp1516b = fmaT2(S5_S9, -u[21], tmp1516b);

  tmp130b = fmaT2(S10_S9, u[21], tmp130b);
  tmp229b = fmaT2(S10_S9, u[26], tmp229b);
  tmp328b = fmaT2(S10_S9, -u[24], tmp328b);
  tmp427b = fmaT2(S10_S9, -u[18], tmp427b);
  tmp526b = fmaT2(S10_S9, u[29], tmp526b);
  tmp625b = fmaT2(S10_S9, u[19], tmp625b);
  tmp724b = fmaT2(S10_S9, -u[28], tmp724b);
  tmp823b = fmaT2(S10_S9, u[22], tmp823b);
  tmp922b = fmaT2(S10_S9, u[23], tmp922b);
  tmp1021b = fmaT2(S10_S9, u[30], tmp1021b);
  tmp1120b = fmaT2(S10_S9, u[16], tmp1120b);
  tmp1219b = fmaT2(S10_S9, u[25], tmp1219b);
  tmp1318b = fmaT2(S10_S9, -u[27], tmp1318b);
  tmp1417b = fmaT2(S10_S9, u[17], tmp1417b);
  tmp1516b = fmaT2(S10_S9, u[20], tmp1516b);

  tmp130b = fmaT2(S6_S9, u[25], tmp130b);
  tmp229b = fmaT2(S6_S9, u[28], tmp229b);
  tmp328b = fmaT2(S6_S9, u[29], tmp328b);
  tmp427b = fmaT2(S6_S9, -u[17], tmp427b);
  tmp526b = fmaT2(S6_S9, -u[26], tmp526b);
  tmp625b = fmaT2(S6_S9, u[30], tmp625b);
  tmp724b = fmaT2(S6_S9, -u[23], tmp724b);
  tmp823b = fmaT2(S6_S9, -u[24], tmp823b);
  tmp922b = fmaT2(S6_S9, u[20], tmp922b);
  tmp1021b = fmaT2(S6_S9, u[18], tmp1021b);
  tmp1120b = fmaT2(S6_S9, u[22], tmp1120b);
  tmp1219b = fmaT2(S6_S9, -u[16], tmp1219b);
  tmp1318b = fmaT2(S6_S9, u[21], tmp1318b);
  tmp1417b = fmaT2(S6_S9, -u[27], tmp1417b);
  tmp1516b = fmaT2(S6_S9, -u[19], tmp1516b);

  tmp130b = fmaT2(S7_S9, u[24], tmp130b);
  tmp229b = fmaT2(S7_S9, -u[19], tmp229b);
  tmp328b = fmaT2(S7_S9, -u[23], tmp328b);
  tmp427b = fmaT2(S7_S9, -u[25], tmp427b);
  tmp526b = fmaT2(S7_S9, -u[20], tmp526b);
  tmp625b = fmaT2(S7_S9, -u[27], tmp625b);
  tmp724b = fmaT2(S7_S9, u[30], tmp724b);
  tmp823b = fmaT2(S7_S9, -u[28], tmp823b);
  tmp922b = fmaT2(S7_S9, -u[18], tmp922b);
  tmp1021b = fmaT2(S7_S9, u[21], tmp1021b);
  tmp1120b = fmaT2(S7_S9, -u[26], tmp1120b);
  tmp1219b = fmaT2(S7_S9, -u[29], tmp1219b);
  tmp1318b = fmaT2(S7_S9, -u[22], tmp1318b);
  tmp1417b = fmaT2(S7_S9, -u[16], tmp1417b);
  tmp1516b = fmaT2(S7_S9, -u[17], tmp1516b);

  tmp130b = fmaT2(S8_S9, u[23], tmp130b);
  tmp229b = fmaT2(S8_S9, u[27], tmp229b);
  tmp328b = fmaT2(S8_S9, u[18], tmp328b);
  tmp427b = fmaT2(S8_S9, u[29], tmp427b);
  tmp526b = fmaT2(S8_S9, u[17], tmp526b);
  tmp625b = fmaT2(S8_S9, -u[22], tmp625b);
  tmp724b = fmaT2(S8_S9, u[21], tmp724b);
  tmp823b = fmaT2(S8_S9, u[30], tmp823b);
  tmp922b = fmaT2(S8_S9, -u[25], tmp922b);
  tmp1021b = fmaT2(S8_S9, u[24], tmp1021b);
  tmp1120b = fmaT2(S8_S9, u[19], tmp1120b);
  tmp1219b = fmaT2(S8_S9, u[20], tmp1219b);
  tmp1318b = fmaT2(S8_S9, u[28], tmp1318b);
  tmp1417b = fmaT2(S8_S9, u[26], tmp1417b);
  tmp1516b = fmaT2(S8_S9, u[16], tmp1516b);

  fma_addsub(u[1], u[30], SIN9, tmp130a, tmp130b);
  fma_addsub(u[2], u[29], SIN9, tmp229a, tmp229b);
  fma_addsub(u[3], u[28], SIN9, tmp328a, tmp328b);
  fma_addsub(u[4], u[27], SIN9, tmp427a, tmp427b);
  fma_addsub(u[5], u[26], SIN9, tmp526a, tmp526b);
  fma_addsub(u[6], u[25], SIN9, tmp625a, tmp625b);
  fma_addsub(u[7], u[24], SIN9, tmp724a, tmp724b);
  fma_addsub(u[8], u[23], SIN9, tmp823a, tmp823b);
  fma_addsub(u[9], u[22], SIN9, tmp922a, tmp922b);
  fma_addsub(u[10], u[21], SIN9, tmp1021a, tmp1021b);
  fma_addsub(u[11], u[20], SIN9, tmp1120a, tmp1120b);
  fma_addsub(u[12], u[19], SIN9, tmp1219a, tmp1219b);
  fma_addsub(u[13], u[18], SIN9, tmp1318a, tmp1318b);
  fma_addsub(u[14], u[17], SIN9, tmp1417a, tmp1417b);
  fma_addsub(u[15], u[16], SIN9, tmp1516a, tmp1516b);
}

#endif


/**************************************************************************/
/*            Similar to above, but for an FFT based on FP32              */
/**************************************************************************/

#if FFT_FP32

// 902 FMA + 90 ADD.  The method of fft11.cl (prime95's gwnum/zr11.mac): pair inputs k and 31-k, then each
// output pair j, 31-j is (u0 + sum_k cos(j*k*tau/31) sum_k) +- sum_k sin(j*k*tau/31) dif_k, with w = e^(2*pi*i/31).  As in fft11.cl, cbase = u0 +
// cos(8*tau/31) * sum zeroes one term per cos row (cos(8*tau/31) = -0.051 is the cosine nearest the middle of the range: every coefficient
// stays below 1.03 in magnitude), and the sin rows are scaled by sin(9*tau/31), the most accurately represented of the 15 sines
// (relative rounding error 1.6e-17).
// Accumulation order as fft13.cl: cos terms by |coef| rank 1,0,3,2,5,4,7,6,9,8,11,10,13,12 (ranks from the smallest), sin terms smallest first.

void OVERLOAD fft31(F2 *u) {
  const float COS8 = -0.05064916883871271227875f; // cos(8*tau/31): subtracted from every cos coefficient, it zeroes one term per cos row
  const float C1_C8 = 1.030179110091207206217f;   // cos(1*tau/31) - cos(8*tau/31)
  const float C2_C8 = 0.9696069804589433414059f;  // cos(2*tau/31) - cos(8*tau/31)
  const float C3_C8 = 0.8714126100459890386423f;  // cos(3*tau/31) - cos(8*tau/31)
  const float C4_C8 = 0.7396160879143992800796f;  // cos(4*tau/31) - cos(8*tau/31)
  const float C5_C8 = 0.5796131791656751696442f;  // cos(5*tau/31) - cos(8*tau/31)
  const float C6_C8 = 0.3979544216835329978206f;  // cos(6*tau/31) - cos(8*tau/31)
  const float C7_C8 = 0.2020769463432893759362f;  // cos(7*tau/31) - cos(8*tau/31)
  const float C9_C8 = -0.2000033634200078270361f; // cos(9*tau/31) - cos(8*tau/31)
  const float C10_C8 = -0.3897449827189215972374f;// cos(10*tau/31) - cos(8*tau/31)
  const float C11_C8 = -0.5614568137089501318680f;// cos(11*tau/31) - cos(8*tau/31)
  const float C12_C8 = -0.7081089538540781896345f;// cos(12*tau/31) - cos(8*tau/31)
  const float C13_C8 = -0.8236974473058694065487f;// cos(13*tau/31) - cos(8*tau/31)
  const float C14_C8 = -0.9034900875613361391971f;// cos(14*tau/31) - cos(8*tau/31)
  const float C15_C8 = -0.9442201545531824340426f;// cos(15*tau/31) - cos(8*tau/31)
  const float SIN9 = 0.9680771188662043051530f;   // sin(9*tau/31): factored out of the sin rows, it makes one term per row a plain add
  const float S1_S9 = 0.2079364506873352525185f;  // sin(1*tau/31) / sin(9*tau/31)
  const float S2_S9 = 0.4073599586520354364685f;  // sin(2*tau/31) / sin(9*tau/31)
  const float S3_S9 = 0.5901061020467584619302f;  // sin(3*tau/31) / sin(9*tau/31)
  const float S4_S9 = 0.7486932322891642376550f;  // sin(4*tau/31) / sin(9*tau/31)
  const float S5_S9 = 0.8766287736339320569721f;  // sin(5*tau/31) / sin(9*tau/31)
  const float S6_S9 = 0.9686750299866192862133f;  // sin(6*tau/31) / sin(9*tau/31)
  const float S7_S9 = 1.021063616797171006429f;   // sin(7*tau/31) / sin(9*tau/31)
  const float S8_S9 = 1.031649739166165645980f;   // sin(8*tau/31) / sin(9*tau/31)
  const float S10_S9 = 0.9274101433388233418963f; // sin(10*tau/31) / sin(9*tau/31)
  const float S11_S9 = 0.8168520064432902518249f; // sin(11*tau/31) / sin(9*tau/31)
  const float S12_S9 = 0.6728518524279333665138f; // sin(12*tau/31) / sin(9*tau/31)
  const float S13_S9 = 0.5013050645174410786296f; // sin(13*tau/31) / sin(9*tau/31)
  const float S14_S9 = 0.3092347883647606770760f; // sin(14*tau/31) / sin(9*tau/31)
  const float S15_S9 = 0.1045044036428821091535f; // sin(15*tau/31) / sin(9*tau/31)

  X2_mul_t4(u[1], u[30]);
  X2_mul_t4(u[2], u[29]);
  X2_mul_t4(u[3], u[28]);
  X2_mul_t4(u[4], u[27]);
  X2_mul_t4(u[5], u[26]);
  X2_mul_t4(u[6], u[25]);
  X2_mul_t4(u[7], u[24]);
  X2_mul_t4(u[8], u[23]);
  X2_mul_t4(u[9], u[22]);
  X2_mul_t4(u[10], u[21]);
  X2_mul_t4(u[11], u[20]);
  X2_mul_t4(u[12], u[19]);
  X2_mul_t4(u[13], u[18]);
  X2_mul_t4(u[14], u[17]);
  X2_mul_t4(u[15], u[16]);

  F2 sum = u[1] + u[2] + u[3] + u[4] + u[5] + u[6] + u[7] + u[8] + u[9] + u[10] + u[11] + u[12] + u[13] + u[14] + u[15];
  F2 cbase = fmaT2(COS8, sum, u[0]);                    // u[0] + COS8 * sum: each cos row below then has one term fewer
  u[0] = u[0] + sum;

  F2 tmp130a = fmaT2(C7_C8, u[7], cbase);
  F2 tmp229a = fmaT2(C7_C8, u[12], cbase);
  F2 tmp328a = fmaT2(C7_C8, u[8], cbase);
  F2 tmp427a = fmaT2(C7_C8, u[6], cbase);
  F2 tmp526a = fmaT2(C7_C8, u[11], cbase);
  F2 tmp625a = fmaT2(C7_C8, u[4], cbase);
  F2 tmp724a = fmaT2(C7_C8, u[1], cbase);
  F2 tmp823a = fmaT2(C7_C8, u[3], cbase);
  F2 tmp922a = fmaT2(C7_C8, u[13], cbase);
  F2 tmp1021a = fmaT2(C7_C8, u[10], cbase);
  F2 tmp1120a = fmaT2(C7_C8, u[5], cbase);
  F2 tmp1219a = fmaT2(C7_C8, u[2], cbase);
  F2 tmp1318a = fmaT2(C7_C8, u[9], cbase);
  F2 tmp1417a = fmaT2(C7_C8, u[15], cbase);
  F2 tmp1516a = fmaT2(C7_C8, u[14], cbase);

  tmp130a = fmaT2(C9_C8, u[9], tmp130a);
  tmp229a = fmaT2(C9_C8, u[11], tmp229a);
  tmp328a = fmaT2(C9_C8, u[3], tmp328a);
  tmp427a = fmaT2(C9_C8, u[10], tmp427a);
  tmp526a = fmaT2(C9_C8, u[8], tmp526a);
  tmp625a = fmaT2(C9_C8, u[14], tmp625a);
  tmp724a = fmaT2(C9_C8, u[12], tmp724a);
  tmp823a = fmaT2(C9_C8, u[5], tmp823a);
  tmp922a = fmaT2(C9_C8, u[1], tmp922a);
  tmp1021a = fmaT2(C9_C8, u[4], tmp1021a);
  tmp1120a = fmaT2(C9_C8, u[2], tmp1120a);
  tmp1219a = fmaT2(C9_C8, u[7], tmp1219a);
  tmp1318a = fmaT2(C9_C8, u[15], tmp1318a);
  tmp1417a = fmaT2(C9_C8, u[6], tmp1417a);
  tmp1516a = fmaT2(C9_C8, u[13], tmp1516a);

  tmp130a = fmaT2(C6_C8, u[6], tmp130a);
  tmp229a = fmaT2(C6_C8, u[3], tmp229a);
  tmp328a = fmaT2(C6_C8, u[2], tmp328a);
  tmp427a = fmaT2(C6_C8, u[14], tmp427a);
  tmp526a = fmaT2(C6_C8, u[5], tmp526a);
  tmp625a = fmaT2(C6_C8, u[1], tmp625a);
  tmp724a = fmaT2(C6_C8, u[8], tmp724a);
  tmp823a = fmaT2(C6_C8, u[7], tmp823a);
  tmp922a = fmaT2(C6_C8, u[11], tmp922a);
  tmp1021a = fmaT2(C6_C8, u[13], tmp1021a);
  tmp1120a = fmaT2(C6_C8, u[9], tmp1120a);
  tmp1219a = fmaT2(C6_C8, u[15], tmp1219a);
  tmp1318a = fmaT2(C6_C8, u[10], tmp1318a);
  tmp1417a = fmaT2(C6_C8, u[4], tmp1417a);
  tmp1516a = fmaT2(C6_C8, u[12], tmp1516a);

  tmp130a = fmaT2(C10_C8, u[10], tmp130a);
  tmp229a = fmaT2(C10_C8, u[5], tmp229a);
  tmp328a = fmaT2(C10_C8, u[7], tmp328a);
  tmp427a = fmaT2(C10_C8, u[13], tmp427a);
  tmp526a = fmaT2(C10_C8, u[2], tmp526a);
  tmp625a = fmaT2(C10_C8, u[12], tmp625a);
  tmp724a = fmaT2(C10_C8, u[3], tmp724a);
  tmp823a = fmaT2(C10_C8, u[9], tmp823a);
  tmp922a = fmaT2(C10_C8, u[8], tmp922a);
  tmp1021a = fmaT2(C10_C8, u[1], tmp1021a);
  tmp1120a = fmaT2(C10_C8, u[15], tmp1120a);
  tmp1219a = fmaT2(C10_C8, u[6], tmp1219a);
  tmp1318a = fmaT2(C10_C8, u[4], tmp1318a);
  tmp1417a = fmaT2(C10_C8, u[14], tmp1417a);
  tmp1516a = fmaT2(C10_C8, u[11], tmp1516a);

  tmp130a = fmaT2(C5_C8, u[5], tmp130a);
  tmp229a = fmaT2(C5_C8, u[13], tmp229a);
  tmp328a = fmaT2(C5_C8, u[12], tmp328a);
  tmp427a = fmaT2(C5_C8, u[9], tmp427a);
  tmp526a = fmaT2(C5_C8, u[1], tmp526a);
  tmp625a = fmaT2(C5_C8, u[6], tmp625a);
  tmp724a = fmaT2(C5_C8, u[14], tmp724a);
  tmp823a = fmaT2(C5_C8, u[11], tmp823a);
  tmp922a = fmaT2(C5_C8, u[4], tmp922a);
  tmp1021a = fmaT2(C5_C8, u[15], tmp1021a);
  tmp1120a = fmaT2(C5_C8, u[8], tmp1120a);
  tmp1219a = fmaT2(C5_C8, u[3], tmp1219a);
  tmp1318a = fmaT2(C5_C8, u[2], tmp1318a);
  tmp1417a = fmaT2(C5_C8, u[7], tmp1417a);
  tmp1516a = fmaT2(C5_C8, u[10], tmp1516a);

  tmp130a = fmaT2(C11_C8, u[11], tmp130a);
  tmp229a = fmaT2(C11_C8, u[10], tmp229a);
  tmp328a = fmaT2(C11_C8, u[14], tmp328a);
  tmp427a = fmaT2(C11_C8, u[5], tmp427a);
  tmp526a = fmaT2(C11_C8, u[4], tmp526a);
  tmp625a = fmaT2(C11_C8, u[7], tmp625a);
  tmp724a = fmaT2(C11_C8, u[6], tmp724a);
  tmp823a = fmaT2(C11_C8, u[13], tmp823a);
  tmp922a = fmaT2(C11_C8, u[15], tmp922a);
  tmp1021a = fmaT2(C11_C8, u[2], tmp1021a);
  tmp1120a = fmaT2(C11_C8, u[1], tmp1120a);
  tmp1219a = fmaT2(C11_C8, u[12], tmp1219a);
  tmp1318a = fmaT2(C11_C8, u[8], tmp1318a);
  tmp1417a = fmaT2(C11_C8, u[3], tmp1417a);
  tmp1516a = fmaT2(C11_C8, u[9], tmp1516a);

  tmp130a = fmaT2(C4_C8, u[4], tmp130a);
  tmp229a = fmaT2(C4_C8, u[2], tmp229a);
  tmp328a = fmaT2(C4_C8, u[9], tmp328a);
  tmp427a = fmaT2(C4_C8, u[1], tmp427a);
  tmp526a = fmaT2(C4_C8, u[7], tmp526a);
  tmp625a = fmaT2(C4_C8, u[11], tmp625a);
  tmp724a = fmaT2(C4_C8, u[5], tmp724a);
  tmp823a = fmaT2(C4_C8, u[15], tmp823a);
  tmp922a = fmaT2(C4_C8, u[3], tmp922a);
  tmp1021a = fmaT2(C4_C8, u[12], tmp1021a);
  tmp1120a = fmaT2(C4_C8, u[6], tmp1120a);
  tmp1219a = fmaT2(C4_C8, u[10], tmp1219a);
  tmp1318a = fmaT2(C4_C8, u[14], tmp1318a);
  tmp1417a = fmaT2(C4_C8, u[13], tmp1417a);
  tmp1516a = fmaT2(C4_C8, u[8], tmp1516a);

  tmp130a = fmaT2(C12_C8, u[12], tmp130a);
  tmp229a = fmaT2(C12_C8, u[6], tmp229a);
  tmp328a = fmaT2(C12_C8, u[4], tmp328a);
  tmp427a = fmaT2(C12_C8, u[3], tmp427a);
  tmp526a = fmaT2(C12_C8, u[10], tmp526a);
  tmp625a = fmaT2(C12_C8, u[2], tmp625a);
  tmp724a = fmaT2(C12_C8, u[15], tmp724a);
  tmp823a = fmaT2(C12_C8, u[14], tmp823a);
  tmp922a = fmaT2(C12_C8, u[9], tmp922a);
  tmp1021a = fmaT2(C12_C8, u[5], tmp1021a);
  tmp1120a = fmaT2(C12_C8, u[13], tmp1120a);
  tmp1219a = fmaT2(C12_C8, u[1], tmp1219a);
  tmp1318a = fmaT2(C12_C8, u[11], tmp1318a);
  tmp1417a = fmaT2(C12_C8, u[8], tmp1417a);
  tmp1516a = fmaT2(C12_C8, u[7], tmp1516a);

  tmp130a = fmaT2(C3_C8, u[3], tmp130a);
  tmp229a = fmaT2(C3_C8, u[14], tmp229a);
  tmp328a = fmaT2(C3_C8, u[1], tmp328a);
  tmp427a = fmaT2(C3_C8, u[7], tmp427a);
  tmp526a = fmaT2(C3_C8, u[13], tmp526a);
  tmp625a = fmaT2(C3_C8, u[15], tmp625a);
  tmp724a = fmaT2(C3_C8, u[4], tmp724a);
  tmp823a = fmaT2(C3_C8, u[12], tmp823a);
  tmp922a = fmaT2(C3_C8, u[10], tmp922a);
  tmp1021a = fmaT2(C3_C8, u[9], tmp1021a);
  tmp1120a = fmaT2(C3_C8, u[11], tmp1120a);
  tmp1219a = fmaT2(C3_C8, u[8], tmp1219a);
  tmp1318a = fmaT2(C3_C8, u[5], tmp1318a);
  tmp1417a = fmaT2(C3_C8, u[2], tmp1417a);
  tmp1516a = fmaT2(C3_C8, u[6], tmp1516a);

  tmp130a = fmaT2(C13_C8, u[13], tmp130a);
  tmp229a = fmaT2(C13_C8, u[9], tmp229a);
  tmp328a = fmaT2(C13_C8, u[6], tmp328a);
  tmp427a = fmaT2(C13_C8, u[11], tmp427a);
  tmp526a = fmaT2(C13_C8, u[15], tmp526a);
  tmp625a = fmaT2(C13_C8, u[3], tmp625a);
  tmp724a = fmaT2(C13_C8, u[7], tmp724a);
  tmp823a = fmaT2(C13_C8, u[10], tmp823a);
  tmp922a = fmaT2(C13_C8, u[2], tmp922a);
  tmp1021a = fmaT2(C13_C8, u[8], tmp1021a);
  tmp1120a = fmaT2(C13_C8, u[4], tmp1120a);
  tmp1219a = fmaT2(C13_C8, u[14], tmp1219a);
  tmp1318a = fmaT2(C13_C8, u[1], tmp1318a);
  tmp1417a = fmaT2(C13_C8, u[12], tmp1417a);
  tmp1516a = fmaT2(C13_C8, u[5], tmp1516a);

  tmp130a = fmaT2(C15_C8, u[15], tmp130a);
  tmp229a = fmaT2(C15_C8, u[8], tmp229a);
  tmp328a = fmaT2(C15_C8, u[5], tmp328a);
  tmp427a = fmaT2(C15_C8, u[4], tmp427a);
  tmp526a = fmaT2(C15_C8, u[3], tmp526a);
  tmp625a = fmaT2(C15_C8, u[13], tmp625a);
  tmp724a = fmaT2(C15_C8, u[11], tmp724a);
  tmp823a = fmaT2(C15_C8, u[2], tmp823a);
  tmp922a = fmaT2(C15_C8, u[12], tmp922a);
  tmp1021a = fmaT2(C15_C8, u[14], tmp1021a);
  tmp1120a = fmaT2(C15_C8, u[7], tmp1120a);
  tmp1219a = fmaT2(C15_C8, u[9], tmp1219a);
  tmp1318a = fmaT2(C15_C8, u[6], tmp1318a);
  tmp1417a = fmaT2(C15_C8, u[10], tmp1417a);
  tmp1516a = fmaT2(C15_C8, u[1], tmp1516a);

  tmp130a = fmaT2(C14_C8, u[14], tmp130a);
  tmp229a = fmaT2(C14_C8, u[7], tmp229a);
  tmp328a = fmaT2(C14_C8, u[15], tmp328a);
  tmp427a = fmaT2(C14_C8, u[12], tmp427a);
  tmp526a = fmaT2(C14_C8, u[9], tmp526a);
  tmp625a = fmaT2(C14_C8, u[8], tmp625a);
  tmp724a = fmaT2(C14_C8, u[2], tmp724a);
  tmp823a = fmaT2(C14_C8, u[6], tmp823a);
  tmp922a = fmaT2(C14_C8, u[5], tmp922a);
  tmp1021a = fmaT2(C14_C8, u[11], tmp1021a);
  tmp1120a = fmaT2(C14_C8, u[10], tmp1120a);
  tmp1219a = fmaT2(C14_C8, u[4], tmp1219a);
  tmp1318a = fmaT2(C14_C8, u[13], tmp1318a);
  tmp1417a = fmaT2(C14_C8, u[1], tmp1417a);
  tmp1516a = fmaT2(C14_C8, u[3], tmp1516a);

  tmp130a = fmaT2(C1_C8, u[1], tmp130a);
  tmp229a = fmaT2(C1_C8, u[15], tmp229a);
  tmp328a = fmaT2(C1_C8, u[10], tmp328a);
  tmp427a = fmaT2(C1_C8, u[8], tmp427a);
  tmp526a = fmaT2(C1_C8, u[6], tmp526a);
  tmp625a = fmaT2(C1_C8, u[5], tmp625a);
  tmp724a = fmaT2(C1_C8, u[9], tmp724a);
  tmp823a = fmaT2(C1_C8, u[4], tmp823a);
  tmp922a = fmaT2(C1_C8, u[7], tmp922a);
  tmp1021a = fmaT2(C1_C8, u[3], tmp1021a);
  tmp1120a = fmaT2(C1_C8, u[14], tmp1120a);
  tmp1219a = fmaT2(C1_C8, u[13], tmp1219a);
  tmp1318a = fmaT2(C1_C8, u[12], tmp1318a);
  tmp1417a = fmaT2(C1_C8, u[11], tmp1417a);
  tmp1516a = fmaT2(C1_C8, u[2], tmp1516a);

  tmp130a = fmaT2(C2_C8, u[2], tmp130a);
  tmp229a = fmaT2(C2_C8, u[1], tmp229a);
  tmp328a = fmaT2(C2_C8, u[11], tmp328a);
  tmp427a = fmaT2(C2_C8, u[15], tmp427a);
  tmp526a = fmaT2(C2_C8, u[12], tmp526a);
  tmp625a = fmaT2(C2_C8, u[10], tmp625a);
  tmp724a = fmaT2(C2_C8, u[13], tmp724a);
  tmp823a = fmaT2(C2_C8, u[8], tmp823a);
  tmp922a = fmaT2(C2_C8, u[14], tmp922a);
  tmp1021a = fmaT2(C2_C8, u[6], tmp1021a);
  tmp1120a = fmaT2(C2_C8, u[3], tmp1120a);
  tmp1219a = fmaT2(C2_C8, u[5], tmp1219a);
  tmp1318a = fmaT2(C2_C8, u[7], tmp1318a);
  tmp1417a = fmaT2(C2_C8, u[9], tmp1417a);
  tmp1516a = fmaT2(C2_C8, u[4], tmp1516a);

  F2 tmp130b = fmaT2(S15_S9, u[16], u[22]);
  F2 tmp229b = fmaT2(S15_S9, -u[23], -u[20]);
  F2 tmp328b = fmaT2(S15_S9, u[26], u[28]);
  F2 tmp427b = fmaT2(S15_S9, -u[27], u[21]);
  F2 tmp526b = fmaT2(S15_S9, u[28], u[23]);
  F2 tmp625b = fmaT2(S15_S9, -u[18], -u[17]);
  F2 tmp724b = fmaT2(S15_S9, u[20], -u[19]);
  F2 tmp823b = fmaT2(S15_S9, -u[29], u[26]);
  F2 tmp922b = fmaT2(S15_S9, u[19], u[30]);
  F2 tmp1021b = fmaT2(S15_S9, -u[17], u[27]);
  F2 tmp1120b = fmaT2(S15_S9, u[24], -u[29]);
  F2 tmp1219b = fmaT2(S15_S9, u[22], -u[24]);
  F2 tmp1318b = fmaT2(S15_S9, -u[25], u[16]);
  F2 tmp1417b = fmaT2(S15_S9, -u[21], -u[25]);
  F2 tmp1516b = fmaT2(S15_S9, u[30], u[18]);

  tmp130b = fmaT2(S1_S9, u[30], tmp130b);
  tmp229b = fmaT2(S1_S9, -u[16], tmp229b);
  tmp328b = fmaT2(S1_S9, -u[21], tmp328b);
  tmp427b = fmaT2(S1_S9, u[23], tmp427b);
  tmp526b = fmaT2(S1_S9, -u[25], tmp526b);
  tmp625b = fmaT2(S1_S9, -u[26], tmp625b);
  tmp724b = fmaT2(S1_S9, u[22], tmp724b);
  tmp823b = fmaT2(S1_S9, u[27], tmp823b);
  tmp922b = fmaT2(S1_S9, u[24], tmp922b);
  tmp1021b = fmaT2(S1_S9, -u[28], tmp1021b);
  tmp1120b = fmaT2(S1_S9, -u[17], tmp1120b);
  tmp1219b = fmaT2(S1_S9, u[18], tmp1219b);
  tmp1318b = fmaT2(S1_S9, u[19], tmp1318b);
  tmp1417b = fmaT2(S1_S9, -u[20], tmp1417b);
  tmp1516b = fmaT2(S1_S9, -u[29], tmp1516b);

  tmp130b = fmaT2(S14_S9, u[17], tmp130b);
  tmp229b = fmaT2(S14_S9, u[24], tmp229b);
  tmp328b = fmaT2(S14_S9, u[16], tmp328b);
  tmp427b = fmaT2(S14_S9, -u[19], tmp427b);
  tmp526b = fmaT2(S14_S9, u[22], tmp526b);
  tmp625b = fmaT2(S14_S9, -u[23], tmp625b);
  tmp724b = fmaT2(S14_S9, u[29], tmp724b);
  tmp823b = fmaT2(S14_S9, -u[25], tmp823b);
  tmp922b = fmaT2(S14_S9, u[26], tmp922b);
  tmp1021b = fmaT2(S14_S9, -u[20], tmp1021b);
  tmp1120b = fmaT2(S14_S9, -u[21], tmp1120b);
  tmp1219b = fmaT2(S14_S9, -u[27], tmp1219b);
  tmp1318b = fmaT2(S14_S9, u[18], tmp1318b);
  tmp1417b = fmaT2(S14_S9, u[30], tmp1417b);
  tmp1516b = fmaT2(S14_S9, u[28], tmp1516b);

  tmp130b = fmaT2(S2_S9, u[29], tmp130b);
  tmp229b = fmaT2(S2_S9, u[30], tmp229b);
  tmp328b = fmaT2(S2_S9, u[20], tmp328b);
  tmp427b = fmaT2(S2_S9, -u[16], tmp427b);
  tmp526b = fmaT2(S2_S9, -u[19], tmp526b);
  tmp625b = fmaT2(S2_S9, -u[21], tmp625b);
  tmp724b = fmaT2(S2_S9, -u[18], tmp724b);
  tmp823b = fmaT2(S2_S9, u[23], tmp823b);
  tmp922b = fmaT2(S2_S9, u[17], tmp922b);
  tmp1021b = fmaT2(S2_S9, -u[25], tmp1021b);
  tmp1120b = fmaT2(S2_S9, u[28], tmp1120b);
  tmp1219b = fmaT2(S2_S9, -u[26], tmp1219b);
  tmp1318b = fmaT2(S2_S9, -u[24], tmp1318b);
  tmp1417b = fmaT2(S2_S9, u[22], tmp1417b);
  tmp1516b = fmaT2(S2_S9, -u[27], tmp1516b);

  tmp130b = fmaT2(S13_S9, u[18], tmp130b);
  tmp229b = fmaT2(S13_S9, -u[22], tmp229b);
  tmp328b = fmaT2(S13_S9, -u[25], tmp328b);
  tmp427b = fmaT2(S13_S9, u[20], tmp427b);
  tmp526b = fmaT2(S13_S9, u[16], tmp526b);
  tmp625b = fmaT2(S13_S9, -u[28], tmp625b);
  tmp724b = fmaT2(S13_S9, -u[24], tmp724b);
  tmp823b = fmaT2(S13_S9, -u[21], tmp823b);
  tmp922b = fmaT2(S13_S9, -u[29], tmp922b);
  tmp1021b = fmaT2(S13_S9, -u[23], tmp1021b);
  tmp1120b = fmaT2(S13_S9, u[27], tmp1120b);
  tmp1219b = fmaT2(S13_S9, u[17], tmp1219b);
  tmp1318b = fmaT2(S13_S9, u[30], tmp1318b);
  tmp1417b = fmaT2(S13_S9, u[19], tmp1417b);
  tmp1516b = fmaT2(S13_S9, u[26], tmp1516b);

  tmp130b = fmaT2(S3_S9, u[28], tmp130b);
  tmp229b = fmaT2(S3_S9, -u[17], tmp229b);
  tmp328b = fmaT2(S3_S9, u[30], tmp328b);
  tmp427b = fmaT2(S3_S9, -u[24], tmp427b);
  tmp526b = fmaT2(S3_S9, u[18], tmp526b);
  tmp625b = fmaT2(S3_S9, -u[16], tmp625b);
  tmp724b = fmaT2(S3_S9, -u[27], tmp724b);
  tmp823b = fmaT2(S3_S9, u[19], tmp823b);
  tmp922b = fmaT2(S3_S9, -u[21], tmp922b);
  tmp1021b = fmaT2(S3_S9, -u[22], tmp1021b);
  tmp1120b = fmaT2(S3_S9, -u[20], tmp1120b);
  tmp1219b = fmaT2(S3_S9, u[23], tmp1219b);
  tmp1318b = fmaT2(S3_S9, u[26], tmp1318b);
  tmp1417b = fmaT2(S3_S9, -u[29], tmp1417b);
  tmp1516b = fmaT2(S3_S9, -u[25], tmp1516b);

  tmp130b = fmaT2(S12_S9, u[19], tmp130b);
  tmp229b = fmaT2(S12_S9, u[25], tmp229b);
  tmp328b = fmaT2(S12_S9, u[27], tmp328b);
  tmp427b = fmaT2(S12_S9, u[28], tmp427b);
  tmp526b = fmaT2(S12_S9, -u[21], tmp526b);
  tmp625b = fmaT2(S12_S9, u[29], tmp625b);
  tmp724b = fmaT2(S12_S9, u[16], tmp724b);
  tmp823b = fmaT2(S12_S9, -u[17], tmp823b);
  tmp922b = fmaT2(S12_S9, -u[22], tmp922b);
  tmp1021b = fmaT2(S12_S9, -u[26], tmp1021b);
  tmp1120b = fmaT2(S12_S9, -u[18], tmp1120b);
  tmp1219b = fmaT2(S12_S9, u[30], tmp1219b);
  tmp1318b = fmaT2(S12_S9, -u[20], tmp1318b);
  tmp1417b = fmaT2(S12_S9, -u[23], tmp1417b);
  tmp1516b = fmaT2(S12_S9, u[24], tmp1516b);

  tmp130b = fmaT2(S4_S9, u[27], tmp130b);
  tmp229b = fmaT2(S4_S9, u[29], tmp229b);
  tmp328b = fmaT2(S4_S9, -u[22], tmp328b);
  tmp427b = fmaT2(S4_S9, u[30], tmp427b);
  tmp526b = fmaT2(S4_S9, u[24], tmp526b);
  tmp625b = fmaT2(S4_S9, u[20], tmp625b);
  tmp724b = fmaT2(S4_S9, u[26], tmp724b);
  tmp823b = fmaT2(S4_S9, -u[16], tmp823b);
  tmp922b = fmaT2(S4_S9, -u[28], tmp922b);
  tmp1021b = fmaT2(S4_S9, -u[19], tmp1021b);
  tmp1120b = fmaT2(S4_S9, u[25], tmp1120b);
  tmp1219b = fmaT2(S4_S9, -u[21], tmp1219b);
  tmp1318b = fmaT2(S4_S9, -u[17], tmp1318b);
  tmp1417b = fmaT2(S4_S9, -u[18], tmp1417b);
  tmp1516b = fmaT2(S4_S9, -u[23], tmp1516b);

  tmp130b = fmaT2(S11_S9, u[20], tmp130b);
  tmp229b = fmaT2(S11_S9, -u[21], tmp229b);
  tmp328b = fmaT2(S11_S9, u[17], tmp328b);
  tmp427b = fmaT2(S11_S9, -u[26], tmp427b);
  tmp526b = fmaT2(S11_S9, -u[27], tmp526b);
  tmp625b = fmaT2(S11_S9, u[24], tmp625b);
  tmp724b = fmaT2(S11_S9, u[25], tmp724b);
  tmp823b = fmaT2(S11_S9, u[18], tmp823b);
  tmp922b = fmaT2(S11_S9, u[16], tmp922b);
  tmp1021b = fmaT2(S11_S9, -u[29], tmp1021b);
  tmp1120b = fmaT2(S11_S9, u[30], tmp1120b);
  tmp1219b = fmaT2(S11_S9, -u[19], tmp1219b);
  tmp1318b = fmaT2(S11_S9, u[23], tmp1318b);
  tmp1417b = fmaT2(S11_S9, u[28], tmp1417b);
  tmp1516b = fmaT2(S11_S9, u[22], tmp1516b);

  tmp130b = fmaT2(S5_S9, u[26], tmp130b);
  tmp229b = fmaT2(S5_S9, -u[18], tmp229b);
  tmp328b = fmaT2(S5_S9, u[19], tmp328b);
  tmp427b = fmaT2(S5_S9, u[22], tmp427b);
  tmp526b = fmaT2(S5_S9, u[30], tmp526b);
  tmp625b = fmaT2(S5_S9, u[25], tmp625b);
  tmp724b = fmaT2(S5_S9, u[17], tmp724b);
  tmp823b = fmaT2(S5_S9, -u[20], tmp823b);
  tmp922b = fmaT2(S5_S9, u[27], tmp922b);
  tmp1021b = fmaT2(S5_S9, -u[16], tmp1021b);
  tmp1120b = fmaT2(S5_S9, -u[23], tmp1120b);
  tmp1219b = fmaT2(S5_S9, u[28], tmp1219b);
  tmp1318b = fmaT2(S5_S9, -u[29], tmp1318b);
  tmp1417b = fmaT2(S5_S9, u[24], tmp1417b);
  tmp1516b = fmaT2(S5_S9, -u[21], tmp1516b);

  tmp130b = fmaT2(S10_S9, u[21], tmp130b);
  tmp229b = fmaT2(S10_S9, u[26], tmp229b);
  tmp328b = fmaT2(S10_S9, -u[24], tmp328b);
  tmp427b = fmaT2(S10_S9, -u[18], tmp427b);
  tmp526b = fmaT2(S10_S9, u[29], tmp526b);
  tmp625b = fmaT2(S10_S9, u[19], tmp625b);
  tmp724b = fmaT2(S10_S9, -u[28], tmp724b);
  tmp823b = fmaT2(S10_S9, u[22], tmp823b);
  tmp922b = fmaT2(S10_S9, u[23], tmp922b);
  tmp1021b = fmaT2(S10_S9, u[30], tmp1021b);
  tmp1120b = fmaT2(S10_S9, u[16], tmp1120b);
  tmp1219b = fmaT2(S10_S9, u[25], tmp1219b);
  tmp1318b = fmaT2(S10_S9, -u[27], tmp1318b);
  tmp1417b = fmaT2(S10_S9, u[17], tmp1417b);
  tmp1516b = fmaT2(S10_S9, u[20], tmp1516b);

  tmp130b = fmaT2(S6_S9, u[25], tmp130b);
  tmp229b = fmaT2(S6_S9, u[28], tmp229b);
  tmp328b = fmaT2(S6_S9, u[29], tmp328b);
  tmp427b = fmaT2(S6_S9, -u[17], tmp427b);
  tmp526b = fmaT2(S6_S9, -u[26], tmp526b);
  tmp625b = fmaT2(S6_S9, u[30], tmp625b);
  tmp724b = fmaT2(S6_S9, -u[23], tmp724b);
  tmp823b = fmaT2(S6_S9, -u[24], tmp823b);
  tmp922b = fmaT2(S6_S9, u[20], tmp922b);
  tmp1021b = fmaT2(S6_S9, u[18], tmp1021b);
  tmp1120b = fmaT2(S6_S9, u[22], tmp1120b);
  tmp1219b = fmaT2(S6_S9, -u[16], tmp1219b);
  tmp1318b = fmaT2(S6_S9, u[21], tmp1318b);
  tmp1417b = fmaT2(S6_S9, -u[27], tmp1417b);
  tmp1516b = fmaT2(S6_S9, -u[19], tmp1516b);

  tmp130b = fmaT2(S7_S9, u[24], tmp130b);
  tmp229b = fmaT2(S7_S9, -u[19], tmp229b);
  tmp328b = fmaT2(S7_S9, -u[23], tmp328b);
  tmp427b = fmaT2(S7_S9, -u[25], tmp427b);
  tmp526b = fmaT2(S7_S9, -u[20], tmp526b);
  tmp625b = fmaT2(S7_S9, -u[27], tmp625b);
  tmp724b = fmaT2(S7_S9, u[30], tmp724b);
  tmp823b = fmaT2(S7_S9, -u[28], tmp823b);
  tmp922b = fmaT2(S7_S9, -u[18], tmp922b);
  tmp1021b = fmaT2(S7_S9, u[21], tmp1021b);
  tmp1120b = fmaT2(S7_S9, -u[26], tmp1120b);
  tmp1219b = fmaT2(S7_S9, -u[29], tmp1219b);
  tmp1318b = fmaT2(S7_S9, -u[22], tmp1318b);
  tmp1417b = fmaT2(S7_S9, -u[16], tmp1417b);
  tmp1516b = fmaT2(S7_S9, -u[17], tmp1516b);

  tmp130b = fmaT2(S8_S9, u[23], tmp130b);
  tmp229b = fmaT2(S8_S9, u[27], tmp229b);
  tmp328b = fmaT2(S8_S9, u[18], tmp328b);
  tmp427b = fmaT2(S8_S9, u[29], tmp427b);
  tmp526b = fmaT2(S8_S9, u[17], tmp526b);
  tmp625b = fmaT2(S8_S9, -u[22], tmp625b);
  tmp724b = fmaT2(S8_S9, u[21], tmp724b);
  tmp823b = fmaT2(S8_S9, u[30], tmp823b);
  tmp922b = fmaT2(S8_S9, -u[25], tmp922b);
  tmp1021b = fmaT2(S8_S9, u[24], tmp1021b);
  tmp1120b = fmaT2(S8_S9, u[19], tmp1120b);
  tmp1219b = fmaT2(S8_S9, u[20], tmp1219b);
  tmp1318b = fmaT2(S8_S9, u[28], tmp1318b);
  tmp1417b = fmaT2(S8_S9, u[26], tmp1417b);
  tmp1516b = fmaT2(S8_S9, u[16], tmp1516b);

  fma_addsub(u[1], u[30], SIN9, tmp130a, tmp130b);
  fma_addsub(u[2], u[29], SIN9, tmp229a, tmp229b);
  fma_addsub(u[3], u[28], SIN9, tmp328a, tmp328b);
  fma_addsub(u[4], u[27], SIN9, tmp427a, tmp427b);
  fma_addsub(u[5], u[26], SIN9, tmp526a, tmp526b);
  fma_addsub(u[6], u[25], SIN9, tmp625a, tmp625b);
  fma_addsub(u[7], u[24], SIN9, tmp724a, tmp724b);
  fma_addsub(u[8], u[23], SIN9, tmp823a, tmp823b);
  fma_addsub(u[9], u[22], SIN9, tmp922a, tmp922b);
  fma_addsub(u[10], u[21], SIN9, tmp1021a, tmp1021b);
  fma_addsub(u[11], u[20], SIN9, tmp1120a, tmp1120b);
  fma_addsub(u[12], u[19], SIN9, tmp1219a, tmp1219b);
  fma_addsub(u[13], u[18], SIN9, tmp1318a, tmp1318b);
  fma_addsub(u[14], u[17], SIN9, tmp1417a, tmp1417b);
  fma_addsub(u[15], u[16], SIN9, tmp1516a, tmp1516b);
}

#endif
