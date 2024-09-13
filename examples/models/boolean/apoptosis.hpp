RULE_SET(APOPTOSIS, 41, "TNF", "GF", "TNFR1", "TRADD", "TRAF", "RIP", "cIAP",
         "TNFR2", "TRAF2", "NIK", "IKK", "A20", "MEKK1", "JNKK", "GFR", "PI3K",
         "PIP2", "PIP3", "AKT", "Mdm2", "FADD", "Cas8", "IkB", "NFkB", "JNK",
         "Cas7", "Cas12", "Cas9", "APC", "Cas3", "Cas3p", "Cas6", "BID", "BclX",
         "BAD", "P53", "Apaf1", "PTEN", "Mito", "IAP", "DNAd")

template <> bool APOPTOSIS::rule<0>(bitset<41> x) {
  int A = x[0];
  int H = 1 - x[0];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[0];
}
template <> vector<ind> APOPTOSIS::depends_on<0>() { return {0}; }

template <> bool APOPTOSIS::rule<1>(bitset<41> x) {
  int A = x[1];
  int H = 1 - x[1];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[1];
}
template <> vector<ind> APOPTOSIS::depends_on<1>() { return {1}; }

template <> bool APOPTOSIS::rule<2>(bitset<41> x) {
  int A = x[0];
  int H = 1 - x[0];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[2];
}
template <> vector<ind> APOPTOSIS::depends_on<2>() { return {0, 2}; }

template <> bool APOPTOSIS::rule<3>(bitset<41> x) {
  int A = x[2];
  int H = 1 - x[2];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[3];
}
template <> vector<ind> APOPTOSIS::depends_on<3>() { return {2, 3}; }

template <> bool APOPTOSIS::rule<4>(bitset<41> x) {
  int A = x[3];
  int H = 1 - x[3];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[4];
}
template <> vector<ind> APOPTOSIS::depends_on<4>() { return {3, 4}; }

template <> bool APOPTOSIS::rule<5>(bitset<41> x) {
  int A = x[3];
  int H = 1 - x[3];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[5];
}
template <> vector<ind> APOPTOSIS::depends_on<5>() { return {3, 5}; }

template <> bool APOPTOSIS::rule<6>(bitset<41> x) {
  int A = x[4];
  int H = 1 - x[4];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[6];
}
template <> vector<ind> APOPTOSIS::depends_on<6>() { return {4, 6}; }

template <> bool APOPTOSIS::rule<7>(bitset<41> x) {
  int A = x[0];
  int H = 1 - x[0];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[7];
}
template <> vector<ind> APOPTOSIS::depends_on<7>() { return {0, 7}; }

template <> bool APOPTOSIS::rule<8>(bitset<41> x) {
  int A = x[5];
  int H = -x[5] + 2 * x[7] + 1;
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[8];
}
template <> vector<ind> APOPTOSIS::depends_on<8>() { return {5, 7, 8}; }

template <> bool APOPTOSIS::rule<9>(bitset<41> x) {
  int A = x[8];
  int H = 1 - x[8];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[9];
}
template <> vector<ind> APOPTOSIS::depends_on<9>() { return {8, 9}; }

template <> bool APOPTOSIS::rule<10>(bitset<41> x) {
  int A = x[18] + x[9];
  int H = 3 * x[11] + (1 - x[18]) * (1 - x[9]);
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[10];
}
template <> vector<ind> APOPTOSIS::depends_on<10>() { return {9, 10, 11, 18}; }

template <> bool APOPTOSIS::rule<11>(bitset<41> x) {
  int A = x[23];
  int H = 1 - x[23];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[11];
}
template <> vector<ind> APOPTOSIS::depends_on<11>() { return {11, 23}; }

template <> bool APOPTOSIS::rule<12>(bitset<41> x) {
  int A = x[4];
  int H = 1 - x[4];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[12];
}
template <> vector<ind> APOPTOSIS::depends_on<12>() { return {4, 12}; }

template <> bool APOPTOSIS::rule<13>(bitset<41> x) {
  int A = x[12];
  int H = 2 * x[18] - x[12] + 1;
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[13];
}
template <> vector<ind> APOPTOSIS::depends_on<13>() { return {12, 13, 18}; }

template <> bool APOPTOSIS::rule<14>(bitset<41> x) {
  int A = x[1];
  int H = 1 - x[1];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[14];
}
template <> vector<ind> APOPTOSIS::depends_on<14>() { return {1, 14}; }

template <> bool APOPTOSIS::rule<15>(bitset<41> x) {
  int A = x[14];
  int H = 1 - x[14];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[15];
}
template <> vector<ind> APOPTOSIS::depends_on<15>() { return {14, 15}; }

template <> bool APOPTOSIS::rule<16>(bitset<41> x) {
  int A = x[14];
  int H = 1 - x[14];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[16];
}
template <> vector<ind> APOPTOSIS::depends_on<16>() { return {14, 16}; }

template <> bool APOPTOSIS::rule<17>(bitset<41> x) {
  int A = x[15] * x[16];
  int H = -x[16] + 2 * x[37] + 1;
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[17];
}
template <> vector<ind> APOPTOSIS::depends_on<17>() { return {15, 16, 17, 37}; }

template <> bool APOPTOSIS::rule<18>(bitset<41> x) {
  int A = x[17];
  int H = 1 - x[17];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[18];
}
template <> vector<ind> APOPTOSIS::depends_on<18>() { return {17, 18}; }

template <> bool APOPTOSIS::rule<19>(bitset<41> x) {
  int A = x[18] + x[35];
  int H = 1 - x[18];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[19];
}
template <> vector<ind> APOPTOSIS::depends_on<19>() { return {18, 19, 35}; }

template <> bool APOPTOSIS::rule<20>(bitset<41> x) {
  int A = x[3];
  int H = 1 - x[3];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[20];
}
template <> vector<ind> APOPTOSIS::depends_on<20>() { return {3, 20}; }

template <> bool APOPTOSIS::rule<21>(bitset<41> x) {
  int A = x[31] + x[20];
  int H = x[6];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[21];
}
template <> vector<ind> APOPTOSIS::depends_on<21>() { return {6, 20, 21, 31}; }

template <> bool APOPTOSIS::rule<22>(bitset<41> x) {
  int A = x[23];
  int H = 2 * x[10];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[22];
}
template <> vector<ind> APOPTOSIS::depends_on<22>() { return {10, 22, 23}; }

template <> bool APOPTOSIS::rule<23>(bitset<41> x) {
  int A = 1 - x[22];
  int H = x[22];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[23];
}
template <> vector<ind> APOPTOSIS::depends_on<23>() { return {22, 23}; }

template <> bool APOPTOSIS::rule<24>(bitset<41> x) {
  int A = x[13];
  int H = 1 - x[13];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[24];
}
template <> vector<ind> APOPTOSIS::depends_on<24>() { return {13, 24}; }

template <> bool APOPTOSIS::rule<25>(bitset<41> x) {
  int A = x[28] + x[21];
  int H = 3 * x[39] + (1 - x[28]) * (1 - x[21]);
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[25];
}
template <> vector<ind> APOPTOSIS::depends_on<25>() { return {21, 25, 28, 39}; }

template <> bool APOPTOSIS::rule<26>(bitset<41> x) {
  int A = x[25];
  int H = 1 - x[25];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[26];
}
template <> vector<ind> APOPTOSIS::depends_on<26>() { return {25, 26}; }

template <> bool APOPTOSIS::rule<27>(bitset<41> x) {
  int A = x[26] + x[29];
  int H = x[18] - x[26] + 2 * x[39] + 1;
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[27];
}
template <> vector<ind> APOPTOSIS::depends_on<27>() {
  return {18, 26, 27, 29, 39};
}

template <> bool APOPTOSIS::rule<28>(bitset<41> x) {
  int A = x[36] * x[27] * x[38];
  int H = -x[36] - x[27] - x[38] + 2 * x[39] + 3;
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[28];
}
template <> vector<ind> APOPTOSIS::depends_on<28>() {
  return {27, 28, 36, 38, 39};
}

template <> bool APOPTOSIS::rule<29>(bitset<41> x) {
  int A = x[28] + x[31] + x[21];
  int H = 4 * x[39] + (1 - x[28]) * (1 - x[31]) * (1 - x[21]);
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[29];
}
template <> vector<ind> APOPTOSIS::depends_on<29>() {
  return {21, 28, 29, 31, 39};
}

template <> bool APOPTOSIS::rule<30>(bitset<41> x) { return x[29]; }
template <> vector<ind> APOPTOSIS::depends_on<30>() { return {29, 30}; }

template <> bool APOPTOSIS::rule<31>(bitset<41> x) {
  int A = x[29];
  int H = -x[29] + 2 * x[39] + 1;
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[31];
}
template <> vector<ind> APOPTOSIS::depends_on<31>() { return {29, 31, 39}; }

template <> bool APOPTOSIS::rule<32>(bitset<41> x) {
  int A = x[21] * x[35] + x[24] * x[35];
  int H = 3 * x[33] - x[35] + 1;
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[32];
}
template <> vector<ind> APOPTOSIS::depends_on<32>() {
  return {21, 24, 32, 33, 35};
}

template <> bool APOPTOSIS::rule<33>(bitset<41> x) {
  int A = x[23];
  int H = 2 * x[34] - x[23] + x[35] + 1;
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[33];
}
template <> vector<ind> APOPTOSIS::depends_on<33>() { return {23, 33, 34, 35}; }

template <> bool APOPTOSIS::rule<34>(bitset<41> x) {
  int A = x[35];
  int H = 2 * x[18] - x[35] + 1;
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[34];
}
template <> vector<ind> APOPTOSIS::depends_on<34>() { return {18, 34, 35}; }

template <> bool APOPTOSIS::rule<35>(bitset<41> x) {
  int A = 3 * x[40] + x[24];
  int H = 2 * x[19] + (1 - x[40]) * (1 - x[24]);
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[35];
}
template <> vector<ind> APOPTOSIS::depends_on<35>() { return {19, 24, 35, 40}; }

template <> bool APOPTOSIS::rule<36>(bitset<41> x) {
  int A = x[35];
  int H = 1 - x[35];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[36];
}
template <> vector<ind> APOPTOSIS::depends_on<36>() { return {35, 36}; }

template <> bool APOPTOSIS::rule<37>(bitset<41> x) {
  int A = x[35];
  int H = 1 - x[35];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[37];
}
template <> vector<ind> APOPTOSIS::depends_on<37>() { return {35, 37}; }

template <> bool APOPTOSIS::rule<38>(bitset<41> x) {
  int A = x[32];
  int H = x[33];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[38];
}
template <> vector<ind> APOPTOSIS::depends_on<38>() { return {32, 33, 38}; }

template <> bool APOPTOSIS::rule<39>(bitset<41> x) {
  int A = x[23];
  int H = x[31] * x[29] + x[38];
  if (A > H)
    return true;
  else if (H > A)
    return false;
  else
    return x[39];
}
template <> vector<ind> APOPTOSIS::depends_on<39>() {
  return {23, 29, 31, 38, 39};
}

template <> bool APOPTOSIS::rule<40>(bitset<41> x) { return x[29] && x[30]; }
template <> vector<ind> APOPTOSIS::depends_on<40>() { return {29, 30, 40}; }
