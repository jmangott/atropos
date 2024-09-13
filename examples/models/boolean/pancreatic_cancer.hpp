RULE_SET(PANCREATICCANCER, 34, "HMGB1", "TLR24", "RAGE", "MYD88", "RAS", "RAC1",
         "IRAKs", "RAF", "MEK", "PI3K", "ERK", "AP1", "TAB1", "PIP3", "AKT",
         "Myc", "INK4a", "IKK", "CyclD", "PTEN", "MDM2", "A20", "E2F", "IkB",
         "RB", "P53", "NFkB", "ARF", "P21", "BAX", "BclXL", "CyclE", "Apop",
         "Prol")

template <> bool PANCREATICCANCER::rule<0>(bitset<34> x) { return x[0]; }
template <> vector<ind> PANCREATICCANCER::depends_on<0>() { return {0}; }

template <> bool PANCREATICCANCER::rule<1>(bitset<34> x) {
  return x[0] || x[1];
}
template <> vector<ind> PANCREATICCANCER::depends_on<1>() { return {0, 1}; }

template <> bool PANCREATICCANCER::rule<2>(bitset<34> x) {
  return x[0] || x[2];
}
template <> vector<ind> PANCREATICCANCER::depends_on<2>() { return {0, 2}; }

template <> bool PANCREATICCANCER::rule<3>(bitset<34> x) {
  return x[3] || x[1];
}
template <> vector<ind> PANCREATICCANCER::depends_on<3>() { return {1, 3}; }

template <> bool PANCREATICCANCER::rule<4>(bitset<34> x) {
  return x[2] || x[4];
}
template <> vector<ind> PANCREATICCANCER::depends_on<4>() { return {2, 4}; }

template <> bool PANCREATICCANCER::rule<5>(bitset<34> x) {
  return x[3] || x[5];
}
template <> vector<ind> PANCREATICCANCER::depends_on<5>() { return {3, 5}; }

template <> bool PANCREATICCANCER::rule<6>(bitset<34> x) {
  return x[6] || x[3];
}
template <> vector<ind> PANCREATICCANCER::depends_on<6>() { return {3, 6}; }

template <> bool PANCREATICCANCER::rule<7>(bitset<34> x) {
  return x[14] || x[7] || x[4];
}
template <> vector<ind> PANCREATICCANCER::depends_on<7>() { return {4, 7, 14}; }

template <> bool PANCREATICCANCER::rule<8>(bitset<34> x) {
  return x[8] || x[7];
}
template <> vector<ind> PANCREATICCANCER::depends_on<8>() { return {7, 8}; }

template <> bool PANCREATICCANCER::rule<9>(bitset<34> x) {
  return x[9] || x[5] || x[4];
}
template <> vector<ind> PANCREATICCANCER::depends_on<9>() { return {4, 5, 9}; }

template <> bool PANCREATICCANCER::rule<10>(bitset<34> x) {
  return x[10] || x[6] || x[8];
}
template <> vector<ind> PANCREATICCANCER::depends_on<10>() {
  return {6, 8, 10};
}

template <> bool PANCREATICCANCER::rule<11>(bitset<34> x) {
  return x[11] || x[10];
}
template <> vector<ind> PANCREATICCANCER::depends_on<11>() { return {10, 11}; }

template <> bool PANCREATICCANCER::rule<12>(bitset<34> x) {
  return x[6] || x[12];
}
template <> vector<ind> PANCREATICCANCER::depends_on<12>() { return {6, 12}; }

template <> bool PANCREATICCANCER::rule<13>(bitset<34> x) {
  return !x[19] && (x[9] || x[13]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<13>() {
  return {9, 13, 19};
}

template <> bool PANCREATICCANCER::rule<14>(bitset<34> x) {
  return x[14] || x[13];
}
template <> vector<ind> PANCREATICCANCER::depends_on<14>() { return {13, 14}; }

template <> bool PANCREATICCANCER::rule<15>(bitset<34> x) {
  return x[10] || x[15] || x[26];
}
template <> vector<ind> PANCREATICCANCER::depends_on<15>() {
  return {10, 15, 26};
}

template <> bool PANCREATICCANCER::rule<16>(bitset<34> x) { return x[16]; }
template <> vector<ind> PANCREATICCANCER::depends_on<16>() { return {16}; }

template <> bool PANCREATICCANCER::rule<17>(bitset<34> x) {
  return !x[21] && (x[14] || x[10] || x[17] || x[12]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<17>() {
  return {10, 12, 14, 17, 21};
}

template <> bool PANCREATICCANCER::rule<18>(bitset<34> x) {
  return !(x[16] || x[28]) && (x[11] || x[18] || x[15] || x[26]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<18>() {
  return {11, 15, 16, 18, 26, 28};
}

template <> bool PANCREATICCANCER::rule<19>(bitset<34> x) {
  return x[25] || x[19];
}
template <> vector<ind> PANCREATICCANCER::depends_on<19>() { return {19, 25}; }

template <> bool PANCREATICCANCER::rule<20>(bitset<34> x) {
  return !x[27] && (x[14] || x[20] || x[25]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<20>() {
  return {14, 20, 25, 27};
}

template <> bool PANCREATICCANCER::rule<21>(bitset<34> x) {
  return x[21] || x[26];
}
template <> vector<ind> PANCREATICCANCER::depends_on<21>() { return {21, 26}; }

template <> bool PANCREATICCANCER::rule<22>(bitset<34> x) {
  return !x[24] && (x[22] || x[15]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<22>() {
  return {15, 22, 24};
}

template <> bool PANCREATICCANCER::rule<23>(bitset<34> x) {
  return !x[17] && (x[23] || x[26]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<23>() {
  return {17, 23, 26};
}

template <> bool PANCREATICCANCER::rule<24>(bitset<34> x) {
  return x[24] && !(x[18] || x[31]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<24>() {
  return {18, 24, 31};
}

template <> bool PANCREATICCANCER::rule<25>(bitset<34> x) {
  return x[25] && !x[20];
}
template <> vector<ind> PANCREATICCANCER::depends_on<25>() { return {20, 25}; }

template <> bool PANCREATICCANCER::rule<26>(bitset<34> x) {
  return x[26] && !x[23];
}
template <> vector<ind> PANCREATICCANCER::depends_on<26>() { return {23, 26}; }

template <> bool PANCREATICCANCER::rule<27>(bitset<34> x) {
  return x[27] || x[22];
}
template <> vector<ind> PANCREATICCANCER::depends_on<27>() { return {22, 27}; }

template <> bool PANCREATICCANCER::rule<28>(bitset<34> x) {
  return x[28] || x[25];
}
template <> vector<ind> PANCREATICCANCER::depends_on<28>() { return {25, 28}; }

template <> bool PANCREATICCANCER::rule<29>(bitset<34> x) {
  return x[29] || x[25];
}
template <> vector<ind> PANCREATICCANCER::depends_on<29>() { return {25, 29}; }

template <> bool PANCREATICCANCER::rule<30>(bitset<34> x) {
  return !x[25] && (x[30] || x[26]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<30>() {
  return {25, 26, 30};
}

template <> bool PANCREATICCANCER::rule<31>(bitset<34> x) {
  return !x[28] && (x[31] || x[22]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<31>() {
  return {22, 28, 31};
}

template <> bool PANCREATICCANCER::rule<32>(bitset<34> x) {
  return !x[30] && (x[32] || x[29]);
}
template <> vector<ind> PANCREATICCANCER::depends_on<32>() {
  return {29, 30, 32};
}

template <> bool PANCREATICCANCER::rule<33>(bitset<34> x) {
  return x[31] && x[33];
}
template <> vector<ind> PANCREATICCANCER::depends_on<33>() { return {31, 33}; }
