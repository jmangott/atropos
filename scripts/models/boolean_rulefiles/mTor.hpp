RULE_SET(orig, 22, "ERK","P27361","RSK","Q15418","MiR1976","MLL","Q03164","HOAX9","MiR196b","TSC12","P31749","PDK1","Q15118","AKT","RHEB","Q15382","GBL","Q9BVC4","mTOR","P42345","RICTOR","Q6R327")

template<> bool orig::rule<0>(bitset<22> x) {
    return x[0];
}
template<> vector<ind> orig::depends_on<0>() {
    return { 0 };
}

template<> bool orig::rule<1>(bitset<22> x) {
    return x[0];
}
template<> vector<ind> orig::depends_on<1>() {
    return { 0 };
}

template<> bool orig::rule<2>(bitset<22> x) {
    return x[2];
}
template<> vector<ind> orig::depends_on<2>() {
    return { 2 };
}

template<> bool orig::rule<3>(bitset<22> x) {
    return x[2];
}
template<> vector<ind> orig::depends_on<3>() {
    return { 2 };
}

template<> bool orig::rule<4>(bitset<22> x) {
    return x[2];
}
template<> vector<ind> orig::depends_on<4>() {
    return { 2 };
}

template<> bool orig::rule<5>(bitset<22> x) {
    return x[5];
}
template<> vector<ind> orig::depends_on<5>() {
    return { 5 };
}

template<> bool orig::rule<6>(bitset<22> x) {
    return x[5] && !x[4];
}
template<> vector<ind> orig::depends_on<6>() {
    return { 4,5 };
}

template<> bool orig::rule<7>(bitset<22> x) {
    return x[6];
}
template<> vector<ind> orig::depends_on<7>() {
    return { 6 };
}

template<> bool orig::rule<8>(bitset<22> x) {
    return x[7];
}
template<> vector<ind> orig::depends_on<8>() {
    return { 7 };
}

template<> bool orig::rule<9>(bitset<22> x) {
    return !x[10] && !(x[1] && x[3]);
}
template<> vector<ind> orig::depends_on<9>() {
    return { 1,3,10 };
}

template<> bool orig::rule<10>(bitset<22> x) {
    return x[13];
}
template<> vector<ind> orig::depends_on<10>() {
    return { 13 };
}

template<> bool orig::rule<11>(bitset<22> x) {
    return x[11];
}
template<> vector<ind> orig::depends_on<11>() {
    return { 11 };
}

template<> bool orig::rule<12>(bitset<22> x) {
    return x[11];
}
template<> vector<ind> orig::depends_on<12>() {
    return { 11 };
}

template<> bool orig::rule<13>(bitset<22> x) {
    return x[19] && x[12] && x[21] && x[17];
}
template<> vector<ind> orig::depends_on<13>() {
    return { 12,17,19,21 };
}

template<> bool orig::rule<14>(bitset<22> x) {
    return !x[9];
}
template<> vector<ind> orig::depends_on<14>() {
    return { 9 };
}

template<> bool orig::rule<15>(bitset<22> x) {
    return x[14];
}
template<> vector<ind> orig::depends_on<15>() {
    return { 14 };
}

template<> bool orig::rule<16>(bitset<22> x) {
    return x[15];
}
template<> vector<ind> orig::depends_on<16>() {
    return { 15 };
}

template<> bool orig::rule<17>(bitset<22> x) {
    return x[16];
}
template<> vector<ind> orig::depends_on<17>() {
    return { 16 };
}

template<> bool orig::rule<18>(bitset<22> x) {
    return x[15];
}
template<> vector<ind> orig::depends_on<18>() {
    return { 15 };
}

template<> bool orig::rule<19>(bitset<22> x) {
    return x[18];
}
template<> vector<ind> orig::depends_on<19>() {
    return { 18 };
}

template<> bool orig::rule<20>(bitset<22> x) {
    return x[15];
}
template<> vector<ind> orig::depends_on<20>() {
    return { 15 };
}

template<> bool orig::rule<21>(bitset<22> x) {
    return x[20] && !x[8];
}
template<> vector<ind> orig::depends_on<21>() {
    return { 8,20 };
}
