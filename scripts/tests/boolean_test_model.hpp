RULE_SET(TEST, 5, "S0","S1","S2","S3","S4")

template<> bool TEST::rule<0>(bitset<5> x) {
    return !x[1];
}
template<> vector<ind> TEST::depends_on<0>() {
    return { 1 };
}

template<> bool TEST::rule<1>(bitset<5> x) {
    return x[0] && (x[2] || x[3]);
}
template<> vector<ind> TEST::depends_on<1>() {
    return { 0, 2, 3 };
}

template<> bool TEST::rule<2>(bitset<5> x) {
    return x[2] || (x[2] && x[4]);
}
template<> vector<ind> TEST::depends_on<2>() {
    return { 2, 4 };
}

template<> bool TEST::rule<3>(bitset<5> x) {
    return (!x[2]) && x[3];
}
template<> vector<ind> TEST::depends_on<3>() {
    return { 2, 3 };
}

template<> bool TEST::rule<4>(bitset<5> x) {
    return x[1] && (!x[4]);
}
template<> vector<ind> TEST::depends_on<4>() {
    return { 1, 4 };
}