#include "tree_class.hpp"

// TODO: Add cast functions `to_internal_node`, `to_external_node`

void cme_internal_node::Initialize(int ncid)
{
    internal_node::Initialize(ncid);

    // initialize A and B for the root
    if (parent == nullptr) {
        for (Index mu = 0; mu < grid.n_reactions; ++mu) {
            coefficients.A[mu](0, 0) = 1.0;
            coefficients.B[mu](0, 0) = 1.0;
            coefficients.A_bar[mu](0, 0) = 1.0;
            coefficients.B_bar[mu](0, 0) = 1.0;
        }
    }
    return;
}

void cme_external_node::Initialize(int ncid)
{
    external_node::Initialize(ncid);

    // read propensity
    propensity = ReadHelpers::ReadPropensity(ncid, grid.n_reactions);
    return;
}

void cme_lr_tree::Read(const std::string fn)
{
    int ncid, grp_ncid;

    NETCDF_CHECK(nc_open(fn.c_str(), NC_NOWRITE, &ncid));
    partition_str = ReadHelpers::ReadPartitionStr(ncid);
    species_names = ReadHelpers::ReadSpeciesNames(ncid);
    grid_parms grid = ReadHelpers::ReadGridParms(ncid);
    std::array<Index, 2> r_out = ReadHelpers::ReadRankOut(ncid);
    root = new cme_internal_node("root", nullptr, grid, 1, r_out, 1);
    root->Initialize(ncid);

    NETCDF_CHECK(nc_inq_ncid(ncid, "0", &grp_ncid));
    root->child[0] = ReadHelpers::ReadNode(grp_ncid, "0", root, root->RankOut()[0]);

    NETCDF_CHECK(nc_inq_ncid(ncid, "1", &grp_ncid));
    root->child[1] = ReadHelpers::ReadNode(grp_ncid, "1", root, root->RankOut()[1]);

    NETCDF_CHECK(nc_close(ncid));
    return;
}

void cme_lr_tree::Write(const std::string fn, const double t, const double tau,
                        const double dm) const
{
    int ncid;

    NETCDF_CHECK(nc_create(fn.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid));
    // `WriteNode` sets netCDF dim `d`, therefore it must be called before
    // `WriteSpeciesNames`
    WriteHelpers::WriteNode(ncid, root);
    WriteHelpers::WritePartitionStr(ncid, partition_str);
    WriteHelpers::WriteSpeciesNames(ncid, species_names);

    int varid_tt, varid_tau, varid_dm;
    NETCDF_CHECK(nc_def_var(ncid, "t", NC_DOUBLE, 0, 0, &varid_tt));
    NETCDF_CHECK(nc_put_var_double(ncid, varid_tt, &t));
    NETCDF_CHECK(nc_def_var(ncid, "tau", NC_DOUBLE, 0, 0, &varid_tau));
    NETCDF_CHECK(nc_put_var_double(ncid, varid_tau, &tau));
    NETCDF_CHECK(nc_def_var(ncid, "dm", NC_DOUBLE, 0, 0, &varid_dm));
    NETCDF_CHECK(nc_put_var_double(ncid, varid_dm, &dm));

    NETCDF_CHECK(nc_close(ncid));
}

void WriteHelpers::WritePartitionStr(int ncid, const std::string partition_str)
{
    size_t len = partition_str.size() + 1;
    NETCDF_CHECK(
        nc_put_att_text(ncid, NC_GLOBAL, "partition_str", len, partition_str.data()));
}

void WriteHelpers::WriteSpeciesNames(int ncid,
                                     const std::vector<std::string> species_names)
{
    int id_d;
    NETCDF_CHECK(nc_inq_dimid(ncid, "d", &id_d));

    std::vector<const char*> species_names_char(species_names.size());
    std::transform(species_names.begin(), species_names.end(),
                   species_names_char.begin(), std::mem_fun_ref(&std::string::c_str));

    int varid_species_names;
    NETCDF_CHECK(
        nc_def_var(ncid, "species_names", NC_STRING, 1, &id_d, &varid_species_names));
    NETCDF_CHECK(
        nc_put_var_string(ncid, varid_species_names, species_names_char.data()));
}

void WriteHelpers::WriteGridParms(int ncid, const grid_parms grid)
{
    // Dimensions
    int id_d, id_dx;
    NETCDF_CHECK(nc_def_dim(ncid, "d", grid.d, &id_d));
    NETCDF_CHECK(nc_def_dim(ncid, "dx", grid.dx, &id_dx));

    std::vector<long long> n(std::begin(grid.n), std::end(grid.n));
    std::vector<long long> binsize(std::begin(grid.binsize), std::end(grid.binsize));

    int varid_n, varid_binsize, varid_liml, varid_species;
    NETCDF_CHECK(nc_def_var(ncid, "n", NC_INT64, 1, &id_d, &varid_n));
    NETCDF_CHECK(nc_put_var_longlong(ncid, varid_n, n.data()));
    NETCDF_CHECK(nc_def_var(ncid, "binsize", NC_INT64, 1, &id_d, &varid_binsize));
    NETCDF_CHECK(nc_put_var_longlong(ncid, varid_binsize, binsize.data()));
    NETCDF_CHECK(nc_def_var(ncid, "liml", NC_DOUBLE, 1, &id_d, &varid_liml));
    NETCDF_CHECK(nc_put_var_double(ncid, varid_liml, grid.liml.data()));
    NETCDF_CHECK(nc_def_var(ncid, "species", NC_INT, 1, &id_d, &varid_species));
    NETCDF_CHECK(nc_put_var_int(ncid, varid_species, grid.species.data()));
}

void WriteHelpers::WriteNode(int ncid, cme_node const* const node)
{
    int id_r_in;
    NETCDF_CHECK(nc_def_dim(ncid, "r_in", node->RankIn(), &id_r_in));
    WriteHelpers::WriteGridParms(ncid, node->grid);

    if (node->IsExternal()) {
        cme_external_node* this_node = (cme_external_node*)node;

        int id_dx;
        NETCDF_CHECK(nc_inq_dimid(ncid, "dx", &id_dx));
        this_node->Write(ncid, id_r_in, id_dx);
    }
    else {
        cme_internal_node* this_node = (cme_internal_node*)node;

        // Write r_out
        std::array<int, 2> id_r_out;
        NETCDF_CHECK(
            nc_def_dim(ncid, "r_out_0", this_node->RankOut()[0], &id_r_out[0]));
        NETCDF_CHECK(
            nc_def_dim(ncid, "r_out_1", this_node->RankOut()[1], &id_r_out[1]));

        this_node->Write(ncid, id_r_in, id_r_out);

        int grp_ncid0, grp_ncid1;
        NETCDF_CHECK(nc_def_grp(ncid, this_node->child[0]->id.c_str(), &grp_ncid0));
        NETCDF_CHECK(nc_def_grp(ncid, this_node->child[1]->id.c_str(), &grp_ncid1));
        WriteNode(grp_ncid0, this_node->child[0]);
        WriteNode(grp_ncid1, this_node->child[1]);
    }
}

using Species = const std::vector<int>;

std::ostream& operator<<(std::ostream& os, Species& species)
{
    os << '[';
    std::copy(std::begin(species), std::end(species),
              std::ostream_iterator<int>(os, ", "));
    os << "\b\b]"; // use two backspace characters '\b' to overwrite final ", "
    return os;
}

// TODO: write virtual print functions and use them in this function
void cme_lr_tree::PrintHelper(std::ostream& os, cme_node const* const node) const
{
    if (node->IsExternal()) {
        os << "external_node, id: " << node->id << ", species: " << node->grid.species
           << ", X.shape(): (" << ((cme_external_node*)node)->X.shape()[0] << ","
           << ((cme_external_node*)node)->X.shape()[1] << ")\n";
    }
    else {
        cme_lr_tree::PrintHelper(os, node->child[0]);
        cme_lr_tree::PrintHelper(os, node->child[1]);
        os << "internal_node, id: " << node->id << ", species: " << node->grid.species
           << ", rank_out: (" << ((cme_internal_node*)node)->RankOut()[0] << ","
           << ((cme_internal_node*)node)->RankOut()[1] << ")\n";
    }
}

void cme_lr_tree::OrthogonalizeHelper(cme_internal_node* const node,
                                      const Ensign::Matrix::blas_ops& blas) const
{
    std::function<double(double*, double*)> ip0;
    std::function<double(double*, double*)> ip1;
    Ensign::multi_array<double, 2> R0({node->RankOut()[0], node->RankOut()[0]});
    Ensign::multi_array<double, 2> R1({node->RankOut()[1], node->RankOut()[1]});
    if (node->child[0]->IsExternal() and node->child[1]->IsExternal()) {
        cme_external_node* node_left = (cme_external_node*)node->child[0];
        cme_external_node* node_right = (cme_external_node*)node->child[1];

        // ip0 = Ensign::inner_product_from_const_weight(node_left->grid.h_mult,
        // node_left->grid.dx); ip1 =
        // Ensign::inner_product_from_const_weight(node_right->grid.h_mult,
        // node_right->grid.dx);

        R0 = node_left->orthogonalize(node_left->grid.h_mult, blas);
        R1 = node_right->orthogonalize(node_right->grid.h_mult, blas);
    }
    else if (node->child[0]->IsExternal() and node->child[1]->IsInternal()) {
        cme_external_node* node_left = (cme_external_node*)node->child[0];
        cme_internal_node* node_right = (cme_internal_node*)node->child[1];

        OrthogonalizeHelper(node_right, blas);

        // ip0 = Ensign::inner_product_from_const_weight(node_left->grid.h_mult,
        // node_left->grid.dx); ip1 = Ensign::inner_product_from_const_weight(1.0,
        // Ensign::prod(node_right->RankOut()));

        R0 = node_left->orthogonalize(node_left->grid.h_mult, blas);
        R1 = node_right->orthogonalize(1.0, blas);
    }
    else if (node->child[0]->IsInternal() and node->child[1]->IsExternal()) {
        cme_internal_node* node_left = (cme_internal_node*)node->child[0];
        cme_external_node* node_right = (cme_external_node*)node->child[1];

        OrthogonalizeHelper(node_left, blas);

        // ip0 = Ensign::inner_product_from_const_weight(1.0,
        // Ensign::prod(node_left->RankOut())); ip1 =
        // Ensign::inner_product_from_const_weight(node_right->grid.h_mult,
        // node_right->grid.dx);

        R0 = node_left->orthogonalize(1.0, blas);
        R1 = node_right->orthogonalize(node_right->grid.h_mult, blas);
    }
    else {
        cme_internal_node* node_left = (cme_internal_node*)node->child[0];
        cme_internal_node* node_right = (cme_internal_node*)node->child[1];

        OrthogonalizeHelper(node_left, blas);
        OrthogonalizeHelper(node_right, blas);

        // ip0 = Ensign::inner_product_from_const_weight(1.0,
        // Ensign::prod(node_left->RankOut())); ip1 =
        // Ensign::inner_product_from_const_weight(1.0,
        // Ensign::prod(node_right->RankOut()));

        R0 = node_left->orthogonalize(1.0, blas);
        R1 = node_right->orthogonalize(1.0, blas);
    }

    for (int j = node->child[0]->n_basisfunctions; j < node->RankOut()[0]; ++j) {
        for (int i = 0; i < node->RankOut()[0]; ++i) {
            R0(i, j) = 0.0;
        }
    }

    for (int j = node->child[1]->n_basisfunctions; j < node->RankOut()[1]; ++j) {
        for (int i = 0; i < node->RankOut()[1]; ++i) {
            R1(i, j) = 0.0;
        }
    }

    Ensign::multi_array<double, 2> tmp1(node->RankOut());
    Ensign::multi_array<double, 2> tmp2(node->RankOut());

    for (Index k = 0; k < node->RankIn(); ++k) {
        for (Index j = 0; j < node->RankOut()[1]; ++j) {
            for (Index i = 0; i < node->RankOut()[0]; ++i) {
                tmp1(i, j) = node->Q(i, j, k);
            }
        }
        blas.matmul(R0, tmp1, tmp2);
        blas.matmul_transb(tmp2, R1, tmp1);
        for (Index j = 0; j < node->RankOut()[1]; ++j) {
            for (Index i = 0; i < node->RankOut()[0]; ++i) {
                node->Q(i, j, k) = tmp1(i, j);
            }
        }
    }
}

void cme_lr_tree::Orthogonalize(const Ensign::Matrix::blas_ops& blas) const
{
    OrthogonalizeHelper(root, blas);
};

void cme_lr_tree::InitializeAB_barHelper(cme_node* const node,
                                         const Ensign::Matrix::blas_ops& blas) const
{
    if (node->IsInternal()) {
        InitializeAB_barHelper(node->child[0], blas);
        InitializeAB_barHelper(node->child[1], blas);
    }
    node->CalculateAB_bar(blas);
}

void cme_lr_tree::InitializeAB_bar(const Ensign::Matrix::blas_ops& blas) const
{
    InitializeAB_barHelper(root, blas);
}

std::vector<double> cme_lr_tree::NormalizeHelper(cme_node const* const node) const
{
    std::vector<double> x_sum(node->RankIn(), 0.0);
    if (node->IsExternal()) {
        cme_external_node* this_node = (cme_external_node*)node;
        for (Index i = 0; i < this_node->RankIn(); ++i) {
            for (Index x = 0; x < this_node->ProblemSize(); ++x) {
                x_sum[i] += this_node->X(x, i) * this_node->grid.h_mult;
            }
        }
    }
    else {
        cme_internal_node* this_node = (cme_internal_node*)node;
        std::vector<double> x0_sum(this_node->RankOut()[0]),
            x1_sum(this_node->RankOut()[1]);
        x0_sum = NormalizeHelper(this_node->child[0]);
        x1_sum = NormalizeHelper(this_node->child[1]);

        for (Index i = 0; i < this_node->RankIn(); ++i) {
            for (Index i0 = 0; i0 < this_node->RankOut()[0]; ++i0) {
                for (Index i1 = 0; i1 < this_node->RankOut()[1]; ++i1) {
                    x_sum[i] += this_node->Q(i0, i1, i) * x0_sum[i0] * x1_sum[i1];
                }
            }
        }
    }
    return x_sum;
}

double cme_lr_tree::Normalize() const
{
    std::vector<double> norm(1);
    norm = NormalizeHelper(root);
    root->Q /= norm[0];
    return norm[0];
}

std::string ReadHelpers::ReadPartitionStr(int ncid)
{
    size_t partition_str_len;
    NETCDF_CHECK(nc_inq_attlen(ncid, NC_GLOBAL, "partition_str", &partition_str_len));

    std::string partition_str(partition_str_len, '\0');
    NETCDF_CHECK(
        nc_get_att_text(ncid, NC_GLOBAL, "partition_str", partition_str.data()));

    return partition_str;
}

std::vector<std::string> ReadHelpers::ReadSpeciesNames(int ncid)
{
    // read dimension
    int id_d;
    NETCDF_CHECK(nc_inq_dimid(ncid, "d", &id_d));

    size_t d_t;
    char tmp[NC_MAX_NAME + 1];
    NETCDF_CHECK(nc_inq_dim(ncid, id_d, tmp, &d_t));

    size_t start[] = {0};
    size_t count[] = {d_t};

    // read variable
    int id_species_names;

    NETCDF_CHECK(nc_inq_varid(ncid, "species_names", &id_species_names));

    std::vector<char*> species_names_char(d_t);

    NETCDF_CHECK(nc_get_vara_string(ncid, id_species_names, start, count,
                                    species_names_char.data()));
    std::vector<std::string> species_names(std::begin(species_names_char),
                                           std::end(species_names_char));

    return species_names;
}

grid_parms ReadHelpers::ReadGridParms(int ncid)
{
    // read dimensions
    int id_d, id_n_reactions;
    NETCDF_CHECK(nc_inq_dimid(ncid, "d", &id_d));
    NETCDF_CHECK(nc_inq_dimid(ncid, "n_reactions", &id_n_reactions));

    size_t d_t, n_reactions_t;
    char tmp[NC_MAX_NAME + 1];
    NETCDF_CHECK(nc_inq_dim(ncid, id_d, tmp, &d_t));
    NETCDF_CHECK(nc_inq_dim(ncid, id_n_reactions, tmp, &n_reactions_t));

    Index d = (Index)d_t;
    Index n_reactions = (Index)n_reactions_t;

    // read variables
    int id_n, id_binsize, id_liml, id_dep, id_nu, id_species;

    NETCDF_CHECK(nc_inq_varid(ncid, "n", &id_n));
    NETCDF_CHECK(nc_inq_varid(ncid, "binsize", &id_binsize));
    NETCDF_CHECK(nc_inq_varid(ncid, "liml", &id_liml));
    NETCDF_CHECK(nc_inq_varid(ncid, "dep", &id_dep));
    NETCDF_CHECK(nc_inq_varid(ncid, "nu", &id_nu));
    NETCDF_CHECK(nc_inq_varid(ncid, "species", &id_species));

    grid_parms grid(d, n_reactions);
    Ensign::multi_array<signed char, 2> dep_int({n_reactions, d});

    NETCDF_CHECK(nc_get_var_long(ncid, id_n, grid.n.data()));
    NETCDF_CHECK(nc_get_var_long(ncid, id_binsize, grid.binsize.data()));
    NETCDF_CHECK(nc_get_var_double(ncid, id_liml, grid.liml.data()));
    NETCDF_CHECK(nc_get_var_schar(ncid, id_dep, dep_int.data()));
    NETCDF_CHECK(nc_get_var_long(ncid, id_nu, grid.nu.data()));
    NETCDF_CHECK(nc_get_var_int(ncid, id_species, grid.species.data()));

    std::copy(std::begin(dep_int), std::end(dep_int), std::begin(grid.dep));
    grid.Initialize();

    return grid;
}

// NOTE: Input .netCDF file is required to store a scalar `r_out`
// (which guarantees that the rank is equal for both childs)
std::array<Index, 2> ReadHelpers::ReadRankOut(int ncid)
{
    // read rank
    int id_r0, id_r1;
    NETCDF_CHECK(nc_inq_dimid(ncid, "r_out0", &id_r0));
    NETCDF_CHECK(nc_inq_dimid(ncid, "r_out1", &id_r1));

    size_t r0_t, r1_t;
    char tmp[NC_MAX_NAME + 1];
    NETCDF_CHECK(nc_inq_dim(ncid, id_r0, tmp, &r0_t));
    NETCDF_CHECK(nc_inq_dim(ncid, id_r1, tmp, &r1_t));

    std::array<Index, 2> r_out = {(Index)r0_t, (Index)r1_t};
    return r_out;
}

Index ReadHelpers::ReadNBasisfunctions(int ncid)
{
    // read number of basisfunctions
    int id_n_basisfunctions;
    NETCDF_CHECK(nc_inq_dimid(ncid, "n_basisfunctions", &id_n_basisfunctions));

    size_t n_basisfunctions_t;
    char tmp[NC_MAX_NAME + 1];
    NETCDF_CHECK(nc_inq_dim(ncid, id_n_basisfunctions, tmp, &n_basisfunctions_t));

    return (Index)n_basisfunctions_t;
}

std::vector<std::vector<double>> ReadHelpers::ReadPropensity(int ncid,
                                                             const Index n_reactions)
{
    std::vector<std::vector<double>> result(n_reactions);

    for (Index mu = 0; mu < n_reactions; ++mu) {
        // read dimension dx_{mu}
        int id_dx_dep;
        NETCDF_CHECK(
            nc_inq_dimid(ncid, ("dx_" + std::to_string(mu)).c_str(), &id_dx_dep));

        size_t dx_dep_t;
        char tmp[NC_MAX_NAME + 1];
        NETCDF_CHECK(nc_inq_dim(ncid, id_dx_dep, tmp, &dx_dep_t));

        result[mu].resize(dx_dep_t);

        // read propensity
        int id_propensity;

        NETCDF_CHECK(nc_inq_varid(ncid, ("propensity_" + std::to_string(mu)).c_str(),
                                  &id_propensity));
        NETCDF_CHECK(nc_get_var_double(ncid, id_propensity, result[mu].data()));
    }

    return result;
}

cme_node* ReadHelpers::ReadNode(int ncid, const std::string id,
                                cme_internal_node* const parent_node, const Index r_in)
{
    int retval0, retval1;
    int grp_ncid0, grp_ncid1;
    grid_parms grid = ReadGridParms(ncid);
    Index n_basisfunctions = ReadNBasisfunctions(ncid);
    cme_node* child_node;

    retval0 = nc_inq_ncid(ncid, (id + "0").c_str(), &grp_ncid0);
    retval1 = nc_inq_ncid(ncid, (id + "1").c_str(), &grp_ncid1);

    if (retval0 or retval1) // NOTE: if retval0 is 1, then retval1 has also to be 1, due
                            // to the binary tree structure of the netCDF file
    {
        child_node =
            new cme_external_node(id, parent_node, grid, r_in, n_basisfunctions);
        child_node->Initialize(ncid);
    }
    else {
        std::array<Index, 2> r_out = ReadRankOut(ncid);
        child_node =
            new cme_internal_node(id, parent_node, grid, r_in, r_out, n_basisfunctions);
        child_node->Initialize(ncid);
        child_node->child[0] =
            ReadNode(grp_ncid0, id + "0", (cme_internal_node*)child_node,
                     ((cme_internal_node*)child_node)->RankOut()[0]);
        child_node->child[1] =
            ReadNode(grp_ncid1, id + "1", (cme_internal_node*)child_node,
                     ((cme_internal_node*)child_node)->RankOut()[1]);
    }
    return child_node;
}

void cme_node::CalculateAB_bar(const Ensign::Matrix::blas_ops& blas)
{
    if (IsExternal()) {
#ifdef __OPENMP__
#pragma omp parallel
#endif
        {
            cme_external_node* this_node = (cme_external_node*)this;
            Ensign::multi_array<double, 2> X_shift(
                {this_node->grid.dx, this_node->RankIn()});
            Ensign::multi_array<double, 1> weight({this_node->grid.dx});

#ifdef __OPENMP__
#pragma omp for
#endif
            for (Index mu = 0; mu < this_node->grid.n_reactions; ++mu) {
                Matrix::ShiftRows<-1>(X_shift, this_node->X, this_node->grid, mu);
                std::vector<Index> vec_index(this_node->grid.d, 0);

                for (Index alpha = 0; alpha < this_node->grid.dx; ++alpha) {
                    Index alpha_dep = VecIndexToDepCombIndex(
                        std::begin(vec_index), std::begin(this_node->grid.n_dep[mu]),
                        std::begin(this_node->grid.idx_dep[mu]),
                        std::end(this_node->grid.idx_dep[mu]));

                    weight(alpha) =
                        this_node->propensity[mu][alpha_dep] * this_node->grid.h_mult;
                    Ensign::IndexFunction::incr_vec_index(std::begin(this_node->grid.n),
                                                          std::begin(vec_index),
                                                          std::end(vec_index));
                }
                coeff(X_shift, this_node->X, weight, coefficients.A_bar[mu], blas);
                coeff(this_node->X, this_node->X, weight, coefficients.B_bar[mu], blas);
            }
        }
    }
    else {
        cme_internal_node* this_node = (cme_internal_node*)this;
        std::array<Index, 2> rank_out = this_node->RankOut();

#ifdef __OPENMP__
#pragma omp parallel
#endif
        {
            Ensign::multi_array<double, 2> QA0_tilde(
                {Ensign::prod(rank_out), this_node->RankIn()});
            Ensign::multi_array<double, 2> QA1_tilde(
                {Ensign::prod(rank_out), this_node->RankIn()});

            Ensign::multi_array<double, 2> Q0_mat(
                {this_node->RankIn() * rank_out[1], rank_out[0]});
            Ensign::multi_array<double, 2> Q1_mat(
                {this_node->RankIn() * rank_out[0], rank_out[1]});
            Ensign::multi_array<double, 2> QA0_mat(Q0_mat.shape());
            Ensign::multi_array<double, 2> QA1_mat(Q1_mat.shape());

            Ensign::multi_array<double, 3> QA0(this_node->Q.shape());
            Ensign::multi_array<double, 3> QA1(this_node->Q.shape());

            Ensign::Tensor::matricize<0>(this_node->Q, Q0_mat);
            Ensign::Tensor::matricize<1>(this_node->Q, Q1_mat);

#ifdef __OPENMP__
#pragma omp for
#endif
            for (Index mu = 0; mu < this_node->grid.n_reactions; ++mu) {
                blas.matmul(Q0_mat, this_node->child[0]->coefficients.A_bar[mu],
                            QA0_mat);
                Ensign::Tensor::tensorize<0>(QA0_mat, QA0);
                Ensign::Tensor::matricize<2>(QA0, QA0_tilde);
                blas.matmul_transb(Q1_mat, this_node->child[1]->coefficients.A_bar[mu],
                                   QA1_mat);
                Ensign::Tensor::tensorize<1>(QA1_mat, QA1);
                Ensign::Tensor::matricize<2>(QA1, QA1_tilde);
                blas.matmul_transa(QA0_tilde, QA1_tilde, coefficients.A_bar[mu]);

                blas.matmul(Q0_mat, this_node->child[0]->coefficients.B_bar[mu],
                            QA0_mat);
                Ensign::Tensor::tensorize<0>(QA0_mat, QA0);
                Ensign::Tensor::matricize<2>(QA0, QA0_tilde);
                blas.matmul_transb(Q1_mat, this_node->child[1]->coefficients.B_bar[mu],
                                   QA1_mat);
                Ensign::Tensor::tensorize<1>(QA1_mat, QA1);
                Ensign::Tensor::matricize<2>(QA1, QA1_tilde);
                blas.matmul_transa(QA0_tilde, QA1_tilde, coefficients.B_bar[mu]);
            }
        }
    }
}

Ensign::multi_array<double, 2> CalculateKDot(const Ensign::multi_array<double, 2>& K,
                                             const cme_external_node* const node,
                                             const Ensign::Matrix::blas_ops& blas)
{
    Ensign::gt::start("CalculateKDot");
    Ensign::multi_array<double, 2> K_dot(K.shape());
    Ensign::Matrix::set_zero(K_dot);

#ifdef __OPENMP__
#pragma omp parallel reduction(+ : K_dot)
#endif
    {
        Ensign::multi_array<double, 2> K_dot_thread(K.shape());
        Ensign::multi_array<double, 2> KA(K.shape());
        Ensign::multi_array<double, 2> KB(K.shape());
        Ensign::multi_array<double, 2> aKA(K.shape());
        Ensign::multi_array<double, 2> aKB(K.shape());
        Ensign::multi_array<double, 2> aKA_shift(K.shape());
        Ensign::multi_array<double, 1> weight({K.shape()[0]});
        Ensign::Matrix::set_zero(K_dot_thread);

#ifdef __OPENMP__
#pragma omp for
#endif
        for (Index mu = 0; mu < node->grid.n_reactions; ++mu) {
            blas.matmul_transb(K, node->coefficients.A[mu], KA);
            blas.matmul_transb(K, node->coefficients.B[mu], KB);

            std::vector<Index> vec_index(node->grid.d, 0);

            for (Index i = 0; i < node->grid.dx; ++i) {
                Index alpha = VecIndexToDepCombIndex(std::begin(vec_index),
                                                     std::begin(node->grid.n_dep[mu]),
                                                     std::begin(node->grid.idx_dep[mu]),
                                                     std::end(node->grid.idx_dep[mu]));
                Ensign::IndexFunction::incr_vec_index(std::begin(node->grid.n),
                                                      std::begin(vec_index),
                                                      std::end(vec_index));

                weight(i) = node->propensity[mu][alpha];
            }

            Ensign::Matrix::ptw_mult_row(KA, weight, aKA);
            Ensign::Matrix::ptw_mult_row(KB, weight, aKB);

            // Shift prod_KC
            Ensign::gt::start("ShiftKDot");
            Matrix::ShiftRows<1>(aKA_shift, aKA, node->grid, mu);
            Ensign::gt::stop("ShiftKDot");

            // Calculate K_dot = shift(C * K) - D * K
            aKA_shift -= aKB;
            K_dot_thread += aKA_shift;
        }
        K_dot += K_dot_thread;
    }

    Ensign::gt::stop("CalculateKDot");
    return K_dot;
}

Ensign::multi_array<double, 2> CalculateSDot(const Ensign::multi_array<double, 2>& S,
                                             const cme_node* const node,
                                             const Ensign::Matrix::blas_ops& blas)
{
    Ensign::multi_array<double, 2> S_dot(S.shape());
    Ensign::Matrix::set_zero(S_dot);

#ifdef __OPENMP__
#pragma omp parallel reduction(+ : S_dot)
#endif
    {
        Ensign::multi_array<double, 2> S_dot_thread(S_dot);
        Ensign::multi_array<double, 2> AS(S_dot);
        Ensign::multi_array<double, 2> ASa(S_dot);
        Ensign::multi_array<double, 2> BS(S_dot);
        Ensign::multi_array<double, 2> BSb(S_dot);

#ifdef __OPENMP__
#pragma omp for
#endif
        for (Index mu = 0; mu < node->grid.n_reactions; ++mu) {
            blas.matmul(node->coefficients.A_bar[mu], S, AS);
            blas.matmul_transb(AS, node->coefficients.A[mu], ASa);
            blas.matmul(node->coefficients.B_bar[mu], S, BS);
            blas.matmul_transb(BS, node->coefficients.B[mu], BSb);
            S_dot_thread += ASa;
            S_dot_thread -= BSb;
        }
        S_dot += S_dot_thread;
    }

    return S_dot;
}

Ensign::multi_array<double, 2> CalculateQDot(const Ensign::multi_array<double, 2>& Qmat,
                                             const cme_internal_node* const node,
                                             const Ensign::Matrix::blas_ops& blas)
{
    Ensign::multi_array<double, 2> Q_dot(Qmat.shape());
    Ensign::Matrix::set_zero(Q_dot);

#ifdef __OPENMP__
#pragma omp parallel reduction(+ : Q_dot)
#endif
    {
        Ensign::multi_array<double, 3> Q_ten(
            {node->RankOut()[0], node->RankOut()[1], node->RankIn()});
        Ensign::multi_array<double, 2> Q_dot_thread(Q_dot);

        Ensign::multi_array<double, 2> Qa(Q_dot);
        Ensign::multi_array<double, 2> Qa_mat(
            {node->RankOut()[1] * node->RankIn(), node->RankOut()[0]});
        Ensign::multi_array<double, 2> QaA0(Qa_mat.shape());
        Ensign::multi_array<double, 2> QaA0_mat(
            {node->RankOut()[0] * node->RankIn(), node->RankOut()[1]});
        Ensign::multi_array<double, 2> QaA0A1(QaA0_mat.shape());

        Ensign::multi_array<double, 2> Qb(Q_dot);
        Ensign::multi_array<double, 2> Qb_mat(
            {node->RankOut()[1] * node->RankIn(), node->RankOut()[0]});
        Ensign::multi_array<double, 2> QbB0(Qb_mat.shape());
        Ensign::multi_array<double, 2> QbB0_mat(
            {node->RankOut()[0] * node->RankIn(), node->RankOut()[1]});
        Ensign::multi_array<double, 2> QbB0B1(QbB0_mat.shape());

#ifdef __OPENMP__
#pragma omp for
#endif
        for (Index mu = 0; mu < node->grid.n_reactions; ++mu) {
            blas.matmul_transb(Qmat, node->coefficients.A[mu], Qa);
            Ensign::gt::start("Mat/Ten");
            Ensign::Tensor::tensorize<2>(Qa, Q_ten);
            Ensign::Tensor::matricize<0>(Q_ten, Qa_mat);
            Ensign::gt::stop("Mat/Ten");
            blas.matmul_transb(Qa_mat, node->child[0]->coefficients.A_bar[mu], QaA0);
            Ensign::gt::start("Mat/Ten");
            Ensign::Tensor::tensorize<0>(QaA0, Q_ten);
            Ensign::Tensor::matricize<1>(Q_ten, QaA0_mat);
            Ensign::gt::stop("Mat/Ten");
            blas.matmul_transb(QaA0_mat, node->child[1]->coefficients.A_bar[mu],
                               QaA0A1);
            Ensign::gt::start("Mat/Ten");
            Ensign::Tensor::tensorize<1>(QaA0A1, Q_ten);
            Ensign::Tensor::matricize<2>(Q_ten, Qa);
            Ensign::gt::stop("Mat/Ten");

            blas.matmul_transb(Qmat, node->coefficients.B[mu], Qb);
            Ensign::gt::start("Mat/Ten");
            Ensign::Tensor::tensorize<2>(Qb, Q_ten);
            Ensign::Tensor::matricize<0>(Q_ten, Qb_mat);
            Ensign::gt::stop("Mat/Ten");
            blas.matmul_transb(Qb_mat, node->child[0]->coefficients.B_bar[mu], QbB0);
            Ensign::gt::start("Mat/Ten");
            Ensign::Tensor::tensorize<0>(QbB0, Q_ten);
            Ensign::Tensor::matricize<1>(Q_ten, QbB0_mat);
            Ensign::gt::stop("Mat/Ten");
            blas.matmul_transb(QbB0_mat, node->child[1]->coefficients.B_bar[mu],
                               QbB0B1);
            Ensign::gt::start("Mat/Ten");
            Ensign::Tensor::tensorize<1>(QbB0B1, Q_ten);
            Ensign::Tensor::matricize<2>(Q_ten, Qb);
            Ensign::gt::stop("Mat/Ten");

            Q_dot_thread += Qa;
            Q_dot_thread -= Qb;
        }
        Q_dot += Q_dot_thread;
    }

    return Q_dot;
}
