#include "tree_class.hpp"

// TODO: Add cast functions `to_internal_node`, `to_external_node`

template<>
void internal_node<double>::Initialize(int ncid)
{
    // read Q
    int id_Q;
    NETCDF_CHECK(nc_inq_varid(ncid, "Q", &id_Q));
    NETCDF_CHECK(nc_get_var_double(ncid, id_Q, Q.data()));

    return;
}

void cme_internal_node::Initialize(int ncid)
{
    internal_node::Initialize(ncid);

    // initialize A, B, E and F for the root
    if (parent == nullptr)
    {
        for (Index mu = 0; mu < grid.n_reactions; ++mu)
        {
            coefficients.A(mu, 0, 0) = 1.0;
            coefficients.B(mu, 0, 0) = 1.0;
        }
        coefficients.E(0, 0, 0, 0) = 1.0;
        coefficients.F(0, 0, 0, 0) = 1.0;
    }
    return;
}

template<>
void external_node<double>::Initialize(int ncid)
{
    // read X
    int id_X;
    NETCDF_CHECK(nc_inq_varid(ncid, "X", &id_X));
    NETCDF_CHECK(nc_get_var_double(ncid, id_X, X.data()));

    return;
}

void cme_external_node::Initialize(int ncid)
{
    external_node::Initialize(ncid);

    // read propensity
    external_coefficients.propensity = ReadHelpers::ReadPropensity(ncid, grid.n_reactions);

    // resize C and D coefficients
    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        external_coefficients.C[mu].resize({grid.dx_dep[mu], RankIn(), RankIn()});
        external_coefficients.D[mu].resize({grid.dx_dep[mu], RankIn(), RankIn()});
    }
    return;
}

template<>
void internal_node<double>::Write(int ncid, int id_r_in, std::array<int, 2> id_r_out) const
{
    int varid_Q;
    // NOTE: netCDF stores arrays in row-major order (but Ensign column-major order)
    int dimids_Q[3] = {id_r_in, id_r_out[1], id_r_out[0]};
    NETCDF_CHECK(nc_def_var(ncid, "Q", NC_DOUBLE, 3, dimids_Q, &varid_Q));
    NETCDF_CHECK(nc_put_var_double(ncid, varid_Q, Q.data()));
}

template<>
void external_node<double>::Write(int ncid, int id_r_in, int id_dx) const
{
    int varid_X;
    // NOTE: netCDF stores arrays in row-major order (but Ensign column-major order)
    int dimids_X[2] = {id_r_in, id_dx};
    NETCDF_CHECK(nc_def_var(ncid, "X", NC_DOUBLE, 2, dimids_X, &varid_X));
    NETCDF_CHECK(nc_put_var_double(ncid, varid_X, X.data()));
}

void cme_lr_tree::Read(const std::string fn)
{
    int ncid, grp_ncid;

    NETCDF_CHECK(nc_open(fn.c_str(), NC_NOWRITE, &ncid));
    partition_str = ReadHelpers::ReadPartitionStr(ncid);
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

void cme_lr_tree::Write(const std::string fn, const double t, const double tau) const
{
    int ncid;

    NETCDF_CHECK(nc_create(fn.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid));
    WriteHelpers::WritePartitionStr(ncid, partition_str);

    int varid_tt, varid_tau;
    NETCDF_CHECK(nc_def_var(ncid, "t", NC_DOUBLE, 0, 0, &varid_tt));
    NETCDF_CHECK(nc_put_var_double(ncid, varid_tt, &t));
    NETCDF_CHECK(nc_def_var(ncid, "tau", NC_DOUBLE, 0, 0, &varid_tau));
    NETCDF_CHECK(nc_put_var_double(ncid, varid_tau, &tau));

    WriteHelpers::WriteNode(ncid, root);

    NETCDF_CHECK(nc_close(ncid));
}

void WriteHelpers::WritePartitionStr(int ncid, const std::string partition_str)
{
    size_t len = partition_str.size() + 1;
    NETCDF_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "partition_str", len, partition_str.data()));
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
    NETCDF_CHECK(nc_put_var_int(ncid, varid_liml, grid.species.data()));
}

void WriteHelpers::WriteNode(int ncid, cme_node const * const node)
{
    int id_r_in;
    NETCDF_CHECK(nc_def_dim(ncid, "r_in", node->RankIn(), &id_r_in));
    WriteHelpers::WriteGridParms(ncid, node->grid);

    if (node->IsExternal())
    {
        cme_external_node* this_node = (cme_external_node*) node;

        int id_dx;
        NETCDF_CHECK(nc_inq_dimid(ncid, "dx", &id_dx));
        this_node->Write(ncid, id_r_in, id_dx);
    }
    else
    {
        cme_internal_node *this_node = (cme_internal_node *)node;

        // Write r_out
        std::array<int, 2> id_r_out;
        NETCDF_CHECK(nc_def_dim(ncid, "r_out_0", this_node->RankOut()[0], &id_r_out[0]));
        NETCDF_CHECK(nc_def_dim(ncid, "r_out_1", this_node->RankOut()[1], &id_r_out[1]));

        this_node->Write(ncid, id_r_in, id_r_out);

        int grp_ncid0, grp_ncid1;
        NETCDF_CHECK(nc_def_grp(ncid, this_node->child[0]->id.c_str(), &grp_ncid0));
        NETCDF_CHECK(nc_def_grp(ncid, this_node->child[1]->id.c_str(), &grp_ncid1));
        WriteNode(grp_ncid0, this_node->child[0]);
        WriteNode(grp_ncid1, this_node->child[1]);
    }
}

using Species = const std::vector<int>;

std::ostream &operator<<(std::ostream &os, Species &species)
{
    os << '[';
    std::copy(std::begin(species), std::end(species), std::ostream_iterator<int>(os, ", "));
    os << "\b\b]"; // use two backspace characters '\b' to overwrite final ", "
    return os;
}

// TODO: write virtual print functions and use them in this function
void cme_lr_tree::PrintHelper(std::ostream &os, cme_node const * const node) const
{
    if (node->IsExternal())
    {
        os << "external_node, id: " << node->id << ", species: " << node->grid.species << ", X.shape(): (" << ((cme_external_node*) node)->X.shape()[0] << "," << ((cme_external_node*) node)->X.shape()[1] << ")\n";
    }
    else
    {
        cme_lr_tree::PrintHelper(os, node->child[0]);
        cme_lr_tree::PrintHelper(os, node->child[1]);
        os << "internal_node, id: " << node->id << ", species: " << node->grid.species << ", rank_out: (" << ((cme_internal_node *)node)->RankOut()[0] << "," << ((cme_internal_node *)node)->RankOut()[1] << ")\n";
    }
}

void cme_lr_tree::OrthogonalizeHelper(cme_internal_node * const node, const blas_ops &blas) const
{
    std::function<double(double *, double *)> ip0;
    std::function<double(double *, double *)> ip1;
    multi_array<double, 2> R0({node->RankOut()[0], node->RankOut()[0]});
    multi_array<double, 2> R1({node->RankOut()[1], node->RankOut()[1]});
    if (node->child[0]->IsExternal() and node->child[1]->IsExternal())
    {
        cme_external_node* node_left = (cme_external_node*) node->child[0];
        cme_external_node* node_right = (cme_external_node*) node->child[1];

        ip0 = inner_product_from_const_weight(node_left->grid.h_mult, node_left->grid.dx);
        ip1 = inner_product_from_const_weight(node_right->grid.h_mult, node_right->grid.dx);
        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }
    else if (node->child[0]->IsExternal() and node->child[1]->IsInternal())
    {
        cme_external_node *node_left = (cme_external_node *)node->child[0];
        cme_internal_node *node_right = (cme_internal_node *)node->child[1];

        OrthogonalizeHelper(node_right, blas);

        ip0 = inner_product_from_const_weight(node_left->grid.h_mult, node_left->grid.dx);
        ip1 = inner_product_from_const_weight(1.0, prod(node_right->RankOut()));

        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }
    else if (node->child[0]->IsInternal() and node->child[1]->IsExternal())
    {
        cme_internal_node *node_left = (cme_internal_node *)node->child[0];
        cme_external_node *node_right = (cme_external_node *)node->child[1];

        OrthogonalizeHelper(node_left, blas);

        ip0 = inner_product_from_const_weight(1.0, prod(node_left->RankOut()));
        ip1 = inner_product_from_const_weight(node_right->grid.h_mult, node_right->grid.dx);

        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }
    else
    {
        cme_internal_node *node_left = (cme_internal_node *)node->child[0];
        cme_internal_node *node_right = (cme_internal_node *)node->child[1];

        OrthogonalizeHelper(node_left, blas);
        OrthogonalizeHelper(node_right, blas);

        ip0 = inner_product_from_const_weight(1.0, prod(node_left->RankOut()));
        ip1 = inner_product_from_const_weight(1.0, prod(node_right->RankOut()));
        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }

    for (int j = node->child[0]->n_basisfunctions; j < node->RankOut()[0]; ++j)
    {
        for (int i = 0; i < node->RankOut()[0]; ++i)
        {
            R0(i, j) = 0.0;
        }
    }

    for (int j = node->child[1]->n_basisfunctions; j < node->RankOut()[1]; ++j)
    {
        for (int i = 0; i < node->RankOut()[1]; ++i)
        {
            R1(i, j) = 0.0;
        }
    }

    multi_array<double, 2> tmp1(node->RankOut());
    multi_array<double, 2> tmp2(node->RankOut());

    for (Index k = 0; k < node->RankIn(); ++k)
    {
        for (Index j = 0; j < node->RankOut()[1]; ++j)
        {
            for (Index i = 0; i < node->RankOut()[0]; ++i)
            {
                tmp1(i, j) = node->Q(i, j, k);
            }
        }
        blas.matmul(R0, tmp1, tmp2);
        blas.matmul_transb(tmp2, R1, tmp1);
        for (Index j = 0; j < node->RankOut()[1]; ++j)
        {
            for (Index i = 0; i < node->RankOut()[0]; ++i)
            {
                node->Q(i, j, k) = tmp1(i, j);
            }
        }
    }
}

void cme_lr_tree::Orthogonalize(const blas_ops& blas) const
{
    OrthogonalizeHelper(root, blas);
};

std::vector<double> cme_lr_tree::NormalizeHelper(cme_node const * const node) const
{
    std::vector<double> x_sum(node->RankIn(), 0.0);
    if (node->IsExternal())
    {
        cme_external_node* this_node = (cme_external_node*) node;
        for (Index i = 0; i < this_node->RankIn(); ++i)
        {
            for (Index x = 0; x < this_node->ProblemSize(); ++x)
            {
                x_sum[i] += this_node->X(x, i) * this_node->grid.h_mult;
            }
        }
    }
    else
    {
        cme_internal_node* this_node = (cme_internal_node*) node;
        std::vector<double> x0_sum(this_node->RankOut()[0]), x1_sum(this_node->RankOut()[1]);
        x0_sum = NormalizeHelper(this_node->child[0]);
        x1_sum = NormalizeHelper(this_node->child[1]);

        for (Index i = 0; i < this_node->RankIn(); ++i)
        {
            for (Index i0 = 0; i0 < this_node->RankOut()[0]; ++i0)
            {
                for (Index i1 = 0; i1 < this_node->RankOut()[1]; ++i1)
                {
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
    NETCDF_CHECK(nc_get_att_text(ncid, NC_GLOBAL, "partition_str", partition_str.data()));

    return partition_str;
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
    multi_array<signed char, 2> dep_int({n_reactions, d});

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
    int id_r;
    NETCDF_CHECK(nc_inq_dimid(ncid, "r_out", &id_r));

    size_t r_t;
    char tmp[NC_MAX_NAME + 1];
    NETCDF_CHECK(nc_inq_dim(ncid, id_r, tmp, &r_t));

    std::array<Index, 2> r_out = {(Index) r_t, (Index) r_t};
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

std::vector<std::vector<double>> ReadHelpers::ReadPropensity(int ncid, const Index n_reactions)
{
    std::vector<std::vector<double>> result(n_reactions);

    for (Index mu = 0; mu < n_reactions; ++mu)
    {
        // read dimension dx_{mu}
        int id_dx_dep;
        NETCDF_CHECK(nc_inq_dimid(ncid, ("dx_" + std::to_string(mu)).c_str(), &id_dx_dep));

        size_t dx_dep_t;
        char tmp[NC_MAX_NAME + 1];
        NETCDF_CHECK(nc_inq_dim(ncid, id_dx_dep, tmp, &dx_dep_t));

        result[mu].resize(dx_dep_t);

        // read propensity
        int id_propensity;

        NETCDF_CHECK(nc_inq_varid(ncid, ("propensity_" + std::to_string(mu)).c_str(), &id_propensity));
        NETCDF_CHECK(nc_get_var_double(ncid, id_propensity, result[mu].data()));
    }

    return result;
}

cme_node* ReadHelpers::ReadNode(int ncid, const std::string id, cme_internal_node * const parent_node, const Index r_in)
{
    int retval0, retval1;
    int grp_ncid0, grp_ncid1;
    grid_parms grid = ReadGridParms(ncid);
    Index n_basisfunctions = ReadNBasisfunctions(ncid);
    cme_node *child_node;

    retval0 = nc_inq_ncid(ncid, (id + "0").c_str(), &grp_ncid0);
    retval1 = nc_inq_ncid(ncid, (id + "1").c_str(), &grp_ncid1);

    if (retval0 or retval1) // NOTE: if retval0 is 1, then retval1 has also to be 1, due to the binary tree structure of the netCDF file
    {
        child_node = new cme_external_node(id, parent_node, grid, r_in, n_basisfunctions);
        child_node->Initialize(ncid);
    }
    else
    {
        std::array<Index, 2> r_out = ReadRankOut(ncid);
        child_node = new cme_internal_node(id, parent_node, grid, r_in, r_out, n_basisfunctions);
        child_node->Initialize(ncid);
        child_node->child[0] = ReadNode(grp_ncid0, id + "0", (cme_internal_node *)child_node, ((cme_internal_node *) child_node)->RankOut()[0]);
        child_node->child[1] = ReadNode(grp_ncid1, id + "1", (cme_internal_node *)child_node, ((cme_internal_node *)child_node)->RankOut()[1]);
    }
    return child_node;
}

// TODO: Rename A_bar and B_bar
// TODO: Reduce number of loops to 5
void CalculateAB_bar(cme_node const * const child_node_init, multi_array<double, 3> &A_bar, multi_array<double, 3> &B_bar, const blas_ops &blas)
{
    std::fill(std::begin(A_bar), std::end(A_bar), 0.0);
    std::fill(std::begin(B_bar), std::end(B_bar), 0.0);

    if (child_node_init->IsExternal())
    {
        cme_external_node *child_node = (cme_external_node *) child_node_init;
        multi_array<double, 2> X_shift({child_node->grid.dx, child_node->RankIn()});
        multi_array<double, 2> tmp_A_bar({child_node->RankIn(), child_node->RankIn()});
        multi_array<double, 2> tmp_B_bar({child_node->RankIn(), child_node->RankIn()});
        multi_array<double, 1> weight({child_node->grid.dx});

        for (Index mu = 0; mu < child_node->grid.n_reactions; ++mu)
        {
            Matrix::ShiftRows<-1>(X_shift, child_node->X, child_node->grid, mu);

            std::vector<Index> vec_index(child_node->grid.d, 0);

            for (Index alpha = 0; alpha < child_node->grid.dx; ++alpha)
            {
                Index alpha_dep = IndexFunction::VecIndexToDepCombIndex(std::begin(vec_index), std::begin(child_node->grid.n_dep[mu]), std::begin(child_node->grid.idx_dep[mu]), std::end(child_node->grid.idx_dep[mu]));

                weight(alpha) = child_node->external_coefficients.propensity[mu][alpha_dep] * child_node->grid.h_mult;
                IndexFunction::IncrVecIndex(std::begin(child_node->grid.n), std::begin(vec_index), std::end(vec_index));
            }
            coeff(X_shift, child_node->X, weight, tmp_A_bar, blas);
            coeff(child_node->X, child_node->X, weight, tmp_B_bar, blas);

            for (Index i = 0; i < child_node->RankIn(); ++i)
            {
                for (Index j = 0; j < child_node->RankIn(); ++j)
                {
                    A_bar(mu, i, j) = tmp_A_bar(i, j);
                    B_bar(mu, i, j) = tmp_B_bar(i, j);
                }
            }
        }
    }
    else
    {
        cme_internal_node *child_node = (cme_internal_node *)child_node_init;
        multi_array<double, 3> A_bar_child0({child_node->grid.n_reactions, child_node->RankOut()[0], child_node->RankOut()[0]});
        multi_array<double, 3> B_bar_child0({child_node->grid.n_reactions, child_node->RankOut()[0], child_node->RankOut()[0]});
        multi_array<double, 3> A_bar_child1({child_node->grid.n_reactions, child_node->RankOut()[1], child_node->RankOut()[1]});
        multi_array<double, 3> B_bar_child1({child_node->grid.n_reactions, child_node->RankOut()[1], child_node->RankOut()[1]});

        CalculateAB_bar(child_node->child[0], A_bar_child0, B_bar_child0, blas);
        CalculateAB_bar(child_node->child[1], A_bar_child1, B_bar_child1, blas);

        for (Index mu = 0; mu < child_node->grid.n_reactions; ++mu)
        {
            for (Index i0 = 0; i0 < child_node->RankOut()[0]; ++i0)
            {
                for (Index j0 = 0; j0 < child_node->RankOut()[0]; ++j0)
                {
                    for (Index i1 = 0; i1 < child_node->RankOut()[1]; ++i1)
                    {
                        for (Index j1 = 0; j1 < child_node->RankOut()[1]; ++j1)
                        {
                            for (Index i = 0; i < child_node->RankIn(); ++i)
                            {
                                for (Index j = 0; j < child_node->RankIn(); ++j)
                                {
                                    A_bar(mu, i, j) += child_node->Q(i0, i1, i) * child_node->Q(j0, j1, j) * A_bar_child0(mu, i0, j0) * A_bar_child1(mu, i1, j1);
                                    B_bar(mu, i, j) += child_node->Q(i0, i1, i) * child_node->Q(j0, j1, j) * B_bar_child0(mu, i0, j0) * B_bar_child1(mu, i1, j1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void cme_external_node::CalculateCD()
{
    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        for (Index alpha = 0; alpha < grid.dx_dep[mu]; ++alpha)
        {
            for (Index i = 0; i < RankIn(); ++i)
            {
                for (Index j = 0; j < RankIn(); ++j)
                {
                    external_coefficients.C[mu](alpha, i, j) = external_coefficients.propensity[mu][alpha] * coefficients.A(mu, i, j);
                    external_coefficients.D[mu](alpha, i, j) = external_coefficients.propensity[mu][alpha] * coefficients.B(mu, i, j);
                }
            }
        }
    }
}

// TODO: Make this more efficient by using coefficients A and B instead of C and D and multiply the propensity in a separate loop

multi_array<double, 2> CalculateKDot(const multi_array<double, 2> &K, const cme_external_node* const node)
{
    get_time::start("CalculateKDot");
    multi_array<double, 2> prod_KC(K.shape());
    multi_array<double, 2> prod_KC_shift(K.shape());
    multi_array<double, 2> prod_KD(K.shape());
    multi_array<double, 2> K_dot(K.shape());
    set_zero(K_dot);

    std::vector<Index> vec_index(node->grid.d);

    for (Index mu = 0; mu < node->grid.n_reactions; ++mu)
    {
        set_zero(prod_KC);
        set_zero(prod_KD);
        std::fill(std::begin(vec_index), std::end(vec_index), 0);

#ifdef __OPENMP__
#pragma omp parallel firstprivate(vec_index)
#endif
        {
            Index alpha;
#ifdef __OPENMP__
            Index chunk_size = IndexFunction::SetVecIndex(std::begin(vec_index), std::end(vec_index), std::begin(grid.n), grid.dx);
#endif

#ifdef __OPENMP__
#pragma omp for schedule(static, chunk_size)
#endif
            for (Index i = 0; i < node->grid.dx; ++i)
            {
                alpha = IndexFunction::VecIndexToDepCombIndex(std::begin(vec_index), std::begin(node->grid.n_dep[mu]), std::begin(node->grid.idx_dep[mu]), std::end(node->grid.idx_dep[mu]));
                IndexFunction::IncrVecIndex(std::begin(node->grid.n), std::begin(vec_index), std::end(vec_index));

                // Calculate matrix-vector multiplication of C2*K and D2*K
                for (Index j = 0; j < node->RankIn(); ++j)
                {
                    for (Index l = 0; l < node->RankIn(); ++l)
                    {
                        prod_KC(i, j) += K(i, l) * node->external_coefficients.C[mu](alpha, j, l);
                        prod_KD(i, j) += K(i, l) * node->external_coefficients.D[mu](alpha, j, l);
                    }
                }
            }
        }
        // Shift prod_KC
        Matrix::ShiftRows<1>(prod_KC_shift, prod_KC, node->grid, mu);

        // Calculate k_dot = shift(C * K) - D * K
        K_dot += prod_KC_shift;
        K_dot -= prod_KD;
    }
    get_time::stop("CalculateKDot");
    return K_dot;
}

void cme_node::CalculateS(const double tau)
{
    multi_array<double, 2> S_dot(S.shape());
    set_zero(S_dot);

#ifdef __OPENMP__
#pragma omp parallel for collapse(2)
#endif
    for (Index i = 0; i < RankIn(); i++)
    {
        for (Index j = 0; j < RankIn(); j++)
        {
            for (Index k = 0; k < RankIn(); k++)
            {
                for (Index l = 0; l < RankIn(); l++)
                {
                    S_dot(i, j) += tau * (coefficients.E(i, j, k, l) - coefficients.F(i, j, k, l)) * S(k, l);
                }
            }
        }
    }
    S -= S_dot;
}

void cme_node::CalculateEF(const blas_ops& blas)
{
    std::fill(std::begin(coefficients.E), std::end(coefficients.E), 0.0);
    std::fill(std::begin(coefficients.F), std::end(coefficients.F), 0.0);

    multi_array<double, 3> A_bar({grid.n_reactions, RankIn(), RankIn()});
    multi_array<double, 3> B_bar({grid.n_reactions, RankIn(), RankIn()});

    CalculateAB_bar(this, A_bar, B_bar, blas);

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        for (Index i = 0; i < RankIn(); ++i)
        {
            for (Index j = 0; j < RankIn(); ++j)
            {
                for (Index k = 0; k < RankIn(); ++k)
                {
                    for (Index l = 0; l < RankIn(); ++l)
                    {
                        coefficients.E(i, j, k, l) += A_bar(mu, i, k) * coefficients.A(mu, j, l);
                        coefficients.F(i, j, k, l) += B_bar(mu, i, k) * coefficients.B(mu, j, l);
                    }
                }
            }
        }
    }
}

// TODO: Reduce number of loops to 5
void cme_internal_node::CalculateGH(const blas_ops& blas)
{
    std::fill(std::begin(internal_coefficients.G), std::end(internal_coefficients.G), 0.0);
    std::fill(std::begin(internal_coefficients.H), std::end(internal_coefficients.H), 0.0);

    multi_array<double, 3> A_bar_child0({grid.n_reactions, RankOut()[0], RankOut()[0]});
    multi_array<double, 3> B_bar_child0({grid.n_reactions, RankOut()[0], RankOut()[0]});
    multi_array<double, 3> A_bar_child1({grid.n_reactions, RankOut()[1], RankOut()[1]});
    multi_array<double, 3> B_bar_child1({grid.n_reactions, RankOut()[1], RankOut()[1]});

    CalculateAB_bar(child[0], A_bar_child0, B_bar_child0, blas);
    CalculateAB_bar(child[1], A_bar_child1, B_bar_child1, blas);

    for (Index mu = 0; mu < grid.n_reactions; ++mu)
    {
        for (Index i = 0; i < RankIn(); ++i)
        {
            for (Index j = 0; j < RankIn(); ++j)
            {
                for (Index i0 = 0; i0 < RankOut()[0]; ++i0)
                {
                    for (Index j0 = 0; j0 < RankOut()[0]; ++j0)
                    {
                        for (Index i1 = 0; i1 < RankOut()[1]; ++i1)
                        {
                            for (Index j1 = 0; j1 < RankOut()[1]; ++j1)
                            {
                                internal_coefficients.G(i, i0, i1, j, j0, j1) += coefficients.A(mu, i, j) * A_bar_child0(mu, i0, j0) * A_bar_child1(mu, i1, j1);
                                internal_coefficients.H(i, i0, i1, j, j0, j1) += coefficients.B(mu, i, j) * B_bar_child0(mu, i0, j0) * B_bar_child1(mu, i1, j1);
                            }
                        }
                    }
                }
            }
        }
    }
}

// TODO: Reduce number of loops to 4
void cme_internal_node::CalculateQ(const double tau)
{
    multi_array<double, 3> Q_dot(Q.shape());
    std::fill(std::begin(Q_dot), std::end(Q_dot), 0.0);

    for (Index i = 0; i < RankIn(); ++i)
    {
        for (Index j = 0; j < RankIn(); ++j)
        {
            for (Index i0 = 0; i0 < RankOut()[0]; ++i0)
            {
                for (Index j0 = 0; j0 < RankOut()[0]; ++j0)
                {
                    for (Index i1 = 0; i1 < RankOut()[1]; ++i1)
                    {
                        for (Index j1 = 0; j1 < RankOut()[1]; ++j1)
                        {
                            Q_dot(i0, i1, i) += tau * (internal_coefficients.G(i, i0, i1, j, j0, j1) - internal_coefficients.H(i, i0, i1, j, j0, j1)) * Q(j0, j1, j);
                        }
                    }
                }
            }
        }
    }
    Q += Q_dot;
}