#include "tree_class.hpp"

template<>
void internal_node<double>::Initialize(int ncid)
{
    int retval;

    // read Q
    int id_Q;
    if ((retval = nc_inq_varid(ncid, "Q", &id_Q)))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_Q, Q.data())))
        NETCDF_ERROR(retval);

    return;
}

void cme_internal_node::Initialize(int ncid)
{
    internal_node::Initialize(ncid);

    // read propensity
    coefficients.propensity = ReadHelpers::ReadPropensity(ncid, grid.n_reactions);
    return;
}

template<>
void external_node<double>::Initialize(int ncid)
{
    int retval;

    // read X
    int id_X;
    if ((retval = nc_inq_varid(ncid, "X", &id_X)))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_X, X.data())))
        NETCDF_ERROR(retval);

    return;
}

void cme_external_node::Initialize(int ncid)
{
    external_node::Initialize(ncid);

    // read propensity
    coefficients.propensity = ReadHelpers::ReadPropensity(ncid, grid.n_reactions);
    return;
}

void cme_lr_tree::Read(std::string fn)
{
    int ncid, retval;

    if ((retval = nc_open(fn.c_str(), NC_NOWRITE, &ncid)))
        NETCDF_ERROR(retval);

    grid_parms grid = ReadHelpers::ReadGridParms(ncid);
    std::array<Index, 2> r_out = ReadHelpers::ReadRankOut(ncid);
    root = new cme_internal_node("root", nullptr, grid, 1, r_out, 1);
    // root->Initialize(ncid);

    int grp_ncid;

    if ((retval = nc_inq_ncid(ncid, "0", &grp_ncid)))
        NETCDF_ERROR(retval);
    root->child[0] = ReadHelpers::ReadNode(grp_ncid, "0", root, root->RankOut()[0]);

    if ((retval = nc_inq_ncid(ncid, "1", &grp_ncid)))
        NETCDF_ERROR(retval);
    root->child[1] = ReadHelpers::ReadNode(grp_ncid, "1", root, root->RankOut()[1]);

    if ((retval = nc_close(ncid)))
        NETCDF_ERROR(retval);
    return;
}

void cme_lr_tree::PrintHelper(node* node)
{
    if (node->IsExternal())
    {
        cout << "external_node, id: " << node->id << ", X.shape(): (" << ((cme_external_node*) node)->X.shape()[0] << "," << ((cme_external_node*) node)->X.shape()[1] << ")" << endl;
    }
    else
    {
        cout << "internal_node, id: " << node->id << ", rank_out: (" << ((cme_internal_node *)node)->RankOut()[0] << "," << ((cme_internal_node *)node)->RankOut()[1] << ")" << endl;
        cme_lr_tree::PrintHelper(node->child[0]);
        cme_lr_tree::PrintHelper(node->child[1]);
    }
}

void cme_lr_tree::Print()
{
    cme_lr_tree::PrintHelper(root);
    return;
}

void cme_lr_tree::OrthogonalizeHelper(cme_internal_node *node)
{
    blas_ops blas;
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

        OrthogonalizeHelper(node_right);

        ip0 = inner_product_from_const_weight(node_left->grid.h_mult, node_left->grid.dx);
        ip1 = inner_product_from_const_weight(1.0, prod(node_right->RankOut()));

        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }
    else if (node->child[0]->IsInternal() and node->child[1]->IsExternal())
    {
        cme_internal_node *node_left = (cme_internal_node *)node->child[0];
        cme_external_node *node_right = (cme_external_node *)node->child[1];

        OrthogonalizeHelper(node_left);

        ip0 = inner_product_from_const_weight(1.0, prod(node_left->RankOut()));
        ip1 = inner_product_from_const_weight(node_right->grid.h_mult, node_right->grid.dx);

        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }
    else
    {
        cme_internal_node *node_left = (cme_internal_node *)node->child[0];
        cme_internal_node *node_right = (cme_internal_node *)node->child[1];

        OrthogonalizeHelper(node_left);
        OrthogonalizeHelper(node_right);

        ip0 = inner_product_from_const_weight(1.0, prod(node_left->RankOut()));
        ip1 = inner_product_from_const_weight(1.0, prod(node_right->RankOut()));
        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }

    for (int j = node->n_basisfunctions; j < node->RankOut()[1]; ++j)
    {
        for (int i = 0; i < node->RankOut()[0]; ++i)
        {
            R0(i, j) = 0.0;
            R1(i, j) = 0.0;
        }
    }

    for (Index k = 0; k < node->RankIn(); ++k)
    {
        for (Index j = 0; j < node->RankOut()[1]; ++j)
        {
            for (Index i = 0; i < node->RankOut()[0]; ++i)
            {
                node->S(i, j) = node->Q(i, j, k);
            }
        }
        blas.matmul_transb(R0, R1, node->S);
        for (Index j = 0; j < node->RankOut()[1]; ++j)
        {
            for (Index i = 0; i < node->RankOut()[0]; ++i)
            {
                node->Q(i, j, k) = node->S(i, j);
            }
        }
    }
}

void cme_lr_tree::Orthogonalize()
{
    OrthogonalizeHelper(root);
};

grid_parms ReadHelpers::ReadGridParms(int ncid)
{
    int retval;

    // read dimensions
    int id_d, id_n_reactions;
    if ((retval = nc_inq_dimid(ncid, "d", &id_d)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dimid(ncid, "n_reactions", &id_n_reactions)))
        NETCDF_ERROR(retval);

    size_t d_t, n_reactions_t;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_d, tmp, &d_t)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dim(ncid, id_n_reactions, tmp, &n_reactions_t)))
        NETCDF_ERROR(retval);

    Index d = (Index)d_t;
    Index n_reactions = (Index)n_reactions_t;

    // read variables
    int id_n, id_binsize, id_liml, id_dep, id_nu;

    if ((retval = nc_inq_varid(ncid, "n", &id_n)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "binsize", &id_binsize)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "liml", &id_liml)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "dep", &id_dep)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "nu", &id_nu)))
        NETCDF_ERROR(retval);

    grid_parms grid(d, n_reactions);
    multi_array<signed char, 2> dep_int({n_reactions, d});

    if ((retval = nc_get_var_long(ncid, id_n, grid.n.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_long(ncid, id_binsize, grid.binsize.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_liml, grid.liml.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_schar(ncid, id_dep, dep_int.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_long(ncid, id_dep, grid.nu.data())))
        NETCDF_ERROR(retval);

    std::copy(dep_int.begin(), dep_int.end(), grid.dep.begin());
    grid.Initialize();

    return grid;
}

// NOTE: Input .netCDF file is required to store a scalar `r_out`
// (which guarantees that the rank is equal for both childs)
std::array<Index, 2> ReadHelpers::ReadRankOut(int ncid)
{
    int retval;

    // read rank
    int id_r;
    if ((retval = nc_inq_dimid(ncid, "r_out", &id_r)))
        NETCDF_ERROR(retval);

    size_t r_t;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_r, tmp, &r_t)))
        NETCDF_ERROR(retval);

    std::array<Index, 2> r_out = {(Index) r_t, (Index) r_t};
    return r_out;
}

Index ReadHelpers::ReadNBasisfunctions(int ncid)
{
    int retval;

    // read number of basisfunctions
    int id_n_basisfunctions;
    if ((retval = nc_inq_dimid(ncid, "n_basisfunctions", &id_n_basisfunctions)))
        NETCDF_ERROR(retval);

    size_t n_basisfunctions_t;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_n_basisfunctions, tmp, &n_basisfunctions_t)))
        NETCDF_ERROR(retval);

    return (Index)n_basisfunctions_t;
}

std::vector<std::vector<double>> ReadHelpers::ReadPropensity(int ncid, const Index n_reactions)
{
    int retval;
    std::vector<std::vector<double>> result(n_reactions);

    for (Index mu = 0; mu < n_reactions; ++mu)
    {
        // read dimension dx_{mu}
        int id_dx_mu;
        if ((retval = nc_inq_dimid(ncid, ("dx_" + std::to_string(mu)).c_str(), &id_dx_mu)))
            NETCDF_ERROR(retval);

        size_t dx_mu_t;
        char tmp[NC_MAX_NAME + 1];
        if ((retval = nc_inq_dim(ncid, id_dx_mu, tmp, &dx_mu_t)))
            NETCDF_ERROR(retval);

        result[mu].resize(dx_mu_t);

        // read propensity
        int id_propensity;

        if ((retval = nc_inq_varid(ncid, ("propensity_" + std::to_string(mu)).c_str(), &id_propensity)))
            NETCDF_ERROR(retval);

        if ((retval = nc_get_var_double(ncid, id_propensity, result[mu].data())))
            NETCDF_ERROR(retval);
    }

    return result;
}

node* ReadHelpers::ReadNode(int ncid, std::string id, cme_internal_node *parent_node, Index r_in)
{
    int retval0, retval1;
    int grp_ncid0, grp_ncid1;
    grid_parms grid = ReadGridParms(ncid);
    Index n_basisfunctions = ReadNBasisfunctions(ncid);
    node *child_node;

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