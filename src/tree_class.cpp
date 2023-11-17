#include "tree_class.hpp"

// template<>
// void root_node<double>::InitializeNode(int ncid)
// {
//     int retval;

//     // read S
//     int id_S;
//     if ((retval = nc_inq_varid(ncid, "S", &id_S)))
//         NETCDF_ERROR(retval);
//     if ((retval = nc_get_var_double(ncid, id_S, S.data())))
//         NETCDF_ERROR(retval);

//     return;
// };

template<>
void internal_node<double>::Initialize(int ncid)
{
    int retval;

    // read Q
    int id_Q;
    if ((retval = nc_inq_varid(ncid, "Q", &id_Q)))
        NETCDF_ERROR(retval);
    // if ((retval = nc_inq_varid(ncid, "S", &id_S)))
    //     NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_Q, Q.data())))
        NETCDF_ERROR(retval);
    // if ((retval = nc_get_var_double(ncid, id_S, S.data())))
    //     NETCDF_ERROR(retval);
    
    return;
};

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
};

void cme_lr_tree::Read(std::string fn)
{
    int ncid, retval;

    if ((retval = nc_open(fn.c_str(), NC_NOWRITE, &ncid)))
        NETCDF_ERROR(retval);

    grid_parms grid = ReadHelpers::ReadGridParms(ncid);
    Index r = ReadHelpers::ReadRank(ncid);
    root = new cme_internal_node("root", nullptr, grid, r, r);

    int grp_ncid;

    if ((retval = nc_inq_ncid(ncid, "0", &grp_ncid)))
        NETCDF_ERROR(retval);
    root->left = ReadHelpers::ReadNode(grp_ncid, "0", root);

    if ((retval = nc_inq_ncid(ncid, "1", &grp_ncid)))
        NETCDF_ERROR(retval);
    root->right = ReadHelpers::ReadNode(grp_ncid, "1", root);

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
        cout << "internal_node, id: " << node->id << ", rank: " << ((cme_internal_node*) node)->Rank() << endl;
        cme_lr_tree::PrintHelper(node->left);
        cme_lr_tree::PrintHelper(node->right);
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
    multi_array<double, 2> R0({node->Rank(), node->Rank()});
    multi_array<double, 2> R1({node->Rank(), node->Rank()});
    multi_array<double, 2> tmp({node->Rank(), node->Rank()});

    if (node->left->IsExternal() and node->right->IsExternal())
    {
        cme_external_node* node_left = (cme_external_node*) node->left;
        cme_external_node* node_right = (cme_external_node*) node->right;

        ip0 = inner_product_from_const_weight(double(node_left->grid.h_mult()), node_left->grid.dx());
        ip1 = inner_product_from_const_weight(double(node_right->grid.h_mult()), node_right->grid.dx());

        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }
    else if (node->left->IsExternal() and node->right->IsInternal())
    {
        cme_external_node *node_left = (cme_external_node *)node->left;
        cme_internal_node *node_right = (cme_internal_node *)node->right;

        OrthogonalizeHelper(node_right);

        ip0 = inner_product_from_const_weight(double(node_left->grid.h_mult()), node_left->grid.dx());
        ip1 = inner_product_from_const_weight(1.0, node_right->Rank() * node_right->Rank());

        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }
    else if (node->left->IsInternal() and node->right->IsExternal())
    {
        cme_internal_node *node_left = (cme_internal_node *)node->left;
        cme_external_node *node_right = (cme_external_node *)node->right;

        OrthogonalizeHelper(node_left);

        ip0 = inner_product_from_const_weight(1.0, node_left->Rank() * node_left->Rank());
        ip1 = inner_product_from_const_weight(double(node_right->grid.h_mult()), node_right->grid.dx());

        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }
    else
    {
        cme_internal_node *node_left = (cme_internal_node *)node->left;
        cme_internal_node *node_right = (cme_internal_node *)node->right;

        OrthogonalizeHelper(node_left);
        OrthogonalizeHelper(node_right);

        ip0 = inner_product_from_const_weight(1.0, node_left->Rank() * node_left->Rank());
        ip0 = inner_product_from_const_weight(1.0, node_right->Rank() * node_right->Rank());

        R0 = node_left->Orthogonalize(ip0, blas);
        R1 = node_right->Orthogonalize(ip1, blas);
    }

    for (int j = node->n_basisfunctions; j < node->Rank(); ++j)
    {
        for (int i = 0; i < node->Rank(); ++i)
        {
            R0(i, j) = 0.0;
            R1(i, j) = 0.0;
        }
    }

    blas.matmul_transb(R0, R1, node->S);

    for (Index k = 0; k < node->ParentRank(); ++k)
    {
        for (Index j = 0; j < node->Rank(); ++j)
        {
            for (Index i = 0; i < node->Rank(); ++i)
            {
                node->S(i, j) = node->Q(i, j, k);
            }
        }
        blas.matmul_transb(R0, R1, node->S);
        for (Index j = 0; j < node->Rank(); ++j)
        {
            for (Index i = 0; i < node->Rank(); ++i)
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
    int id_d, id_mu;
    if ((retval = nc_inq_dimid(ncid, "d", &id_d)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dimid(ncid, "mu", &id_mu)))
        NETCDF_ERROR(retval);

    size_t d_t, mu_t;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_d, tmp, &d_t)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_dim(ncid, id_mu, tmp, &mu_t)))
        NETCDF_ERROR(retval);

    Index d = (Index)d_t;
    Index mu = (Index)mu_t;

    // read variables
    int id_n, id_binsize, id_liml, id_dep;

    if ((retval = nc_inq_varid(ncid, "n", &id_n)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "binsize", &id_binsize)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "liml", &id_liml)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "dep", &id_dep)))
        NETCDF_ERROR(retval);

    grid_parms grid(d, mu);
    multi_array<signed char, 2> dep_int({mu, d});

    if ((retval = nc_get_var_long(ncid, id_n, grid.n.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_long(ncid, id_binsize, grid.binsize.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_liml, grid.liml.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_schar(ncid, id_dep, dep_int.data())))
        NETCDF_ERROR(retval);

    std::copy(dep_int.begin(), dep_int.end(), grid.dep.begin());

    return grid;
}

Index ReadHelpers::ReadRank(int ncid)
{
    int retval;

    // read rank
    int id_r;
    if ((retval = nc_inq_dimid(ncid, "r", &id_r)))
        NETCDF_ERROR(retval);

    size_t r_t;
    char tmp[NC_MAX_NAME + 1];
    if ((retval = nc_inq_dim(ncid, id_r, tmp, &r_t)))
        NETCDF_ERROR(retval);

    return (Index)r_t;
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

node* ReadHelpers::ReadNode(int ncid, std::string id, cme_internal_node *parent_node)
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
        child_node = new cme_external_node(id, parent_node, grid, n_basisfunctions);
        child_node->Initialize(ncid);
    }
    else
    {
        Index r = ReadRank(ncid);
        child_node = new cme_internal_node(id, parent_node, grid, r, n_basisfunctions);
        child_node->Initialize(ncid);
        child_node->left = ReadNode(grp_ncid0, id + "0", (cme_internal_node *)child_node);
        child_node->right = ReadNode(grp_ncid1, id + "1", (cme_internal_node *)child_node);
    }
    return child_node;
}