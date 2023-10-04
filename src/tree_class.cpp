#include "tree_class.hpp"

void root_node::InitializeNode(int ncid)
{
    int retval;

    // read S
    int id_S;
    if ((retval = nc_inq_varid(ncid, "S", &id_S)))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_S, S.data())))
        NETCDF_ERROR(retval);

    return;
};

void internal_node::InitializeNode(int ncid)
{
    int retval;

    // read Q and S
    int id_Q, id_S;
    if ((retval = nc_inq_varid(ncid, "Q", &id_Q)))
        NETCDF_ERROR(retval);
    if ((retval = nc_inq_varid(ncid, "S", &id_S)))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_Q, Q.data())))
        NETCDF_ERROR(retval);
    if ((retval = nc_get_var_double(ncid, id_S, S.data())))
        NETCDF_ERROR(retval);
    
    return;
};

void external_node::InitializeNode(int ncid)
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

struct cme_lr_tree::ReadTreeHelpers
{
    grid_parms ReadGridParms(int ncid)
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

    Index ReadRank(int ncid)
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

    node* ReadNode(int ncid, std::string id, node *parent_node)
    {
        int retval0, retval1;
        int grp_ncid0, grp_ncid1;
        grid_parms grid = ReadGridParms(ncid);
        node *child_node;

        retval0 = nc_inq_ncid(ncid, (id + "0").c_str(), &grp_ncid0);
        retval1 = nc_inq_ncid(ncid, (id + "1").c_str(), &grp_ncid1);

        if (retval0 or retval1) // NOTE: if retval0 is 1, then retval1 has also to be 1, due to the binary tree structure of the netCDF file
        {
            if (parent_node->IsRoot())
            {
                child_node = new cme_external_node(id, (cme_root_node *)parent_node, grid);
            }
            else
            {
                child_node = new cme_external_node(id, (cme_internal_node *)parent_node, grid);
            }
            child_node->InitializeNode(ncid);
        }
        else
        {
            Index r = ReadRank(ncid);
            if (parent_node->IsRoot())
            {
                child_node = new cme_internal_node(id, (cme_root_node *)parent_node, grid, r);
            }
            else
            {
                child_node = new cme_internal_node(id, (cme_internal_node *)parent_node, grid, r);
            }
            child_node->InitializeNode(ncid);
            child_node->left = ReadNode(grp_ncid0, id + "0", child_node);
            child_node->right = ReadNode(grp_ncid1, id + "1", child_node);
        }
        return child_node;
    }
};

void cme_lr_tree::ReadTree(std::string fn)
{
    int ncid, retval;

    if ((retval = nc_open(fn.c_str(), NC_NOWRITE, &ncid)))
        NETCDF_ERROR(retval);

    grid_parms grid = pReadTreeHelpers->ReadGridParms(ncid);
    Index r = pReadTreeHelpers->ReadRank(ncid);
    root = new cme_root_node(grid, r);

    int grp_ncid;

    if ((retval = nc_inq_ncid(ncid, "0", &grp_ncid)))
        NETCDF_ERROR(retval);
    root->left = pReadTreeHelpers->ReadNode(grp_ncid, "0", root);

    if ((retval = nc_inq_ncid(ncid, "1", &grp_ncid)))
        NETCDF_ERROR(retval);
    root->right = pReadTreeHelpers->ReadNode(grp_ncid, "1", root);

    if ((retval = nc_close(ncid)))
        NETCDF_ERROR(retval);
    return;
}

void cme_lr_tree::PrintTreeHelper(node* node)
{
    if (node->IsExternal())
    {
        cout << "external_node, id: " << node->id << "X.shape(): " << endl;
    }
    else
    {
        if (node->IsRoot())
        {
            cout << "root_node, id: " << node->id << endl;
        }
        else
        {
            cout << "internal_node, id: " << node->id << endl;
        }
        cme_lr_tree::PrintTreeHelper(node->left);
        cme_lr_tree::PrintTreeHelper(node->right);
    }
}

void cme_lr_tree::PrintTree()
{
    cme_lr_tree::PrintTreeHelper(root);
    return;
}