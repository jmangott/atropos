#ifndef NETCDF_ERROR_HPP
#define NETCDF_ERROR_HPP

#define NETCDF_CHECK(e)                                                     \
    {                                                                       \
        int res = e;                                                        \
        if(res != NC_NOERR)                                                 \
        {                                                                   \
            cout << "NetCDF Error " << __FILE__ << ":" << __LINE__ << ": "  \
             << nc_strerror(res) << endl;                                   \
            exit(1);                                                        \
        }                                                                   \
    }

#define NETCDF_ERROR(e)                                                         \
    {                                                                           \
        printf("NetCDF Error %s:%i: %s\n", __FILE__, __LINE__, nc_strerror(e)); \
        exit(1);                                                                \
    }

#endif