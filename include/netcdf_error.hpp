#ifndef NETCDF_ERROR_HPP
#define NETCDF_ERROR_HPP

#define NETCDF_ERROR(e)                                                         \
    {                                                                           \
        printf("NetCDF Error %s:%i: %s\n", __FILE__, __LINE__, nc_strerror(e)); \
        exit(1);                                                                \
    }

#endif