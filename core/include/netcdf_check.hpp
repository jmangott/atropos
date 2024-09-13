#ifndef NETCDF_CHECK_HPP
#define NETCDF_CHECK_HPP

#include <cstdlib>

#define NETCDF_CHECK(e)                                                                \
    {                                                                                  \
        int res = e;                                                                   \
        if (res != NC_NOERR) {                                                         \
            std::cout << "NetCDF Error " << __FILE__ << ":" << __LINE__ << ": "        \
                      << nc_strerror(res) << std::endl;                                \
            exit(EXIT_FAILURE);                                                        \
        }                                                                              \
    }

#endif
