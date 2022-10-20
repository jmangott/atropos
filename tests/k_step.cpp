#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "index_functions.hpp"
#include "reactions.hpp"

TEST_CASE("K-Step", "[k_step]")
{
    SECTION("VecIndexToState")
    {
        multi_array<Index, 1> vec_index({10});
        multi_array<Index, 1> interval({10});
        multi_array<double, 1> limit({10});
        multi_array<double, 1> state_vec({10});
        multi_array<double, 1> comparison_vec({10});
        for (auto i = 0; i < 10; i++)
        {
            vec_index(i) = i;
            interval(i) = 11;
            limit(i) = 20.0;
            comparison_vec(i) = i * 2.0;
        }

        state_vec = VecIndexToState(vec_index, interval, limit);
        REQUIRE(bool(state_vec == comparison_vec));
    }
}