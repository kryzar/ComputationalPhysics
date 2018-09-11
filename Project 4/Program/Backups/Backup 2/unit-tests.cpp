//
//  unit-tests.cpp
//  Program
//
//  Created by Antoine Hugounet on 06/11/2017.
//  Copyright © 2017 Antoine Hugounet and Ethel Villeneuve. All rights reserved.
//

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "lattice.hpp"

using namespace std;

//  doesn't work for catch 2.0

int run_unittest(int argc, const char* argv[])
{
    return (Catch::Session().run(argc, argv));
}

TEST_CASE("Unit-tests on a few methods", "[lattice]")
{
    SECTION("neighbors_spins_sum")
    {
        lattice system(4, 1., 1.);
        REQUIRE(system.neighbors_spins_sum(0, 0) == 4);
        REQUIRE(system.neighbors_spins_sum(3, 3) == 4);
        REQUIRE(system.neighbors_spins_sum(1, 2) == 4);
        system.change(1, 1);
        REQUIRE(system.neighbors_spins_sum(0, 1) == 2);
        system.negativize();
        REQUIRE(system.neighbors_spins_sum(0, 0) == -4);
        REQUIRE(system.neighbors_spins_sum(3, 3) == -4);
    }
    
    SECTION("energie")
    {
        lattice system(4, 3., 1.);
        REQUIRE(system.energie(0, 0) == -12.);
        REQUIRE(system.energie(3, 1) == -12.);
        REQUIRE(system.energie(1, 1) == -12.);
        system.change(0, 0);
        REQUIRE(system.energie(1, 1) == -12.);
        REQUIRE(system.energie(0, 1) == -6.);
        system.change(1, 1);
        REQUIRE(system.energie(1, 0) == 0.);
        system.negativize();
        REQUIRE(system.energie(0, 0) == -12.);
        system.change(1,2);
        REQUIRE(system.energie(1, 1) == -6.);
        system.change(0, 1);
        system.change(2, 1);
        REQUIRE(system.energie(1, 1) == 6.);
    }
    
    SECTION("energie_delta")
    {
        lattice system(10, 4., 1.);
        system.change(1,1);
        REQUIRE(system.energie_delta(1, 1) == 32.);
        REQUIRE(system.energie_delta(0, 0) == -32.);
        REQUIRE(system.energie_delta(0, 1) == -16.);
        system.negativize();
        REQUIRE(system.energie_delta(0, 0) == -32.);
        system.change(0, 0);
        REQUIRE(system.energie_delta(0, 0) == 32.);
        system.change(1, 9);
        REQUIRE(system.energie_delta(1, 0) == 0.);
        system.change(2, 0);
        REQUIRE(system.energie_delta(1, 0) == 16.);
    }
}


