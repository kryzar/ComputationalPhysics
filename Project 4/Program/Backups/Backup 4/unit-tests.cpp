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
#include <math.h>

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
        lattice system(4, 1.);
        REQUIRE(system.neighbors_spins_sum(0, 0) == 4);
        REQUIRE(system.neighbors_spins_sum(3, 3) == 4);
        REQUIRE(system.neighbors_spins_sum(1, 2) == 4);
        system.change(1, 1);
        REQUIRE(system.neighbors_spins_sum(0, 1) == 2);
        system.negativize();
        REQUIRE(system.neighbors_spins_sum(0, 0) == -4);
        REQUIRE(system.neighbors_spins_sum(3, 3) == -4);
    }
    
    SECTION("energy")
    {
        lattice system(4, 42.);
        REQUIRE(system.energy(0, 0) == -4*J);
        REQUIRE(system.energy(3, 1) == -4*J);
        REQUIRE(system.energy(1, 1) == -4*J);
        system.change(0, 0);
        REQUIRE(system.energy(1, 1) == -4*J);
        REQUIRE(system.energy(0, 1) == -2*J);
        system.change(1, 1);
        REQUIRE(system.energy(1, 0) == 0.);
        system.negativize();
        REQUIRE(system.energy(0, 0) == -4*J);
        system.change(1,2);
        REQUIRE(system.energy(1, 1) == -2*J);
        system.change(0, 1);
        system.change(2, 1);
        REQUIRE(system.energy(1, 1) == 2*J);
    }
    
    SECTION("energy_delta")
    {
        lattice system(10, 42.);
        system.change(1,1);
        REQUIRE(system.energy_delta(1, 1) == 8*J);
        REQUIRE(system.energy_delta(0, 0) == -8*J);
        REQUIRE(system.energy_delta(0, 1) == -4*J);
        system.negativize();
        REQUIRE(system.energy_delta(0, 0) == -8*J);
        system.change(0, 0);
        REQUIRE(system.energy_delta(0, 0) == 8*J);
        system.change(1, 9);
        REQUIRE(system.energy_delta(1, 0) == 0.);
        system.change(2, 0);
        REQUIRE(system.energy_delta(1, 0) == 4*J);
    }
    
    SECTION("exp_energy_delta")
    {
        lattice system(10, 42.);
        system.change(1,1);
        REQUIRE(system.exp_energy_delta(1, 1) == exp(-(8 * J) / (K * 42.)));
        REQUIRE(system.exp_energy_delta(0, 0) == exp((8 * J) / (K * 42.)));
        REQUIRE(system.exp_energy_delta(0, 1) == exp((4 * J) / (K * 42.)));
        system.negativize();
        REQUIRE(system.exp_energy_delta(0, 0) == exp((8 * J) / (K * 42.)));
        system.change(0, 0);
        REQUIRE(system.exp_energy_delta(0, 0) == exp(-(8 * J) / (K * 42.)));
        system.change(1, 9);
        REQUIRE(system.exp_energy_delta(1, 0) == 1.);
        system.change(2, 0);
        REQUIRE(system.exp_energy_delta(1, 0) == exp(-(4 * J) / (K * 42.)));
    }
}


