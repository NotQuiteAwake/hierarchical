#include "doctest/doctest.h"

#include "Distribution.hpp"

namespace sim {

namespace dist {

TEST_CASE("Test Distribution") {
    double mean = 5;
    double spread = 5;
    int n = int(1e5);
    double eps = 0.01;
    
    // do some Monte Carlo...
    SUBCASE("Test MakeNormalMass") {
        double one_sig = mean + spread * 1;
        double two_sig = mean + spread * 2;
        int one_sig_cnt = 0;
        int two_sig_cnt = 0;
        double one_sig_exp = 0.341;
        double two_sig_exp = 0.477;

        std::vector<Particle> par_list = MakeNormalMass(mean, spread, n);

        for (const Particle& par : par_list) {
            CHECK(par.GetMass() >= GetMassCutoff());
            double mass = par.GetMass();
            if (mass >= mean && mass <= two_sig) {
                two_sig_cnt++;
                if (mass <= one_sig) {
                    one_sig_cnt++;
                }
            }
        }

        CHECK(double(two_sig_cnt) / n
                == doctest::Approx(two_sig_exp).epsilon(eps));
        CHECK(double(one_sig_cnt) / n
                == doctest::Approx(one_sig_exp).epsilon(eps));

    }

    SUBCASE("Test MakeUniformMass") {
        double l = 0, r = 10;
        int intervals = int(r - l + 0.5);
        double exp = 1.0 / intervals;
        double statistics[int(r)];

        std::vector<Particle> par_list = MakeUniformMass(mean, spread, n);
        for (const Particle& par : par_list) {
            CHECK(par.GetMass() >= GetMassCutoff());
            CHECK(par.GetMass() <= r);
            int val_interval = int(par.GetMass());
            statistics[val_interval]++;
        }
        
        for (int i = 0; i < intervals; i++) {
            CHECK(exp == doctest::Approx(1.0 * statistics[i] / n).epsilon(eps));
        }
    }

    SUBCASE("Test SetSeed") {
        dist::SetSeed(1);
        const Particle par1 = MakeNormalMass(mean, spread, 1)[0];

        dist::SetSeed(1);
        const Particle par2 = MakeNormalMass(mean, spread, 1)[0];

        dist::SetSeed(2);
        const Particle par3 = MakeNormalMass(mean, spread, 1)[0];

        CHECK(par1.GetMass() == par2.GetMass());
        CHECK(par1.GetMass() != par3.GetMass());
    }

}
    

}

}

