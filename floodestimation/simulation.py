# -*- coding: utf-8 -*-

# Copyright (c) 2014  Florenz A.P. Hollebrandse <f.a.p.hollebrandse@protonmail.ch>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import math
import lmoments3 as lm
import lmoments3.distr as distr


class H2Simulation(object):

    def __init__(self, growth_curve_analysis):
        self.growth_curve_analysis = growth_curve_analysis
        if not growth_curve_analysis.donor_catchments:
            growth_curve_analysis.donor_catchments = growth_curve_analysis.find_donors()
        self.donor_weights = np.array([float(len(donor.amax_records)) for donor in self.growth_curve_analysis.donor_catchments])
        self.donor_weights /= sum(self.donor_weights)

    def simulated_mean_dev(self):
        # Pooling group L-moments
        t2_pool, t3_pool = self.growth_curve_analysis._var_and_skew(self.growth_curve_analysis.donor_catchments)
        # Fitted kappa distribution function paras
        kap_para = lm.pelkap([1, t2_pool, t3_pool, 0])
        n = 500

        # V2 observed:
        v2 = []
        for donor_index, donor in enumerate(self.growth_curve_analysis.donor_catchments):
            record = np.array([record.flow for record in donor.amax_records if record.flag == 0])
            record /= np.median(record)
            l1, l2, t3 = lm.samlmu(record, nmom=3)
            t2 = l2 / l1
            v2.append(self.donor_weights[donor_index] * ((t2 - t2_pool) ** 2 + (t3 - t3_pool) ** 2))
        v2_obs = np.sqrt(sum(v2))
        print("v2 obs: {}".format(v2_obs))


        v2s = np.empty(n)
        total_years = sum([len(donor.amax_records) for donor in self.growth_curve_analysis.donor_catchments])
        # Generated records using pooling group kappa distribution for all donors at once
        kappa_distr = distr.Kappa(loc=kap_para[0], scale=kap_para[1], k=kap_para[2], h=kap_para[3])

        total_sim_record = kappa_distr.ppf(np.random.random(total_years * n))
        record_start = 0
        for sim_index in range(n):
            donor_t2s = np.empty(len(self.growth_curve_analysis.donor_catchments))
            donor_t3s = np.empty(len(self.growth_curve_analysis.donor_catchments))

            for donor_index, donor in enumerate(self.growth_curve_analysis.donor_catchments):
                record_length = len(donor.amax_records)
                sim_record = total_sim_record[record_start : record_start + record_length]
                # simulated sample L-moment ratios
                l1, l2, donor_t3s[donor_index] = lm.samlmu(sim_record, nmom=3)
                donor_t2s[donor_index] = l2 / l1
                record_start += record_length
            v2_sim = self.donor_weights * ((donor_t2s - t2_pool) ** 2 + (donor_t3s - t3_pool) ** 2)
            v2s[sim_index] = math.sqrt(sum(v2_sim))

        h2 = (v2_obs - np.median(v2s)) / np.std(v2s)
        print("H2: {}".format(h2))



