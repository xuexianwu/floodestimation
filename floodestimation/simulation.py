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

    def __init__(self, rec_lengths, lmom_p):
        self.rec_lengths = rec_lengths
        try:
            if len(lmom_p) == 4:
                self.lmom_p = lmom_p
            else:
                raise ValueError("Pooled L-moments must be a list of 4 values.")
        except TypeError:
            raise ValueError("Pooled L-moments must be a list of 4 values.")

    def simulated_mean_dev(self):
        n = 500
        n_donors = len(self.rec_lengths)
        print("Donors: {}".format(n_donors))
        total_years = sum(self.rec_lengths)
        print("Pooled records: {}".format(total_years))
        donor_weights = np.asarray(self.rec_lengths, dtype=np.float32) / total_years
        print("Donor weights: {}".format(donor_weights))

        # V2 observed:
        v2 = []
        # donor_t2s = np.empty(n_donors)
        # donor_t3s = np.empty(n_donors)
        # for donor_index, donor in enumerate(self.growth_curve_analysis.donor_catchments):
        #     record = np.array([record.flow for record in donor.amax_records if record.flag == 0])
        #     record /= np.median(record)
        #     l1, l2, donor_t3s[donor_index] = lm.samlmu(record, nmom=3)
        #     donor_t2s[donor_index] = l2 / l1
        # v2_obs = math.sqrt(sum(donor_weights * ((donor_t2s - t2_pool) ** 2 + (donor_t3s - t3_pool) ** 2)))
        # print("v2 obs: {}".format(v2_obs))


        # Generated records using pooling group kappa distribution for all donors at once
        kap_para = lm.pelkap(self.lmom_p)
        print("Kappa paras: {}".format(kap_para))
        kappa_distr = distr.Kappa(loc=kap_para[0], scale=kap_para[1], k=kap_para[2], h=kap_para[3])
        total_sim_record = kappa_distr.ppf(np.random.random(total_years * n))
        print("Simulated records: {}".format(len(total_sim_record)))

        record_start = 0
        simulated_v2s = np.empty(n)
        for sim_index in range(n):
            donor_t2s = np.empty(n_donors)
            donor_t3s = np.empty(n_donors)

            for donor_index in range(n_donors):
                sim_record = total_sim_record[record_start : record_start + self.rec_lengths[donor_index]]
                # simulated sample L-moment ratios
                l1, l2, t3 = lm.samlmu(sim_record, nmom=3)
                donor_t2s[donor_index] = l2 / l1
                donor_t3s[donor_index] = t3
                record_start += self.rec_lengths[donor_index]

            simulated_v2s[sim_index] = np.sqrt(
                np.sum(donor_weights * ((donor_t2s - self.lmom_p[1]) ** 2 + (donor_t3s - self.lmom_p[2]) ** 2)))

        print("median: {}".format(np.median(simulated_v2s)))
        print("st dev: {}".format(np.std(simulated_v2s)))



