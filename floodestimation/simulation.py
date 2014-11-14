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
        # Number of simulations
        n_sim = 500
        # Number of donor catchments in pooling group
        n_d = len(self.rec_lengths)
        # Total number of years in pooling group
        n_p = sum(self.rec_lengths)
        # Donor weights
        weight_d = np.asarray(self.rec_lengths, dtype=np.float64) / n_p

        # V2 observed:
        v2 = []
        # donor_t2s = np.empty(n_donors)
        # donor_t3s = np.empty(n_donors)
        # for donor_index, donor in enumerate(self.growth_curve_analysis.donor_catchments):
        # record = np.array([record.flow for record in donor.amax_records if record.flag == 0])
        #     record /= np.median(record)
        #     l1, l2, donor_t3s[donor_index] = lm.samlmu(record, nmom=3)
        #     donor_t2s[donor_index] = l2 / l1
        # v2_obs = math.sqrt(sum(donor_weights * ((donor_t2s - t2_pool) ** 2 + (donor_t3s - t3_pool) ** 2)))
        # print("v2 obs: {}".format(v2_obs))


        # Generated records using pooling group kappa distribution for all donors at once
        kap_para = lm.pelkap(self.lmom_p)
        kappa_distr = distr.Kappa(loc=kap_para[0], scale=kap_para[1], k=kap_para[2], h=kap_para[3])
        record_sim_all = kappa_distr.ppf(np.random.random(n_p * n_sim))

        record_start = 0
        # Simulated test statistic: variability V2
        variabilities_sim = np.empty(n_sim)

        # Loop through all simulations
        for i_sim in range(n_sim):
            # Second and third sample L-moment ratios for the **simulated** donor record
            t2s_d = np.empty(n_d)
            t3s_d = np.empty(n_d)

            # Loop through all donors
            for i_d in range(n_d):
                # Simulated donor record taken from large simulated record. All donors use the same distribution as the
                # pooling group because null hypothesis is that all donors have same distribution as the pooling group.
                record_d_sim = record_sim_all[record_start: record_start + self.rec_lengths[i_d]]
                # Sample L-moment ratios from simulated record
                l1, l2, t3 = lm.samlmu(record_d_sim, nmom=3)
                t2s_d[i_d] = l2 / l1
                t3s_d[i_d] = t3
                # Next time, take the next sequence of simulated records
                record_start += self.rec_lengths[i_d]

            # The test statistic is V2, root of squared errors of second and third L-moment ratios, weighted by the
            # record length of each donor.
            variabilities_sim[i_sim] = math.sqrt(
                np.sum(weight_d * ((t2s_d - self.lmom_p[1]) ** 2 + (t3s_d - self.lmom_p[2]) ** 2)))

        return np.median(variabilities_sim), np.std(variabilities_sim)



