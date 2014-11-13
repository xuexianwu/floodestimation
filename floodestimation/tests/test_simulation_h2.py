# -*- coding: utf-8 -*-

import unittest
import os
from urllib.request import pathname2url
from datetime import datetime

from floodestimation.simulation import H2Simulation
from floodestimation.collections import CatchmentCollections
from floodestimation.loaders import load_catchment
from floodestimation.analysis import GrowthCurveAnalysis
from floodestimation import settings
from floodestimation import db


class TestH2Simulation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        settings.OPEN_HYDROLOGY_JSON_URL = 'file:' + pathname2url(os.path.abspath('./floodestimation/fehdata_test.json'))
        cls.db_session = db.Session()
        gauged_catchments = CatchmentCollections(cls.db_session)
        catchment = load_catchment('floodestimation/tests/data/37017.CD3')


        analysis = GrowthCurveAnalysis(catchment, gauged_catchments)
        analysis.find_donor_catchments()
        cls.rec_lengths = [len(donor.amax_records) for donor in analysis.donor_catchments]
        t2, t3 = analysis._var_and_skew(analysis.donor_catchments)
        cls.lmom_p = [1, t2, t3, 0]  # For just now assume t4=0, need to check how to calculate pooled t4
        print("Pooled L-moments: {}".format(cls.lmom_p))

    def setUp(self):
        self.start_time = datetime.now()

    def tearDown(self):
        duration = datetime.now() - self.start_time
        print("Test ran for {} s.".format(duration.total_seconds()))
        self.db_session.rollback()

    @classmethod
    def tearDownClass(cls):
        cls.db_session.close()

    def test_h2(self):
        sim = H2Simulation(self.rec_lengths, self.lmom_p)
        sim.simulated_mean_dev()