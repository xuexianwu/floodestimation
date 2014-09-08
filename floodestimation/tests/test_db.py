import unittest
from datetime import date

from floodestimation.db import Session, create_new_db
from floodestimation.entities import Catchment, AmaxRecord


class TestDatabase(unittest.TestCase):
    def setUp(self):
        create_new_db()
        self.db_session = Session()
        # TODO: assert something

    def test_add_catchment(self):
        catchment = Catchment(location="Aberdeen", watercourse="River Dee")
        self.db_session.add(catchment)
        self.db_session.commit()
        # TODO: assert something

    def test_add_catchment_with_amax(self):
        catchment = Catchment("Aberdeen", "River Dee")
        catchment.amax_records = [AmaxRecord(date(1999, 12, 31), 3.0, 0.5),
                                  AmaxRecord(date(2000, 12, 31), 2.0, 0.5),
                                  AmaxRecord(date(2001, 12, 31), 1.0, 0.5)]
        self.db_session.add(catchment)
        self.db_session.commit()
        # TODO: assert something