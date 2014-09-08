import unittest
from floodestimation.entities import Catchment

class TestCatchmentObject(unittest.TestCase):
    def test_catchment_location_watercourse(self):
        catchment = Catchment("Aberdeen", "River Dee")
        self.assertEqual(catchment.location, "Aberdeen")
        self.assertEqual(catchment.watercourse, "River Dee")

    def test_catchment_location_only(self):
        catchment = Catchment("Aberdeen")
        self.assertEqual(catchment.location, "Aberdeen")
        self.assertEqual(catchment.watercourse, None)

    def test_catchment_distance_no_country(self):
        catchment_1 = Catchment("Aberdeen", "River Dee")
        catchment_1.descriptors.centroid_ngr = (0, 0)

        catchment_2 = Catchment("Dundee", "River Tay")
        catchment_2.descriptors.centroid_ngr = (3, 4)

        self.assertEqual(catchment_1.distance_to(catchment_2), 5)

    def test_catchment_distance_same_country(self):
        catchment_1 = Catchment("Aberdeen", "River Dee")
        catchment_1.descriptors.centroid_ngr = (0, 0)
        catchment_1.country = 'gb'

        catchment_2 = Catchment("Dundee", "River Tay")
        catchment_2.descriptors.centroid_ngr = (3, 4)
        catchment_2.country = 'gb'

        self.assertEqual(catchment_1.distance_to(catchment_2), 5)

    def test_catchment_distance_different_country(self):
        catchment_1 = Catchment("Aberdeen", "River Dee")
        catchment_1.descriptors.centroid_ngr = (0, 0)

        catchment_2 = Catchment("Belfast")
        catchment_2.descriptors.centroid_ngr = (3, 4)
        catchment_2.country = 'ni'

        self.assertEqual(catchment_1.distance_to(catchment_2), float('inf'))

    def test_descriptor_urbext(self):
        catchment_1 = Catchment("Aberdeen", "River Dee")
        catchment_1.descriptors.urbext = 0.2

        catchment_2 = Catchment("Dundee", "River Tay")
        catchment_2.descriptors.urbext = 0.5

        self.assertEqual(0.2, catchment_1.descriptors.urbext)
        self.assertEqual(0.5, catchment_2.descriptors.urbext)