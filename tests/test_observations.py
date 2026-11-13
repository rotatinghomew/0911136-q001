import unittest
from datetime import datetime, timezone

from rainnet.observations import Reading, valid_readings


class ObservationTests(unittest.TestCase):
    def test_quality_and_missing_values_are_not_measurements(self):
        values = [
            Reading("s1", "v1", datetime.now(timezone.utc), 2.0),
            Reading("s2", "v1", datetime.now(timezone.utc), 9.0, "suspect"),
            Reading("s3", "v2", datetime.now(timezone.utc), None),
        ]
        self.assertEqual(["s1"], [item.station_id for item in valid_readings(values)])
