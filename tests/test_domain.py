"""验证现有领域摘要的业务边界。"""
import copy
import json
import unittest
from pathlib import Path

from rainnet.analysis import station_totals


class DomainTests(unittest.TestCase):
    def setUp(self):
        self.tables = json.loads((Path(__file__).resolve().parents[1]/'fixtures/dataset.json').read_text())['tables']

    def test_missing_is_not_zero_and_suspect_is_not_valid(self):
        tables = copy.deepcopy(self.tables)
        tables['readings'] = [dict(tables['readings'][0],quality='missing',rain_mm=None), dict(tables['readings'][0],reading_id='second',quality='suspect',rain_mm=80)]
        result = station_totals(tables)
        self.assertEqual(0, next(iter(result.values()))['valid_records'])
        self.assertEqual(2, next(iter(result.values()))['missing_records'])

    def test_position_reference_must_resolve(self):
        tables = copy.deepcopy(self.tables)
        tables['readings'][0]['position_version'] = 999
        with self.assertRaises(ValueError):
            station_totals(tables)

    def test_negative_valid_rain_is_rejected(self):
        tables = copy.deepcopy(self.tables)
        tables['readings'][0]['rain_mm'] = -1
        with self.assertRaises(ValueError):
            station_totals(tables)

    def test_moved_station_positions_remain_separate(self):
        result = station_totals(self.tables)
        self.assertIn('s0:1', result)
        self.assertIn('s0:2', result)
