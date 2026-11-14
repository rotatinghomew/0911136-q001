"""验证村域面积分析的空间支撑、口径分离与站网变动说明。"""

import copy
import datetime
import json
import unittest
from pathlib import Path

from rainnet.analysis import analyze_area, station_totals, summarize
from rainnet.observations import ObservationError
from rainnet.spatial import rectangle_area, thiessen_weights

ROOT = Path(__file__).resolve().parents[1]


def shifted(value, hours):
    return (datetime.datetime.fromisoformat(value) + datetime.timedelta(hours=hours)).isoformat()


class AreaAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads((ROOT / "fixtures/dataset.json").read_text())
        self.tables = self.data["tables"]

    def analysis(self, tables=None):
        return analyze_area(tables if tables is not None else self.tables)

    def test_point_summary_keeps_original_interface(self):
        summary = summarize(self.tables)
        self.assertEqual({"point_totals", "village_count", "scope"}, set(summary))
        self.assertEqual(station_totals(self.tables), summary["point_totals"])

    def test_processing_basis_is_relative_and_explained(self):
        analysis = self.analysis()
        basis = analysis["basis"]
        self.assertEqual(0.0, basis["start_offset_seconds"])
        self.assertEqual(3600.0, basis["duration_seconds"])
        self.assertIn("note", basis)
        self.assertTrue(basis["alternatives"])
        chosen = [item for item in basis["alternatives"]
                  if item["start_offset_seconds"] == basis["start_offset_seconds"]
                  and item["duration_seconds"] == basis["duration_seconds"]]
        self.assertEqual(1, len(chosen))

    def test_coverage_status_marks_single_station_and_absent_village(self):
        analysis = self.analysis()
        for identifier in ("v0", "v1", "v2"):
            village = analysis["villages"][identifier]
            self.assertEqual("estimated_single_station", village["status"])
            self.assertLess(village["window"]["coverage_min"], 1.0)
            self.assertGreaterEqual(village["window"]["coverage_min"], 0.5)
        absent = analysis["villages"]["v3"]
        self.assertEqual("not_inferable", absent["status"])
        self.assertIsNone(absent["window"])
        self.assertIsNone(absent["peak_area"])
        self.assertTrue(any("没有任何位置版本" in reason for reason in absent["reasons"]))
        self.assertEqual(0, absent["spatial_support"]["station_location_count"])

    def test_uncovered_area_gets_no_value(self):
        analysis = self.analysis()
        for village in analysis["villages"].values():
            for process in village["processes"]:
                covered = sum(location["area_m2"] for location in process["locations"]
                              if location["usable"])
                self.assertEqual(process["supported_area_m2"], round(covered, 3))
                self.assertLessEqual(process["coverage"], 1.0)
                if process["coverage_status"] in {"not_inferable", "below_coverage_threshold"}:
                    self.assertIsNone(process["covered_area_mean_mm"])
                    self.assertIsNone(process["interval95_mm"])

    def test_missing_station_does_not_become_zero_rain(self):
        tables = copy.deepcopy(self.tables)
        # s1 的原生区间正好是 1 小时口径；把它改成缺测。
        for reading in tables["readings"]:
            if reading["process_id"] == "p1" and reading["station_id"] == "s1":
                reading["quality"] = "missing"
                reading["rain_mm"] = None
        analysis = self.analysis(tables)
        process = next(item for item in analysis["villages"]["v0"]["processes"]
                       if item["process_id"] == "p1")
        self.assertEqual("not_inferable", process["coverage_status"])
        self.assertIsNone(process["covered_area_mean_mm"])
        self.assertIsNone(process["interval95_mm"])
        missing = [station for location in process["locations"]
                   for station in location["unusable_stations"]
                   if station["position"] == "s1:1"]
        self.assertEqual("missing_observation", missing[0]["reason"])
        village = analysis["villages"]["v0"]["window"]
        self.assertEqual(23, village["process_count_included"])
        self.assertEqual(24, village["process_count_total"])

    def test_unequal_observation_intervals_are_not_added_together(self):
        tables = copy.deepcopy(self.tables)
        for reading in tables["readings"]:
            if reading["station_id"] == "s0" and reading["process_id"] == "p0":
                reading["end"] = shifted(reading["start"], 1)
        analysis = self.analysis(tables)
        point = analysis["point_scope"]["positions"]["s0:1"]
        self.assertEqual([3600.0, 7200.0], point["valid_interval_durations_seconds"])
        self.assertFalse(point["intervals_consistent"])
        # 口径仍取 1 小时；s0 在 p0 的 1 小时段可用，但不能把 1 小时值当 2 小时用。
        selected = (analysis["basis"]["start_offset_seconds"],
                    analysis["basis"]["duration_seconds"])
        self.assertEqual((0.0, 3600.0), selected)

    def test_interval_that_overlaps_basis_start_is_rejected(self):
        tables = copy.deepcopy(self.tables)
        reading = tables["readings"][0]
        reading["start"] = shifted(reading["start"], -1)
        tables["stations"][0]["start"] = shifted(tables["stations"][0]["start"], -1)
        analysis = self.analysis(tables)
        # 观测跨过本过程起点：该站不能截取子区间，只能退出该口径。
        process = next(item for item in analysis["villages"]["v0"]["processes"]
                       if item["process_id"] == "p0")
        station = next(location for location in process["locations"]
                       if location["location"] == [0.0, 0.0])
        self.assertFalse(station["usable"])
        self.assertEqual("interval_mismatch", station["unusable_stations"][0]["reason"])

    def test_station_move_is_reported_not_blended(self):
        analysis = self.analysis()
        network = analysis["network"]
        moved = network["station_changes"]["moved_stations"]
        self.assertEqual(1, len(moved))
        self.assertEqual(["v0", "v1"], moved[0]["affected_villages"])
        self.assertIn("s0:1", network["station_changes"]["stopped_positions"][0]["position"])
        v0 = analysis["villages"]["v0"]
        v1 = analysis["villages"]["v1"]
        self.assertTrue(any(position["position"] == "s0:1" for position in v0["positions"]))
        self.assertTrue(any(position["position"] == "s0:2" for position in v1["positions"]))
        self.assertFalse(any(position["position"] == "s0:2" for position in v0["positions"]))

    def test_network_sensitivity_reports_stable_and_unstable_villages(self):
        analysis = self.analysis()
        network = analysis["network"]
        self.assertEqual(["s1:1", "s2:1", "s3:1", "s4:1", "s5:1"], network["stable_positions"])
        self.assertEqual("primary_not_estimated", network["comparison"]["v3"]["verdict"])
        self.assertIn(network["decision"]["verdict"],
                      {"station_network_robust", "station_network_changes_ranking",
                       "station_network_changes_conclusions"})
        v1_leave = {entry["station_id"]: entry["verdict"]
                    for entry in network["leave_one_station_out"]["v1"]}
        self.assertEqual("unestimable_without_station", v1_leave["s3"])

    def test_network_change_can_flip_a_village_conclusion(self):
        tables = copy.deepcopy(self.tables)
        # 统一到 2 小时口径，并把迁移后的 s0:2 雨量放大，使 v1 结论依赖变动站网。
        for reading in tables["readings"]:
            if reading["station_id"] != "s0":
                reading["end"] = shifted(reading["start"], 2)
            if reading["position_version"] == 2:
                reading["rain_mm"] = reading["rain_mm"] * 3
        analysis = self.analysis(tables)
        self.assertEqual(7200.0, analysis["basis"]["duration_seconds"])
        comparison = analysis["network"]["comparison"]["v1"]
        self.assertEqual("conclusion_changed", comparison["verdict"])
        self.assertGreater(comparison["relative_change"], 0.10)
        self.assertEqual("station_network_changes_conclusions",
                         analysis["network"]["decision"]["verdict"])
        self.assertIn("v1", analysis["network"]["decision"]["explanation"])

    def test_dependence_explains_why_stations_are_not_independent(self):
        analysis = self.analysis()
        dependence = analysis["dependence"]
        self.assertIn("不是独立样本", dependence["why_not_independent"])
        self.assertGreater(dependence["degrees_of_freedom"], 0)
        self.assertGreater(dependence["spatial_sd_mm"], 0)
        self.assertTrue(dependence["pairwise_correlations"])
        self.assertIn("length_scale_m", dependence["correlation_model"])
        for pair in dependence["pairwise_correlations"]:
            self.assertGreater(pair["distance_m"], 0)
            self.assertLessEqual(pair["correlation"], 1)

    def test_two_station_village_has_correlation_adjusted_design_effect(self):
        tables = copy.deepcopy(self.tables)
        # s3 原生 1 小时，与 s2 的 1 小时口径一致；v1 内出现两个可用位置。
        for reading in tables["readings"]:
            if reading["station_id"] == "s2":
                reading["end"] = shifted(reading["start"], 1)
        analysis = self.analysis(tables)
        self.assertEqual(3600.0, analysis["basis"]["duration_seconds"])
        village = analysis["villages"]["v1"]
        self.assertEqual("estimated_partial_coverage", village["status"])
        interval = next(process["interval"] for process in village["processes"]
                        if process["interval"] is not None)
        self.assertGreater(interval["design_effect"], 0)
        self.assertGreater(interval["design_effect"], interval["independent_design_effect"])
        self.assertIsNotNone(interval["mean_pairwise_correlation"])
        self.assertGreater(interval["mean_pairwise_correlation"], 0)

    def test_thiessen_weights_cover_rectangle_exactly(self):
        village = {"xmin": 0, "ymin": 0, "xmax": 600, "ymax": 800}
        points = [(0, 0), (300, 200)]
        weights = thiessen_weights(points, village)
        self.assertAlmostEqual(rectangle_area(village), sum(weights), places=6)
        self.assertEqual(2, len(weights))

    def test_identical_positions_are_combined_not_doubly_weighted(self):
        tables = copy.deepcopy(self.tables)
        for position in tables["stations"]:
            if position["station_id"] == "s1":
                position["x"], position["y"] = 0, 0
        analysis = self.analysis(tables)
        process = next(item for item in analysis["villages"]["v0"]["processes"]
                       if item["process_id"] == "p1")
        location = next(item for item in process["locations"]
                        if item["location"] == [0.0, 0.0])
        self.assertEqual(2, len(location["stations"]))
        self.assertAlmostEqual(1.0, sum(item["area_fraction"]
                                        for item in process["locations"]), places=6)
        # 同一点合并后只计一个位置权重，村域覆盖没有被重复计权。
        self.assertAlmostEqual(position_area_fraction(process), 1.0, places=6)

    def test_suspect_observations_are_not_used(self):
        tables = copy.deepcopy(self.tables)
        for reading in tables["readings"]:
            if reading["process_id"] == "p1" and reading["station_id"] == "s1":
                reading["quality"] = "suspect"
        analysis = self.analysis(tables)
        process = next(item for item in analysis["villages"]["v0"]["processes"]
                       if item["process_id"] == "p1")
        self.assertEqual("not_inferable", process["coverage_status"])
        self.assertIsNone(process["covered_area_mean_mm"])
        suspect = [station for location in process["locations"]
                   for station in location["unusable_stations"]
                   if station["position"] == "s1:1"]
        self.assertEqual("suspect_observation", suspect[0]["reason"])

    def test_no_readings_yields_explicit_not_inferable(self):
        tables = copy.deepcopy(self.tables)
        tables["readings"] = []
        analysis = self.analysis(tables)
        self.assertIsNone(analysis["basis"])
        for village in analysis["villages"].values():
            self.assertEqual("not_inferable", village["status"])
            self.assertIsNone(village["window"])

    def test_empty_village_has_no_distribution(self):
        tables = copy.deepcopy(self.tables)
        tables["villages"] = []
        analysis = self.analysis(tables)
        self.assertEqual({}, analysis["villages"])
        self.assertIsNone(analysis["regional"]["window_distribution"])

    def test_report_is_deterministic(self):
        first = self.analysis()
        second = self.analysis()
        self.assertEqual(json.dumps(first, sort_keys=True, ensure_ascii=False),
                         json.dumps(second, sort_keys=True, ensure_ascii=False))

    def test_every_included_process_has_interval_and_scope(self):
        analysis = self.analysis()
        for village in analysis["villages"].values():
            for process in village["processes"]:
                self.assertIn("interpretation", process)
                if process["covered_area_mean_mm"] is not None:
                    low, high = process["interval95_mm"]
                    self.assertLessEqual(low, process["covered_area_mean_mm"])
                    self.assertGreaterEqual(high, process["covered_area_mean_mm"])
                else:
                    self.assertIsNone(process["interval95_mm"])
            if village["window"] is not None:
                self.assertTrue(village["status_label"])
                self.assertTrue(village["window"]["process_count_included"] <=
                                village["window"]["process_count_total"])
                self.assertIsNotNone(village["window"]["interval"])

    def test_process_windows_do_not_borrow_readings_from_each_other(self):
        tables = {
            "stations": [
                {"station_id": "a", "position_version": 1, "x": 10, "y": 10,
                 "start": "2025-01-01T00:00:00+00:00", "end": None},
                {"station_id": "b", "position_version": 1, "x": 50, "y": 50,
                 "start": "2025-01-01T00:00:00+00:00", "end": None},
            ],
            "readings": [
                {"reading_id": "r1", "station_id": "a", "position_version": 1,
                 "start": "2025-06-01T12:00:00+00:00", "end": "2025-06-01T13:00:00+00:00",
                 "rain_mm": 5.0, "quality": "valid", "process_id": "p1"},
                {"reading_id": "r2", "station_id": "b", "position_version": 1,
                 "start": "2025-06-01T12:30:00+00:00", "end": "2025-06-01T13:30:00+00:00",
                 "rain_mm": 6.0, "quality": "valid", "process_id": "p2"},
                {"reading_id": "r3", "station_id": "a", "position_version": 1,
                 "start": "2025-06-01T14:00:00+00:00", "end": "2025-06-01T15:00:00+00:00",
                 "rain_mm": 4.0, "quality": "valid", "process_id": "p3"},
            ],
            "villages": [{"village_id": "v", "xmin": 0, "ymin": 0, "xmax": 100, "ymax": 100}],
        }
        analysis = analyze_area(tables)
        processes = {item["process_id"]: item
                     for item in analysis["villages"]["v"]["processes"]}
        # p1 与 p2 在绝对时间上部分重叠；每个过程只能看到自己的观测。
        self.assertEqual("estimated_single_station", processes["p2"]["coverage_status"])
        self.assertEqual(6.0, processes["p2"]["covered_area_mean_mm"])
        for identifier in ("p1", "p3"):
            self.assertIn(processes[identifier]["coverage_status"],
                          {"not_inferable", "below_coverage_threshold"})
            self.assertIsNone(processes[identifier]["covered_area_mean_mm"])
        self.assertEqual(["p1", "p3"],
                         sorted(analysis["villages"]["v"]["processes_not_included"]))
        self.assertEqual(1, analysis["villages"]["v"]["window"]["process_count_included"])

    def test_overlapping_process_windows_suppress_total(self):
        tables = {
            "stations": [
                {"station_id": "a", "position_version": 1, "x": 10, "y": 10,
                 "start": "2025-01-01T00:00:00+00:00", "end": None},
                {"station_id": "b", "position_version": 1, "x": 50, "y": 50,
                 "start": "2025-01-01T00:00:00+00:00", "end": None},
            ],
            "readings": [
                {"reading_id": "r1", "station_id": "a", "position_version": 1,
                 "start": "2025-06-01T12:00:00+00:00", "end": "2025-06-01T13:00:00+00:00",
                 "rain_mm": 5.0, "quality": "valid", "process_id": "p1"},
                {"reading_id": "r1b", "station_id": "b", "position_version": 1,
                 "start": "2025-06-01T12:00:00+00:00", "end": "2025-06-01T13:00:00+00:00",
                 "rain_mm": 7.0, "quality": "valid", "process_id": "p1"},
                {"reading_id": "r2", "station_id": "a", "position_version": 1,
                 "start": "2025-06-01T12:30:00+00:00", "end": "2025-06-01T13:30:00+00:00",
                 "rain_mm": 6.0, "quality": "valid", "process_id": "p2"},
                {"reading_id": "r2b", "station_id": "b", "position_version": 1,
                 "start": "2025-06-01T12:30:00+00:00", "end": "2025-06-01T13:30:00+00:00",
                 "rain_mm": 8.0, "quality": "valid", "process_id": "p2"},
            ],
            "villages": [{"village_id": "v", "xmin": 0, "ymin": 0, "xmax": 100, "ymax": 100}],
        }
        analysis = analyze_area(tables)
        window = analysis["villages"]["v"]["window"]
        self.assertTrue(window["process_count_included"] >= 1)
        if window["overlapping_process_windows"]:
            self.assertFalse(window["total_available"])
            self.assertIsNone(window["covered_area_total_mm"])

    def test_area_values_do_not_change_when_processes_are_repeated(self):
        """过程重复不改变空间支撑：同一口径下各村半宽一致，均值不因过程多而收窄。"""
        analysis = self.analysis()
        village = analysis["villages"]["v1"]
        halves = {process["interval"]["half_width_mm"] for process in village["processes"]
                  if process["interval"]}
        self.assertEqual(1, len(halves))
        half = halves.pop()
        window = village["window"]["interval"]
        mean = village["window"]["covered_area_mean_mm"]
        self.assertAlmostEqual(mean - half, window["mean_interval95_mm"][0], places=3)
        self.assertAlmostEqual(mean + half, window["mean_interval95_mm"][1], places=3)


def position_area_fraction(process):
    return sum(location["area_fraction"] for location in process["locations"]
               if location["usable"])
