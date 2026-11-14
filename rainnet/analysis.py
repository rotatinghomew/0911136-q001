"""在不均匀站网下，把按站摘要扩展成有空间支撑说明的村域降雨分析。

三个统计口径必须分开：

* **点观测**：某个位置版本在其原生观测区间内的累积雨量。不同位置的区间长度
  可能不同，不能直接比较、相加或换算成面雨量。
* **面积口径**：在同一"过程相对口径"下，只使用恰好铺满该口径的 valid 观测链，
  以村域内最近站分区的面积作权重，得到有站覆盖面积上的面积加权均值。未覆盖
  面积不作外推，缺测不当零雨，可疑观测不参与，区间不齐不截取也不折算。
* **推断范围**：只覆盖有可用观测的村域面积。没有村内站、观测不齐、只有单站或
  覆盖不足的村，结果必须显式标注，不能把单站值说成整村均匀降雨。

站与站之间、同一过程内的多个读数都不是独立样本：它们同属一次降雨过程，并在
空间上随距离相关。区间估计因此采用"空间方差 × 设计效应"的空间支撑近似，而
不是把站点数当样本量；过程数再多也不会把这种站网代表性误差按时间收窄。
"""

import math
import statistics
from collections import Counter, defaultdict
from datetime import timedelta

from .observations import ObservationIndex, iso
from .spatial import (max_corner_distance, nearest_neighbor_spacing,
                      point_in_rectangle, rectangle_area, thiessen_weights)

COVERAGE_THRESHOLD = 0.5
RELATIVE_CHANGE_TOLERANCE = 0.10
MAX_BASIS_CANDIDATES = 24
MAX_PAIRS_PER_PROCESS = 200
MIN_PAIR_PROCESSES = 5
T95 = ((1, 12.706), (2, 4.303), (3, 3.182), (4, 2.776), (5, 2.571),
       (6, 2.447), (7, 2.365), (8, 2.306), (9, 2.262), (10, 2.228),
       (12, 2.179), (15, 2.131), (20, 2.086), (24, 2.064), (30, 2.042),
       (40, 2.021), (60, 2.000), (120, 1.980))


def station_totals(tables):
    """现行按站位置累积摘要，保持原有输出与边界说明。"""
    positions = {(row["station_id"], row["position_version"]): row for row in tables["stations"]}
    result = defaultdict(lambda: {"rain_mm": 0.0, "valid_records": 0, "missing_records": 0})
    for row in tables["readings"]:
        key = (row["station_id"], row["position_version"])
        if key not in positions:
            raise ValueError("观测引用了不存在的位置版本")
        output = result[row["station_id"] + ":" + str(row["position_version"])]
        if row["quality"] == "valid" and row["rain_mm"] is not None:
            if row["rain_mm"] < 0:
                raise ValueError("雨量不能为负")
            output["rain_mm"] += row["rain_mm"]
            output["valid_records"] += 1
        else:
            output["missing_records"] += 1
    return dict(sorted(result.items()))


def summarize(tables):
    """现行按站摘要：保持原有字段与边界说明不变，另给出面积分析入口。"""
    return {"point_totals": station_totals(tables),
            "village_count": len(ObservationIndex(tables).villages),
            "scope": "按站位置累积摘要，不表示村域面积平均"}


def analyze_area(tables):
    """对一份同代三表资料生成村域面积口径分析。"""
    return area_analysis(ObservationIndex(tables))


def _round(value, digits=3):
    if value is None:
        return None
    return round(float(value), digits)


def _position_label(position):
    return position.station_id + ":" + str(position.position_version)


def _position_label_key(key):
    return key[0] + ":" + str(key[1])


def _t95(degrees):
    """双侧 95% t 分位数的保守近似：取不大于实际自由度的表值。"""
    if degrees < 1:
        return None
    for limit, value in T95:
        if degrees <= limit:
            return value
    return 1.960


def _distance_point_to_rectangle(x, y, rectangle):
    dx = max(rectangle["xmin"] - x, 0.0, x - rectangle["xmax"])
    dy = max(rectangle["ymin"] - y, 0.0, y - rectangle["ymax"])
    return math.hypot(dx, dy)


def _village_rectangle(village):
    return {"xmin": village.xmin, "ymin": village.ymin,
            "xmax": village.xmax, "ymax": village.ymax}


def _candidate_bases(index):
    """从各过程相对起点的观测边界收集候选累积口径。

    口径必须能跨站比较，所以用每个过程自己的起点偏移，而不是绝对时间。
    """
    counts = Counter()
    for process in index.processes:
        offsets = {process.base}
        for reading in process.readings:
            offsets.add(reading.start)
            offsets.add(reading.end)
        relative = sorted({offset - process.base for offset in offsets if offset > process.base})
        relative = [timedelta(0)] + relative
        pairs = [(start, end) for position, start in enumerate(relative)
                 for end in relative[position + 1:] if end > start]
        if len(pairs) > MAX_PAIRS_PER_PROCESS:
            pairs = sorted(pairs, key=lambda pair: (pair[1] - pair[0], pair[0]))[:MAX_PAIRS_PER_PROCESS]
        for pair in pairs:
            counts[pair] += 1
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0][1] - item[0][0], item[0][0]))
    return [basis for basis, _ in ranked[:MAX_BASIS_CANDIDATES]]


class _Analysis:
    def __init__(self, index):
        self.index = index
        self.value_cache = {}
        self.weight_cache = {}
        self.rectangles = {identifier: _village_rectangle(village)
                           for identifier, village in index.villages.items()}

    def position_value(self, key, start, end, process_id):
        cache_key = (key, start, end, process_id)
        if cache_key not in self.value_cache:
            self.value_cache[cache_key] = self.index.position_value(key, start, end, process_id)
        return self.value_cache[cache_key]

    def positions_in_village(self, village):
        rectangle = self.rectangles[village.village_id]
        return [position for position in self.index.positions.values()
                if point_in_rectangle(position.x, position.y, rectangle)]

    def village_for(self, position):
        for village in self.index.villages.values():
            if point_in_rectangle(position.x, position.y, self.rectangles[village.village_id]):
                return village.village_id
        return None

    def weights_for(self, village, coordinates):
        """同一组坐标在同一村的分区面积只算一次。"""
        cache_key = (village.village_id, tuple(coordinates))
        if cache_key not in self.weight_cache:
            self.weight_cache[cache_key] = thiessen_weights(
                coordinates, self.rectangles[village.village_id])
        return self.weight_cache[cache_key]

    def evaluate(self, village, process, basis, allowed=None):
        """一个村、一个过程、一个口径下的空间支撑与面积加权均值。"""
        start = process.base + basis[0]
        end = process.base + basis[1]
        considered = [position for position in self.positions_in_village(village)
                      if allowed is None or position.key in allowed]
        active = [position for position in considered if self.index.active(position, start, end)]
        inactive = [position for position in considered if position not in active]
        # 重合坐标先合并；它们的最近站分区相同，雨量取可用位置的平均。
        groups = defaultdict(list)
        for position in active:
            groups[(position.x, position.y)].append(position)
        coordinates = sorted(groups)
        weights = self.weights_for(village, coordinates)
        area = rectangle_area(self.rectangles[village.village_id])
        details = []
        usable_locations = []
        covered = 0.0
        for coordinate, weight in zip(coordinates, weights):
            positions = sorted(groups[coordinate], key=lambda item: item.key)
            values = []
            station_entries = []
            unusable = []
            for position in positions:
                value, reason = self.position_value(position.key, start, end, process.process_id)
                station_entries.append({"position": _position_label(position),
                                        "value_mm": _round(value, 3)})
                if value is None:
                    unusable.append({"position": _position_label(position), "reason": reason})
                else:
                    values.append(value)
            location_value = statistics.mean(values) if values else None
            if location_value is not None:
                covered += weight
                usable_locations.append((coordinate, location_value, weight))
            details.append({
                "location": [_round(coordinate[0], 3), _round(coordinate[1], 3)],
                "area_m2": _round(weight, 3),
                "area_fraction": _round(weight / area, 4),
                "usable": location_value is not None,
                "value_mm": _round(location_value, 3),
                "stations": station_entries,
                "unusable_stations": unusable,
            })
        coverage = covered / area if area > 0 else 0.0
        covered_mean = (sum(value * weight for _, value, weight in usable_locations) / covered
                        if covered > 0 else None)
        status = self._status(coverage, len(usable_locations))
        record = {
            "process_id": process.process_id,
            "start": iso(start),
            "end": iso(end),
            "coverage": _round(coverage, 4),
            "supported_area_m2": _round(covered, 3),
            "village_area_m2": _round(area, 3),
            "usable_location_count": len(usable_locations),
            "single_station_support": len(usable_locations) == 1,
            "coverage_status": status,
            "interpretation": _INTERPRETATION.get(status),
            "covered_area_mean_mm": _round(covered_mean, 3) if status in _ESTIMATED else None,
            "locations": details,
            "not_effective_positions": [_position_label(position) for position in inactive],
        }
        record["_usable"] = usable_locations
        record["_coverage_raw"] = coverage
        record["_covered_mean_raw"] = covered_mean
        record["_covered_raw"] = covered
        return record

    @staticmethod
    def _status(coverage, locations):
        if coverage <= 0:
            return "not_inferable"
        if coverage + 1e-9 < COVERAGE_THRESHOLD:
            return "below_coverage_threshold"
        if locations == 1:
            return "estimated_single_station"
        if coverage < 1.0 - 1e-9:
            return "estimated_partial_coverage"
        return "estimated"

    @staticmethod
    def public(record):
        return {key: value for key, value in record.items() if not key.startswith("_")}


_ESTIMATED = {"estimated", "estimated_single_station", "estimated_partial_coverage"}
_STATUS_LABELS = {
    "not_inferable": "无法推断：没有可用观测或覆盖面积为 0",
    "below_coverage_threshold": "无法推断：可用观测覆盖面积低于阈值",
    "estimated_single_station": "仅单站支撑：面积值等于该站所辖面积上的点值外推，村内差异未知",
    "estimated_partial_coverage": "部分面积支撑：只用有可用观测的覆盖面积，未覆盖部分不推断",
    "estimated": "完整面积支撑：村域被可用观测的位置分区覆盖",
    "estimated_partial_periods": "部分过程可估计：其余过程因覆盖或观测原因未纳入汇总",
}

_INTERPRETATION = {
    "estimated": "area_weighted_mean",
    "estimated_single_station": "point_value_extrapolated_to_station_catchment",
    "estimated_partial_coverage": "area_weighted_mean_on_covered_area",
    "below_coverage_threshold": "area_weighted_mean_below_threshold_not_reported",
    "not_inferable": "not_inferable",
}


def _basis_score(analysis, basis):
    """优先支撑面积；并列时取覆盖更多村-过程单元的口径。"""
    covered_total = 0.0
    usable_cells = 0
    for process in analysis.index.processes:
        for village in analysis.index.villages.values():
            record = analysis.evaluate(village, process, basis)
            covered_total += record["_coverage_raw"]
            if record["_usable"]:
                usable_cells += 1
    return covered_total, usable_cells


def _choose_basis(analysis):
    candidates = _candidate_bases(analysis.index)
    alternatives = []
    best = None
    best_key = None
    best_score = None
    for basis in candidates:
        covered, cells = _basis_score(analysis, basis)
        alternative = {"start_offset_seconds": basis[0].total_seconds(),
                       "duration_seconds": (basis[1] - basis[0]).total_seconds(),
                       "covered_fraction_sum": _round(covered, 6),
                       "usable_cells": cells}
        alternatives.append(alternative)
        key = (covered, cells, -(basis[1] - basis[0]).total_seconds(), -basis[0].total_seconds())
        if best_key is None or key > best_key:
            best_key, best, best_score = key, basis, (covered, cells)
    alternatives.sort(key=lambda item: (-item["covered_fraction_sum"], -item["usable_cells"],
                                        item["duration_seconds"], item["start_offset_seconds"]))
    basis_output = None
    if best is not None:
        basis_output = {
            "start_offset_seconds": best[0].total_seconds(),
            "duration_seconds": (best[1] - best[0]).total_seconds(),
            "rule": "选择可支撑村域覆盖面积最大的过程相对口径；面积并列时取时长短、起点早者",
            "score": {"covered_fraction_sum": _round(best_score[0], 6),
                      "usable_cells": best_score[1]},
            "alternatives": alternatives[:8],
            "note": ("该口径只用于跨站可比；与按站原生区间的点观测累积不能直接比较。"
                     "alternatives 展示若改用其他口径，可支撑的覆盖面积如何变化。"),
        }
    return best, basis_output, alternatives


def _spatial_statistics(analysis, basis):
    """同过程内站间差异估计空间方差；站间相关按距离拟合指数模型。"""
    per_process = []
    for process in analysis.index.processes:
        start = process.base + basis[0]
        end = process.base + basis[1]
        values = {}
        for position in analysis.index.positions.values():
            if not analysis.index.active(position, start, end):
                continue
            value, _ = analysis.position_value(position.key, start, end, process.process_id)
            if value is not None:
                values[position.key] = value
        per_process.append(values)
    sums = 0.0
    degrees = 0
    processes = 0
    for values in per_process:
        if len(values) >= 2:
            mean = statistics.mean(values.values())
            sums += sum((value - mean) ** 2 for value in values.values())
            degrees += len(values) - 1
            processes += 1
    variance = sums / degrees if degrees > 0 else None

    pairs = []
    keys = sorted({key for values in per_process for key in values})
    for position, first in enumerate(keys):
        for second in keys[position + 1:]:
            first_position = analysis.index.positions[first]
            second_position = analysis.index.positions[second]
            distance = math.hypot(first_position.x - second_position.x,
                                  first_position.y - second_position.y)
            if distance <= 0:
                continue
            common = [values for values in per_process if first in values and second in values]
            if len(common) < MIN_PAIR_PROCESSES:
                continue
            x = [values[first] for values in common]
            y = [values[second] for values in common]
            mean_x, mean_y = statistics.mean(x), statistics.mean(y)
            divisor = math.sqrt(sum((item - mean_x) ** 2 for item in x)
                                * sum((item - mean_y) ** 2 for item in y))
            if divisor <= 0:
                continue
            correlation = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y)) / divisor
            pairs.append({"station_a": _position_label(first_position),
                          "station_b": _position_label(second_position),
                          "distance_m": _round(distance, 3),
                          "correlation": _round(correlation, 4),
                          "common_processes": len(common)})

    model = None
    usable_pairs = [(item["distance_m"], item["correlation"]) for item in pairs
                    if item["correlation"] > 0 and item["distance_m"] > 0]
    if usable_pairs:
        numerator = sum(distance * math.log(correlation) for distance, correlation in usable_pairs)
        denominator = sum(distance * distance for distance, _ in usable_pairs)
        inverse = -numerator / denominator if denominator > 0 else 0.0
        if inverse > 0:
            model = {"form": "correlation = exp(-distance / length_scale_m)",
                     "length_scale_m": _round(1.0 / inverse, 3),
                     "pairs_used": len(usable_pairs),
                     "note": "由有限站对拟合，仅作近似；没有可用站对时按完全相关处理"}
    return {
        "method": "同过程内站间差异估计空间方差；站间相关按距离用指数模型拟合",
        "pooled_spatial_variance_mm2": _round(variance, 4) if variance is not None else None,
        "spatial_sd_mm": _round(math.sqrt(variance), 4) if variance is not None else None,
        "degrees_of_freedom": degrees,
        "processes_with_multiple_stations": processes,
        "pairwise_correlations": pairs,
        "correlation_model": model,
        "explanation": ("同一 process_id 的站间观测是同一次降雨过程的同步采样，不是独立样本；"
                        "站与站之间还随距离相关。因此不能用站数作样本量，也不能因为过程数多"
                        "就把区间按时间收窄。站点重合时先合并，避免把同一点重复计权。"),
    }


def _design_effect(weights, coordinates, model):
    """归一化面积权重的设计效应。

    返回 (design, independent, mean_correlation)：
    design 为按距离相关的加权方差相对单个独立位置的倍数；
    independent 为等权独立站位（不做相关修正）的对应倍数；
    mean_correlation 是站位之间按面积权重加权的平均相关系数。
    """
    normalized = [weight / sum(weights) for weight in weights]
    design = 0.0
    independent = sum(weight * weight for weight in normalized)
    correlation_weight = 0.0
    pair_weight = 0.0
    for first, (coordinate_a, weight_a) in enumerate(zip(coordinates, normalized)):
        for second, (coordinate_b, weight_b) in enumerate(zip(coordinates, normalized)):
            if first == second:
                correlation = 1.0
            elif model is None:
                correlation = 1.0
            else:
                distance = math.hypot(coordinate_a[0] - coordinate_b[0],
                                      coordinate_a[1] - coordinate_b[1])
                correlation = math.exp(-distance / model["length_scale_m"])
            design += weight_a * weight_b * correlation
            if first < second:
                correlation_weight += weight_a * weight_b * correlation
                pair_weight += weight_a * weight_b
    mean_correlation = correlation_weight / pair_weight if pair_weight > 0 else None
    return design, independent, mean_correlation


def _interval(record, coordinates, weights, spatial, model):
    if spatial["pooled_spatial_variance_mm2"] is None or not weights:
        return None
    design, independent, mean_correlation = _design_effect(weights, coordinates, model)
    standard_error = math.sqrt(spatial["pooled_spatial_variance_mm2"] * design)
    quantile = _t95(spatial["degrees_of_freedom"])
    if quantile is None:
        return None
    half = quantile * standard_error
    return {
        "method": "空间方差 × 设计效应，t 区间",
        "standard_error_mm": _round(standard_error, 3),
        "half_width_mm": _round(half, 3),
        "interval95_mm": [_round(max(0.0, record["_covered_mean_raw"] - half), 3),
                          _round(record["_covered_mean_raw"] + half, 3)],
        "design_effect": _round(design, 4),
        "independent_design_effect": _round(independent, 4),
        "mean_pairwise_correlation": _round(mean_correlation, 4),
        "effective_station_count": _round(1.0 / design, 3) if design > 0 else None,
        "variance_degrees_of_freedom": spatial["degrees_of_freedom"],
        "basis": ("空间代表性近似区间：假设该村覆盖面积的次站网空间方差与全网相近。"
                  "同一站网的偏差持续存在，不随过程数收窄；单站位置无法给出村内代表性区间。"),
    }


def _add_interval(record, spatial, model):
    if record["covered_area_mean_mm"] is None:
        record["interval95_mm"] = None
        record["interval"] = None
        return
    coordinates = [coordinate for coordinate, _, _ in record["_usable"]]
    weights = [weight for _, _, weight in record["_usable"]]
    if record["usable_location_count"] == 1:
        weights = [1.0]
    interval = _interval(record, coordinates, weights, spatial, model)
    if interval is not None and record["usable_location_count"] == 1:
        interval["single_location_note"] = ("只有一个可用位置，设计效应为 1；区间只反映全网空间方差，"
                                            "不包含村内空间代表性，不能当作村域真值区间")
    record["interval95_mm"] = interval["interval95_mm"] if interval else None
    record["interval"] = interval


def _village_series(analysis, village, basis, spatial, model, allowed=None):
    records = []
    for process in analysis.index.processes:
        record = analysis.evaluate(village, process, basis, allowed)
        _add_interval(record, spatial, model)
        records.append(record)
    included = [record for record in records if record["coverage_status"] in _ESTIMATED]
    summary = {"processes": [analysis.public(record) for record in records],
               "processes_included": [record["process_id"] for record in included],
               "processes_not_included": [record["process_id"] for record in records
                                          if record not in included]}
    if not included:
        summary["window"] = None
        summary["peak_area"] = None
        summary["coverage_status"] = "not_inferable"
    else:
        means = [record["covered_area_mean_mm"] for record in included]
        mean = statistics.mean(means)
        # 过程窗口若在绝对时间上重叠，把它们相加会重复计雨；只给均值并明确标注。
        windows = sorted((record["start"], record["end"], record["process_id"])
                          for record in included)
        overlapping = [[previous[2], current[2]]
                       for previous, current in zip(windows, windows[1:])
                       if current[0] < previous[1]]
        total = None if overlapping else sum(means)
        half_widths = [record["interval"]["half_width_mm"] for record in included
                       if record["interval"] and record["interval"].get("half_width_mm") is not None]
        representative_half = max(half_widths) if half_widths else None
        window_interval = None
        if representative_half is not None:
            window_interval = {
                "mean_interval95_mm": [_round(max(0.0, mean - representative_half), 3),
                                       _round(mean + representative_half, 3)],
                "total_interval95_mm": (None if total is None else
                                        [_round(max(0.0, total - representative_half * len(included)), 3),
                                         _round(total + representative_half * len(included), 3)]),
                "half_width_source": ("各过程中最大空间半宽；过程均值不因平均而收窄，"
                                      "总量按过程数线性放大，假设站网偏差持续存在"),
                "note": "空间支撑近似区间，不是时间抽样置信区间",
            }
        peak = max(included, key=lambda record: (record["covered_area_mean_mm"], record["process_id"]))
        summary["window"] = {
            "covered_area_mean_mm": _round(mean, 3),
            "covered_area_total_mm": _round(total, 3) if total is not None else None,
            "total_available": total is not None,
            "overlapping_process_windows": overlapping,
            "process_count_included": len(included),
            "process_count_total": len(records),
            "coverage_min": _round(min(record["coverage"] for record in included), 4),
            "coverage_max": _round(max(record["coverage"] for record in included), 4),
            "single_station_processes": sum(1 for record in included
                                            if record["coverage_status"] == "estimated_single_station"),
            "interval": window_interval,
            "complete_period": len(included) == len(records),
        }
        summary["peak_area"] = {
            "process_id": peak["process_id"],
            "covered_area_mean_mm": peak["covered_area_mean_mm"],
            "coverage": peak["coverage"],
            "interval95_mm": peak["interval95_mm"],
            "scope": "该村有站覆盖面积上的过程最大值，不是村内任意点或全村干部的极值",
        }
        summary["coverage_status"] = _village_status(included, records)
    summary["reasons"] = _village_reasons(analysis, village, records)
    return summary, records


def _village_status(included, records):
    if len(included) < len(records):
        return "estimated_partial_periods"
    if any(record["coverage_status"] == "estimated_single_station" for record in included):
        return "estimated_single_station"
    if any(record["coverage_status"] == "estimated_partial_coverage" for record in included):
        return "estimated_partial_coverage"
    return "estimated"


def _village_reasons(analysis, village, records):
    positions = analysis.positions_in_village(village)
    reasons = []
    if not positions:
        reasons.append("村域内没有任何位置版本；不能给出面积推断，未使用村外站外推")
        nearest = _nearest_outside(analysis, village)
        if nearest:
            reasons.append("最近站 " + nearest["position"] + " 在村外 "
                           + str(nearest["distance_m"]) + " 米，仅作距离参照")
        return reasons
    if not any(record["coverage_status"] in _ESTIMATED for record in records):
        reason_counts = Counter()
        for record in records:
            for location in record["locations"]:
                for entry in location["unusable_stations"]:
                    reason_counts[entry["reason"]] += 1
        if reason_counts:
            reasons.append("没有可用观测：" + "、".join(sorted(reason_counts)))
    if any(record["coverage_status"] == "below_coverage_threshold" for record in records):
        reasons.append("部分过程覆盖面积低于阈值 " + str(COVERAGE_THRESHOLD) + "，未纳入窗口汇总")
    if any(record["coverage_status"] == "estimated_single_station" for record in records):
        reasons.append("部分过程只有单站可用；面积值实为该站所辖面积上的点值外推，村内空间变化未知")
    if any(record["coverage_status"] == "estimated_partial_coverage" for record in records):
        reasons.append("部分村域面积没有可用观测；估计只覆盖有站面积")
    if positions and any(record["not_effective_positions"] for record in records):
        reasons.append("部分位置版本在过程口径内未生效或因迁移/停测被排除")
    return reasons


def _nearest_outside(analysis, village):
    candidates = []
    rectangle = analysis.rectangles[village.village_id]
    for position in analysis.index.positions.values():
        if point_in_rectangle(position.x, position.y, rectangle):
            continue
        distance = _distance_point_to_rectangle(position.x, position.y, rectangle)
        candidates.append((distance, position.station_id, position.position_version, position))
    if not candidates:
        return None
    distance, _, _, position = min(candidates)
    return {"position": _position_label(position), "distance_m": _round(distance, 3),
            "scope": "仅作位置参照，未作为村域面积观测"}


def _station_inventory(analysis, village):
    entries = [{"position": _position_label(position), "x": position.x, "y": position.y,
                "effective_start": iso(position.start), "effective_end": iso(position.end)}
               for position in analysis.positions_in_village(village)]
    return sorted(entries, key=lambda item: (item["effective_start"], item["position"]))


def _relative_change(primary, alternative):
    if primary is None or alternative is None:
        return None
    return abs(alternative - primary) / max(abs(primary), 0.5)


def _stable_positions(analysis, basis):
    return {position.key for position in analysis.index.positions.values()
            if all(analysis.index.active(position, process.base + basis[0],
                                         process.base + basis[1])
                   for process in analysis.index.processes)}


def _station_changes(analysis, window_start, window_end, basis):
    versions = defaultdict(list)
    for position in analysis.index.positions.values():
        versions[position.station_id].append(position)
    moved = []
    for station_id, entries in sorted(versions.items()):
        if len(entries) < 2:
            continue
        entries.sort(key=lambda position: position.start)
        assigned = [{"position": _position_label(position), "x": position.x, "y": position.y,
                     "village_id": analysis.village_for(position),
                     "effective_start": iso(position.start), "effective_end": iso(position.end)}
                    for position in entries]
        moved.append({"station_id": station_id, "versions": assigned,
                      "analysis_window_overlap": any(
                          position.start < window_end and (position.end is None or position.end > window_start)
                          for position in entries),
                      "affected_villages": sorted({entry["village_id"] for entry in assigned
                                                   if entry["village_id"] is not None})})
    stopped = []
    for position in analysis.index.positions.values():
        if position.end is None:
            continue
        stopped.append({"position": _position_label(position),
                        "effective_end": iso(position.end),
                        "village_id": analysis.village_for(position),
                        "within_analysis_window": position.end <= window_end and position.end >= window_start,
                        "replacement": any(other.station_id == position.station_id
                                           and other.start >= position.end
                                           for other in versions[position.station_id])})
    late = [{"position": _position_label(position), "effective_start": iso(position.start),
             "village_id": analysis.village_for(position)}
            for position in analysis.index.positions.values() if position.start > window_start]
    return {
        "moved_stations": moved,
        "stopped_positions": sorted(stopped, key=lambda item: (item["effective_end"], item["position"])),
        "positions_starting_after_window": sorted(late, key=lambda item: (item["effective_start"], item["position"])),
        "window_start": iso(window_start),
        "window_end": iso(window_end),
        "basis_duration_seconds": (basis[1] - basis[0]).total_seconds(),
        "method": ("位置版本按有效期分开处理；迁移前后不合并成同一序列，也不跨版本补值。"
                   "只有位置有效期完整包住口径区间时才参与面积权重；"
                   "停测和晚启用的位置不参加该过程，也不在窗口内补造观测。"),
    }


def _network_sensitivity(analysis, basis, spatial, model):
    stable = _stable_positions(analysis, basis)
    window_start = min((process.base for process in analysis.index.processes), default=None)
    window_end = max((process.base + basis[1] for process in analysis.index.processes), default=None)
    comparison = {}
    rankings_primary = {}
    rankings_stable = {}
    for village in analysis.index.villages.values():
        primary_summary, _ = _village_series(analysis, village, basis, spatial, model)
        stable_summary, _ = _village_series(analysis, village, basis, spatial, model, stable)
        primary = primary_summary["window"]["covered_area_mean_mm"] if primary_summary["window"] else None
        alternative = stable_summary["window"]["covered_area_mean_mm"] if stable_summary["window"] else None
        change = _relative_change(primary, alternative)
        if primary is None:
            verdict = "primary_not_estimated"
        elif alternative is None:
            verdict = "depends_on_changed_network"
        elif change is not None and change > RELATIVE_CHANGE_TOLERANCE:
            verdict = "conclusion_changed"
        else:
            verdict = "conclusion_stable"
        comparison[village.village_id] = {
            "primary_covered_area_mean_mm": primary,
            "stable_network_covered_area_mean_mm": alternative,
            "relative_change": _round(change, 4),
            "verdict": verdict,
            "change_tolerance": RELATIVE_CHANGE_TOLERANCE,
            "primary_coverage_min": primary_summary["window"]["coverage_min"] if primary_summary["window"] else None,
            "stable_coverage_min": stable_summary["window"]["coverage_min"] if stable_summary["window"] else None,
        }
        if primary is not None:
            rankings_primary[village.village_id] = primary
        if alternative is not None:
            rankings_stable[village.village_id] = alternative
    shared = sorted(set(rankings_primary) & set(rankings_stable))
    shared_primary = sorted(shared, key=lambda identifier: (-rankings_primary[identifier], identifier))
    shared_stable = sorted(shared, key=lambda identifier: (-rankings_stable[identifier], identifier))
    leave_one_out = {}
    for village in analysis.index.villages.values():
        base_summary, _ = _village_series(analysis, village, basis, spatial, model)
        base = base_summary["window"]["covered_area_mean_mm"] if base_summary["window"] else None
        entries = []
        station_ids = sorted({position.station_id for position in analysis.positions_in_village(village)})
        for station_id in station_ids:
            allowed = {position.key for position in analysis.index.positions.values()
                       if position.station_id != station_id}
            alternative_summary, _ = _village_series(analysis, village, basis, spatial, model, allowed)
            alternative = alternative_summary["window"]["covered_area_mean_mm"] if alternative_summary["window"] else None
            change = _relative_change(base, alternative)
            entries.append({"station_id": station_id,
                            "covered_area_mean_without_station_mm": alternative,
                            "relative_change": _round(change, 4),
                            "verdict": ("unestimable_without_station" if alternative is None else
                                        "over_tolerance" if change is not None and change > RELATIVE_CHANGE_TOLERANCE
                                        else "within_tolerance")})
        leave_one_out[village.village_id] = entries
    changes = _station_changes(analysis, window_start, window_end, basis)
    return {
        "station_changes": changes,
        "stable_positions": sorted(_position_label_key(key) for key in stable),
        "comparison": comparison,
        "ranking": {"villages_compared": shared, "primary_order": shared_primary,
                    "stable_order": shared_stable, "stable": shared_primary == shared_stable,
                    "variance_rule": "排名并列时按村编号排序",
                    "note": "仅在两种站网下都能估计的村之间比较；排序变化说明站网变动影响相对研判"},
        "leave_one_station_out": leave_one_out,
        "decision": _network_decision(comparison, shared_primary, shared_stable),
        "note": ("稳定站网指在全部过程口径内都有效的位置版本；迁移、停测、晚启用的位置被排除。"
                 "结论变化按相对差 " + str(RELATIVE_CHANGE_TOLERANCE) + " 判定。"),
    }


def _network_decision(comparison, primary_order, stable_order):
    changed = [identifier for identifier, entry in comparison.items()
               if entry["verdict"] in {"conclusion_changed", "depends_on_changed_network"}]
    if changed:
        verdict = "station_network_changes_conclusions"
        explanation = ("排除迁移、停测和晚启用位置后，以下村的结论发生变化或不能再估计："
                       + "、".join(sorted(changed)) + "。对外使用前应说明站网依赖。")
    elif primary_order != stable_order:
        verdict = "station_network_changes_ranking"
        explanation = "各村数值在容差内，但相对排序变化；不能把排序当成稳健结论。"
    else:
        verdict = "station_network_robust"
        explanation = "在稳定站网子集下，可估计村的数值在容差内且排序不变。"
    return {"verdict": verdict, "explanation": explanation,
            "policy": ("失效率和停测位置不参与稳定站网比较；不稳定的结论不得包装成"
                       "单纯站网变化前后的同比。")}


def _point_extreme(index):
    candidates = [reading for reading in index.readings
                  if reading.quality == "valid" and reading.rain_mm is not None]
    if not candidates:
        return None
    peak = max(candidates, key=lambda reading: (reading.rain_mm, reading.reading_id))
    return {"reading_id": peak.reading_id, "position": _position_label(peak),
            "process_id": peak.process_id, "start": iso(peak.start), "end": iso(peak.end),
            "rain_mm": _round(peak.rain_mm, 3),
            "scope": ("单站原生区间点观测极值；区间长度可能不同，"
                      "不能当作村域面积极值或瞬时强度")}


def _regional(index, villages, basis, spatial, model):
    point_extreme = _point_extreme(index)
    if basis is None:
        return {"per_process": [], "window_distribution": None, "area_extreme": None,
                "point_extreme": point_extreme,
                "note": "没有可用的过程口径，不能形成区域分布"}
    per_process = []
    for process in index.processes:
        included = []
        excluded = []
        for village in index.villages.values():
            record = villages[village.village_id]["processes_by_id"][process.process_id]
            if record["covered_area_mean_mm"] is not None:
                included.append({"village_id": village.village_id,
                                 "covered_area_mean_mm": record["covered_area_mean_mm"],
                                 "coverage": record["coverage"],
                                 "single_station": record["coverage_status"] == "estimated_single_station",
                                 "interval95_mm": record["interval95_mm"]})
            else:
                excluded.append({"village_id": village.village_id,
                                 "coverage_status": record["coverage_status"]})
        values = [entry["covered_area_mean_mm"] for entry in included]
        per_process.append({
            "process_id": process.process_id,
            "start": iso(process.base + basis[0]),
            "end": iso(process.base + basis[1]),
            "villages_included": included,
            "villages_not_included": excluded,
            "min_mm": _round(min(values), 3) if values else None,
            "median_mm": _round(statistics.median(values), 3) if values else None,
            "max_mm": _round(max(values), 3) if values else None,
            "range_mm": _round(max(values) - min(values), 3) if values else None,
            "range_scope": ("只在有可用面积估计的村之间计算；不是全县真值范围，"
                            "也未计各区间宽度重叠"),
        })
    window_values = []
    window_villages = []
    for village in index.villages.values():
        window = villages[village.village_id]["window"]
        if window is not None:
            window_values.append(window["covered_area_mean_mm"])
            window_villages.append(village.village_id)
    window_distribution = None
    if window_values:
        widest_halves = [villages[village.village_id]["window"]["interval"]["mean_interval95_mm"]
                         for village in index.villages.values()
                         if villages[village.village_id]["window"]
                         and villages[village.village_id]["window"]["interval"]]
        widest = None
        if widest_halves:
            widest = max((item for item in widest_halves if item),
                         key=lambda item: item[1] - item[0], default=None)
        window_distribution = {
            "villages_included": window_villages,
            "villages_not_included": [village.village_id for village in index.villages.values()
                                      if village.village_id not in window_villages],
            "min_mm": _round(min(window_values), 3),
            "median_mm": _round(statistics.median(window_values), 3),
            "max_mm": _round(max(window_values), 3),
            "range_mm": _round(max(window_values) - min(window_values), 3),
            "single_station_villages": [village.village_id for village in index.villages.values()
                                        if villages[village.village_id]["coverage_status"]
                                        == "estimated_single_station"],
            "widest_mean_interval95_mm": widest,
            "note": ("只在有可用面积估计的村之间汇总；未推断的村不参与，因此不是全县真值分布。"
                     "单站村的值是点值外推，村内差异未知；区间重叠时不能仅凭点值排名断定谁更多。"),
        }
    peaks = [(villages[village.village_id]["peak_area"], village.village_id)
             for village in index.villages.values()
             if villages[village.village_id]["peak_area"] is not None]
    area_extreme = None
    if peaks:
        peak, village_id = max(peaks, key=lambda item: (item[0]["covered_area_mean_mm"], item[1]))
        area_extreme = {"village_id": village_id, "process_id": peak["process_id"],
                        "covered_area_mean_mm": peak["covered_area_mean_mm"],
                        "coverage": peak["coverage"], "interval95_mm": peak["interval95_mm"],
                        "scope": "有站覆盖面积上的村域过程均值极值；不是村内任意点极值"}
    return {"per_process": per_process, "window_distribution": window_distribution,
            "area_extreme": area_extreme, "point_extreme": point_extreme,
            "note": ("区域分布只在可估计的村域覆盖面积上形成；未覆盖面积和未推断村不在其中。"
                     "局地点极值与面积极值是不同口径。")}


def _basis_sensitivity(analysis, spatial, model, alternatives):
    rows = []
    for alternative in alternatives:
        basis = (timedelta(seconds=alternative["start_offset_seconds"]),
                 timedelta(seconds=alternative["start_offset_seconds"] + alternative["duration_seconds"]))
        statuses = Counter()
        means = []
        for village in analysis.index.villages.values():
            summary, _ = _village_series(analysis, village, basis, spatial, model)
            if summary["window"] is None:
                statuses["not_inferable"] += 1
            else:
                statuses["estimated"] += 1
                means.append(summary["window"]["covered_area_mean_mm"])
        rows.append({"start_offset_seconds": alternative["start_offset_seconds"],
                     "duration_seconds": alternative["duration_seconds"],
                     "villages_estimated": statuses["estimated"],
                     "villages_not_inferable": statuses["not_inferable"],
                     "window_range_mm": ([_round(min(means), 3), _round(max(means), 3)]
                                         if means else None),
                     "covered_fraction_sum": alternative["covered_fraction_sum"]})
    return {
        "alternatives": rows,
        "note": ("改变跨站累积口径会改变哪些站可用、哪些村可推断，从而可能改变村域结论。"
                 "表中只用同一代资料，未混用其他代际。"),
    }


def _spatial_support(analysis, village):
    positions = analysis.positions_in_village(village)
    coordinates = [(position.x, position.y) for position in positions]
    spacing = nearest_neighbor_spacing(coordinates)
    corner = max_corner_distance(coordinates, analysis.rectangles[village.village_id])
    return {
        "station_location_count": len(positions),
        "village_area_m2": _round(rectangle_area(analysis.rectangles[village.village_id]), 3),
        "typical_nearest_neighbor_spacing_m": _round(spacing, 3),
        "max_corner_distance_to_nearest_location_m": _round(corner, 3),
        "resident_station_count": sum(1 for position in positions if position.end is None),
        "meaning": ("典型间距越大、角落离最近站越远，村内未观测部分的空间不确定性越大。"
                    "这些是几何距离，不是降雨相关范围；不能据此把缺少观测的角落说成有推断。"),
    }


def _village_output(analysis, village, summary):
    output = {
        "status": summary["coverage_status"],
        "status_label": _STATUS_LABELS.get(summary["coverage_status"]),
        "window": summary["window"],
        "peak_area": summary["peak_area"],
        "spatial_support": _spatial_support(analysis, village),
        "reasons": summary["reasons"],
        "positions": _station_inventory(analysis, village),
        "processes": summary["processes"],
        "_by_id": {record["process_id"]: record for record in summary["processes"]},
    }
    if summary["processes_not_included"]:
        output["processes_not_included"] = summary["processes_not_included"]
    return output


def _point_scope(index):
    """点观测口径的单独说明：原生区间可能不等长，不能与面积口径直接比较。"""
    positions = {}
    for position in index.positions.values():
        readings = index.readings_by_position.get(position.key, ())
        valid = [reading for reading in readings
                 if reading.quality == "valid" and reading.rain_mm is not None]
        durations = {(reading.end - reading.start).total_seconds() for reading in valid}
        ordered = sorted(valid, key=lambda reading: reading.start)
        contiguous = bool(ordered) and all(
            previous.end == current.start for previous, current in zip(ordered, ordered[1:]))
        total = sum(reading.rain_mm for reading in valid)
        positions[_position_label(position)] = {
            "village_id": None,
            "native_total_mm": _round(total, 3),
            "valid_records": len(valid),
            "missing_records": sum(1 for reading in readings if reading.quality == "missing"),
            "suspect_records": sum(1 for reading in readings if reading.quality == "suspect"),
            "valid_interval_durations_seconds": sorted(durations),
            "intervals_consistent": len(durations) <= 1,
            "valid_intervals_contiguous": contiguous,
            "valid_span": ([iso(ordered[0].start), iso(ordered[-1].end)] if ordered else None),
            "first_observation": iso(min((reading.start for reading in readings), default=None)) if readings else None,
            "last_observation": iso(max((reading.end for reading in readings), default=None)) if readings else None,
        }
    for village in index.villages.values():
        for position in index.positions.values():
            if point_in_rectangle(position.x, position.y, _village_rectangle(village)):
                positions[_position_label(position)]["village_id"] = village.village_id
    return {
        "positions": positions,
        "scope": ("点观测按各位置原生区间累积，保留现有按站摘要口径。"
                  "intervals_consistent 为 false 表示同一站点的有效观测区间不等长；"
                  "valid_intervals_contiguous 为 false 表示有效观测之间存在缺口。"
                  "这两种情况都不能把 native_total_mm 直接与别站比较或换算成面雨量。"
                  "missing_records 与 suspect_records 不计入 native_total_mm，缺测不是零雨。"),
    }


def _limitations():
    return [
        "面积估计只覆盖有可用观测的村域面积；coverage 小于 1 时，未覆盖部分没有推断值。",
        "只有单站可用的村，面积值等于该站所辖面积上的点值外推；村内空间差异未知，不能当作均匀降雨。",
        "点观测使用各自原生区间累积，区间长度不同，不能直接相加、比较或换算成面积口径。",
        "缺测、可疑、观测区间不齐、停测和位置迁移都不折算成零雨，也不插值、不跨版本补值。",
        "同一 process_id 的站间观测属于同一次降雨过程，并随距离相关；站点数不是独立样本量。",
        "区间估计是空间支撑近似区间；同一站网的偏差持续存在，增加过程数不会按时间收窄。",
        "村内没有位置版本时不做外推；最近站只作为距离参照。",
    ]


def area_analysis(index):
    analysis = _Analysis(index)
    basis, basis_output, alternatives = _choose_basis(analysis)
    villages = {}
    if basis is None:
        for village in index.villages.values():
            villages[village.village_id] = {
                "status": "not_inferable", "status_label": _STATUS_LABELS["not_inferable"],
                "window": None, "peak_area": None,
                "reasons": ["没有观测记录，无法定义共同累积口径"],
                "positions": _station_inventory(analysis, village), "processes": [], "_by_id": {},
            }
        return {"scope": "村域面积口径分析", "basis": None, "villages": villages,
                "point_scope": _point_scope(index),
                "regional": _regional(index, {"per_process": []}, None, None, None),
                "dependence": None, "network": None, "basis_sensitivity": None,
                "limitations": _limitations(),
                "thresholds": {"coverage": COVERAGE_THRESHOLD,
                               "relative_change": RELATIVE_CHANGE_TOLERANCE}}
    spatial = _spatial_statistics(analysis, basis)
    model = spatial["correlation_model"]
    for village in index.villages.values():
        summary, _ = _village_series(analysis, village, basis, spatial, model)
        villages[village.village_id] = _village_output(analysis, village, summary)
    regional_villages = {identifier: {"processes_by_id": value.pop("_by_id"),
                                      "window": value["window"],
                                      "peak_area": value["peak_area"],
                                      "coverage_status": value["status"]}
                         for identifier, value in villages.items()}
    regional = _regional(index, regional_villages, basis, spatial, model)
    network = _network_sensitivity(analysis, basis, spatial, model)
    basis_sensitivity = _basis_sensitivity(analysis, spatial, model, alternatives)
    design_examples = {}
    for village in index.villages.values():
        candidates = [record for record in villages[village.village_id]["processes"]
                      if record["interval"] is not None]
        representative = max(candidates,
                             key=lambda record: (record["usable_location_count"], record["coverage"]),
                             default=None)
        design_examples[village.village_id] = {
            "process_id": representative["process_id"] if representative else None,
            "nominal_location_count": representative["usable_location_count"] if representative else 0,
            "design_effect_for_that_process": representative["interval"]["design_effect"] if representative else None,
            "effective_location_count_for_that_process": representative["interval"]["effective_station_count"] if representative else None,
            "correlation_inflation_vs_independent": (
                _round(representative["interval"]["design_effect"]
                       / representative["interval"]["independent_design_effect"], 3)
                if representative and representative["interval"]["independent_design_effect"] else None),
            "note": ("design_effect > 1 表示按距离相关的站位等效独立样本少于名义站数；"
                     "correlation_inflation_vs_independent 是相对等权独立平均的方差放大倍数。"
                     "取该村可用位置最多的过程作示例，实际逐过程计算。"),
        }
    dependence = dict(spatial)
    dependence["design_effect_by_village"] = design_examples
    dependence["why_not_independent"] = spatial["explanation"]
    return {
        "scope": ("村域面积口径：只对村域内有可用观测的覆盖面积形成估计；"
                  "未覆盖面积不作外推。点观测与面积口径分开呈现。"),
        "basis": basis_output,
        "point_scope": _point_scope(index),
        "villages": villages,
        "regional": regional,
        "dependence": dependence,
        "network": network,
        "basis_sensitivity": basis_sensitivity,
        "limitations": _limitations(),
        "thresholds": {"coverage": COVERAGE_THRESHOLD,
                       "relative_change": RELATIVE_CHANGE_TOLERANCE},
    }
