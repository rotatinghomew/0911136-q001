"""同代观测资料的索引、校验与累积区间取值。

本模块只处理一个 generation 内的完整记录：位置版本、观测区间与村域边界必须
来自同一次供数快照，跨代混用应在读取快照时即被拒绝。

取值规则刻意保守：只有一条从区间起点连续铺到终点、且每段都是 valid 的观测链
才会得到一个累积值。区间不齐、重叠、换测或停测都不会被折算、插值或当作 0；
缺测的记录不会贡献任何雨量，也不会被当作零雨。
"""

import math
from dataclasses import dataclass
from datetime import datetime, timezone


class ObservationError(ValueError):
    """资料不满足领域契约。"""


QUALITIES = {"valid", "suspect", "missing"}


def timestamp(value, field):
    if not isinstance(value, str) or not value:
        raise ObservationError(f"{field} 必须是 ISO 时间字符串")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise ObservationError(f"{field} 不是有效 ISO 时间：{value!r}") from error
    if parsed.tzinfo is None:
        raise ObservationError(f"{field} 必须带时区：{value!r}")
    return parsed.astimezone(timezone.utc)


def iso(value):
    return None if value is None else value.isoformat()


def number(value, field):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ObservationError(f"{field} 必须是有限数值")
    return float(value)


@dataclass(frozen=True)
class Position:
    station_id: str
    position_version: int
    x: float
    y: float
    start: datetime
    end: datetime | None

    @property
    def key(self):
        return (self.station_id, self.position_version)


@dataclass(frozen=True)
class Reading:
    reading_id: str
    station_id: str
    position_version: int
    start: datetime
    end: datetime
    rain_mm: float | None
    quality: str
    process_id: str

    @property
    def position_key(self):
        return (self.station_id, self.position_version)


@dataclass(frozen=True)
class Process:
    process_id: str
    base: datetime
    readings: tuple


@dataclass(frozen=True)
class Village:
    village_id: str
    xmin: float
    ymin: float
    xmax: float
    ymax: float

    @property
    def area(self):
        return (self.xmax - self.xmin) * (self.ymax - self.ymin)


class ObservationIndex:
    """把三张同代表整理成按位置、过程索引的只读结构。"""

    MAX_ROWS = 2_000_000

    def __init__(self, tables):
        if not isinstance(tables, dict):
            raise ObservationError("供数表必须是对象")
        missing = {"stations", "readings", "villages"} - set(tables)
        if missing:
            raise ObservationError("缺少领域表：" + "、".join(sorted(missing)))
        total = sum(len(tables[name]) for name in ("stations", "readings", "villages"))
        if total > self.MAX_ROWS:
            raise ObservationError("同代资料记录数超过分析上限，请先按需缩小资料范围")
        self.positions = self._positions(tables["stations"])
        self.villages = self._villages(tables["villages"])
        self.readings_by_position, self.readings = self._readings(tables["readings"])
        self.processes = self._processes(self.readings)
        self.process_ids = [process.process_id for process in self.processes]

    @staticmethod
    def _positions(rows):
        positions = {}
        for row in rows:
            try:
                key = (row["station_id"], row["position_version"])
                x = number(row["x"], "stations.x")
                y = number(row["y"], "stations.y")
                start = timestamp(row["start"], "stations.start")
                end = None if row["end"] is None else timestamp(row["end"], "stations.end")
            except KeyError as error:
                raise ObservationError("站点位置缺少字段：" + str(error)) from error
            if not isinstance(key[0], str) or not key[0]:
                raise ObservationError("station_id 不能为空")
            if isinstance(key[1], bool) or not isinstance(key[1], int):
                raise ObservationError("position_version 必须是整数")
            if key in positions:
                raise ObservationError("站点位置版本重复：" + repr(key))
            if end is not None and end <= start:
                raise ObservationError("站点位置有效期不能为空：" + repr(key))
            positions[key] = Position(key[0], key[1], x, y, start, end)
        by_station = {}
        for position in positions.values():
            by_station.setdefault(position.station_id, []).append(position)
        for station_id, versions in by_station.items():
            versions.sort(key=lambda position: position.start)
            for previous, current in zip(versions, versions[1:]):
                if previous.end is None or current.start < previous.end:
                    raise ObservationError("同一站点的位置版本有效期重叠：" + station_id)
        return positions

    @staticmethod
    def _villages(rows):
        villages = {}
        for row in rows:
            try:
                identifier = row["village_id"]
                xmin = number(row["xmin"], "villages.xmin")
                ymin = number(row["ymin"], "villages.ymin")
                xmax = number(row["xmax"], "villages.xmax")
                ymax = number(row["ymax"], "villages.ymax")
            except KeyError as error:
                raise ObservationError("村域缺少字段：" + str(error)) from error
            if not isinstance(identifier, str) or not identifier:
                raise ObservationError("village_id 不能为空")
            if identifier in villages:
                raise ObservationError("村域编号重复：" + identifier)
            if xmax <= xmin or ymax <= ymin:
                raise ObservationError("村域矩形必须非空：" + identifier)
            villages[identifier] = Village(identifier, xmin, ymin, xmax, ymax)
        ordered = sorted(villages.values(), key=lambda village: village.village_id)
        for index, first in enumerate(ordered):
            for second in ordered[index + 1:]:
                if (first.xmin < second.xmax and second.xmin < first.xmax
                        and first.ymin < second.ymax and second.ymin < first.ymax):
                    raise ObservationError("村域矩形不能重叠：" + first.village_id + "、" + second.village_id)
        return villages

    def _readings(self, rows):
        readings = []
        by_position = {}
        identifiers = set()
        for row in rows:
            try:
                identifier = row["reading_id"]
                key = (row["station_id"], row["position_version"])
                start = timestamp(row["start"], "readings.start")
                end = timestamp(row["end"], "readings.end")
                quality = row["quality"]
                process_id = row["process_id"]
            except KeyError as error:
                raise ObservationError("观测记录缺少字段：" + str(error)) from error
            if not isinstance(identifier, str) or not identifier:
                raise ObservationError("reading_id 不能为空")
            if identifier in identifiers:
                raise ObservationError("观测记录编号重复：" + identifier)
            identifiers.add(identifier)
            if key not in self.positions:
                raise ObservationError("观测引用了不存在的位置版本：" + repr(key))
            if not isinstance(quality, str) or quality not in QUALITIES:
                raise ObservationError("quality 非法：" + repr(quality))
            if not isinstance(process_id, str) or not process_id:
                raise ObservationError("process_id 不能为空")
            if end <= start:
                raise ObservationError("观测区间必须非空：" + identifier)
            position = self.positions[key]
            if start < position.start or (position.end is not None and end > position.end):
                raise ObservationError("观测超出位置有效期：" + identifier)
            rain = row["rain_mm"]
            if quality == "valid":
                if rain is None:
                    raise ObservationError("valid 观测必须有雨量：" + identifier)
                rain = number(rain, "readings.rain_mm")
                if rain < 0:
                    raise ObservationError("雨量不能为负：" + identifier)
            elif quality == "missing":
                if rain is not None:
                    raise ObservationError("缺测观测的 rain_mm 必须为 null：" + identifier)
                rain = None
            else:
                rain = None if rain is None else number(rain, "readings.rain_mm")
                if rain is not None and rain < 0:
                    raise ObservationError("雨量不能为负：" + identifier)
            reading = Reading(identifier, key[0], key[1], start, end, rain, quality, process_id)
            readings.append(reading)
            by_position.setdefault(key, []).append(reading)
        for entries in by_position.values():
            entries.sort(key=lambda reading: (reading.start, reading.end, reading.reading_id))
        readings.sort(key=lambda reading: (reading.start, reading.end, reading.reading_id))
        return by_position, readings

    @staticmethod
    def _processes(readings):
        grouped = {}
        for reading in readings:
            grouped.setdefault(reading.process_id, []).append(reading)
        processes = []
        for identifier, entries in grouped.items():
            entries.sort(key=lambda reading: (reading.start, reading.end, reading.reading_id))
            base = min(reading.start for reading in entries)
            processes.append(Process(identifier, base, tuple(entries)))
        processes.sort(key=lambda process: (process.base, process.process_id))
        return processes

    def active(self, position, start, end):
        return position.start <= start and (position.end is None or end <= position.end)

    def position_value(self, key, start, end, process_id=None):
        """返回 (累积雨量或 None, 不可用原因或 None)。

        只有恰好铺满 [start, end) 的连续 valid 观测链才得到数值；更长、更短、
        有缺口或有缺测/可疑记录的链都返回 None，绝不折算成部分累积。

        process_id 给出时只在该过程内寻找观测链：不同过程即使绝对时间相邻或
        重叠，也各是各的降雨过程，不能相互拼接或干扰。
        """
        readings = self.readings_by_position.get(key, ())
        if process_id is not None:
            readings = [reading for reading in readings if reading.process_id == process_id]
        relevant = [reading for reading in readings if reading.start < end and reading.end > start]
        if not relevant:
            return None, "no_observation"
        for previous, current in zip(relevant, relevant[1:]):
            if current.start < previous.end:
                return None, "overlapping_observations"
        if relevant[0].start < start or any(reading.end > end for reading in relevant):
            return None, "interval_mismatch"
        total = 0.0
        cursor = start
        for reading in relevant:
            if reading.end <= cursor:
                continue
            if reading.start > cursor:
                return None, "interval_gap"
            if reading.start < cursor:
                return None, "overlapping_observations"
            if reading.quality == "missing":
                return None, "missing_observation"
            if reading.quality == "suspect":
                return None, "suspect_observation"
            if reading.rain_mm is None:
                return None, "missing_observation"
            total += reading.rain_mm
            cursor = reading.end
            if cursor >= end:
                break
        if cursor != end:
            return None, "interval_gap"
        return total, None
