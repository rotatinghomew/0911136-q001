"""按站位置与观测区间整理有效雨量，保持点位口径。"""
from collections import defaultdict


def station_totals(tables):
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
    return {"point_totals": station_totals(tables), "village_count": len(tables["villages"]),
            "scope": "按站位置累积摘要，不表示村域面积平均"}
