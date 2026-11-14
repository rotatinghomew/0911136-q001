"""独立核对：扩展分析是否覆盖需求中的每个判断点。

不是单元测试的替代，而是把最终输出逐条对照需求，避免只验证实现内部结构。
运行：python3 -m tests.check_requirements
"""

import json
import sys
from pathlib import Path

from rainnet.analysis import analyze_area, summarize
from rainnet.datasource import Snapshot
from rainnet.reports import package_area_report

ROOT = Path(__file__).resolve().parents[1]


def main():
    tables = json.loads((ROOT / "fixtures/dataset.json").read_text())["tables"]
    analysis = analyze_area(tables)
    checks = []

    def check(name, condition, detail=""):
        checks.append((name, bool(condition), detail))

    villages = analysis["villages"]
    # 1. 点观测与面积口径分开。
    check("点/面口径分开", "point_scope" in analysis and "basis" in analysis)
    # 2. 有站区域 vs 不能推断的村。
    check("无站村不可推断", villages["v3"]["status"] == "not_inferable")
    check("可估计村标注覆盖状态",
          all(v["status"] in {"estimated", "estimated_single_station",
                              "estimated_partial_coverage", "estimated_partial_periods"}
              for k, v in villages.items() if k != "v3"))
    # 3. 单站不等同整村。
    check("单站支撑有独立状态",
          all(v["status"] != "estimated" for k, v in villages.items()
              if any(p["single_station_support"] for p in v["processes"])))
    # 4. 缺测不是零雨。
    missing = [(v, p) for v, item in villages.items() for p in item["processes"]
               if any(s["reason"] == "missing_observation"
                      for loc in p["locations"] for s in loc["unusable_stations"])]
    check("缺测明确为不可用", all(p["coverage_status"] not in {"estimated"}
                                 or p["coverage"] < 1.0 for _, p in missing),
          f"{len(missing)} 个村-过程含缺测")
    # 5. 站点迁移不被合并。
    moved = analysis["network"]["station_changes"]["moved_stations"]
    check("迁移站分版本",
          moved and all(len(item["versions"]) >= 2 for item in moved))
    # 6. 区间不齐有候选口径与说明。
    check("口径候选有支撑比较",
          len(analysis["basis"]["alternatives"]) >= 1
          and all("covered_fraction_sum" in a for a in analysis["basis"]["alternatives"]))
    check("改口径的可用性已量化",
          len(analysis["basis_sensitivity"]["alternatives"]) ==
          len(analysis["basis"]["alternatives"]))
    # 7. 站网变动是否改变结论。
    decision = analysis["network"]["decision"]
    check("站网变动有结论", decision["verdict"] in {
        "station_network_robust", "station_network_changes_ranking",
        "station_network_changes_conclusions"}, decision["verdict"])
    # 8. 站非独立样本。
    dependence = analysis["dependence"]
    check("非独立说明",
          "不是独立样本" in dependence["explanation"]
          and dependence["correlation_model"] is not None)
    check("设计效应可解释", dependence["design_effect_by_village"])
    # 9. 局地极值、区域分布、区间估计。
    check("区域分布存在", analysis["regional"]["window_distribution"] is not None)
    check("面积极值有边界",
          analysis["regional"]["area_extreme"] is None
          or "不是村内任意点" in analysis["regional"]["area_extreme"]["scope"])
    check("点极值单列", analysis["regional"]["point_extreme"] is not None)
    check("逐过程有区间",
          all(p["interval95_mm"] is not None or p["covered_area_mean_mm"] is None
              for v in villages.values() for p in v["processes"]))
    # 10. 现有按站摘要仍可用。
    check("按站摘要保持原字段",
          set(summarize(tables)) == {"point_totals", "village_count", "scope"})
    # 11. 分析报告与摘要共用同一代快照和来源证据。
    generation = json.loads((ROOT / "fixtures/dataset.json").read_text())["generation"]
    snapshot = Snapshot(generation, tables, {"stations": {"count": len(tables["stations"])}}, {})
    area_report = package_area_report(snapshot, analysis)
    check("面积报告带同代来源",
          area_report["generation"] == generation
          and area_report["sources"] == snapshot.evidence)

    width = max(len(item[0]) for item in checks)
    failed = 0
    for name, ok, detail in checks:
        print(f"{name:<{width}} {'通过' if ok else '未通过'} {detail}")
        failed += 0 if ok else 1
    print()
    print(f"{len(checks) - failed}/{len(checks)} 项通过")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
