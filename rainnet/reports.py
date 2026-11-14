"""原子保存分析摘要和来源证据。"""

import json
import os
import tempfile
from pathlib import Path

from .datasource import canonical, digest


def package_report(snapshot, summary):
    content = {"generation": snapshot.generation, "sources": snapshot.evidence,
               "summary": summary, "metadata": snapshot.metadata}
    return {"report_id": digest(content), **content}


def package_area_report(snapshot, summary):
    """面积口径报告：单独打包，现有按站摘要报告保持原样可用。"""
    content = {"kind": "village_area_analysis", "generation": snapshot.generation,
               "sources": snapshot.evidence, "analysis": summary,
               "metadata": snapshot.metadata}
    return {"report_id": digest(content), **content}


def verify_report(report):
    body = {key: value for key, value in report.items() if key != "report_id"}
    if digest(body) != report.get("report_id"):
        raise ValueError("摘要内容已变化")
    return report


def save_report(path, report):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".report-", dir=target.parent)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(canonical(report) + b"\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def load_report(path):
    report = json.loads(Path(path).read_text(encoding="utf-8"))
    identifier = report.pop("report_id")
    if digest(report) != identifier:
        raise ValueError("摘要内容已变化")
    return {"report_id": identifier, **report}


def load_area_report(path):
    report = json.loads(Path(path).read_text(encoding="utf-8"))
    if report.get("kind") != "village_area_analysis":
        raise ValueError("不是村域面积分析报告")
    return verify_report(report)
