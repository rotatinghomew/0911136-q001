"""运行领域摘要、村域面积分析或启动本地供数服务。"""

import argparse
import json
from http.server import ThreadingHTTPServer
from pathlib import Path

from .analysis import analyze_area, summarize
from .datasource import SourceClient
from .reports import package_area_report, package_report, save_report
from .service import handler_for, load_repository, serve


def main():
    parser = argparse.ArgumentParser(description="领域资料摘要与村域降雨分析")
    parser.add_argument("--dataset", default="fixtures/dataset.json")
    parser.add_argument("--source-url")
    parser.add_argument("--output", default="/tmp/domain-summary.json")
    parser.add_argument("--area-output",
                        help="把村域面积分析另存为独立 JSON；默认摘要仍保存在 --output")
    parser.add_argument("--serve", action="store_true")
    parser.add_argument("--port", type=int, default=8765)
    arguments = parser.parse_args()
    schema_path = Path(__file__).with_name("schema.json")
    if arguments.serve:
        server = ThreadingHTTPServer(("127.0.0.1", arguments.port), handler_for(load_repository(arguments.dataset, schema_path)))
        try:
            print(json.dumps({"listening": server.server_address}, ensure_ascii=False), flush=True)
            server.serve_forever()
        finally:
            server.server_close()
        return
    if arguments.source_url:
        snapshot = SourceClient(arguments.source_url).read()
    else:
        with serve(load_repository(arguments.dataset, schema_path)) as base_url:
            snapshot = SourceClient(base_url).read()
    # 一次快照同时供现行摘要与面积分析使用，避免不同代资料混入同一次输出。
    summary = summarize(snapshot.tables)
    report = package_report(snapshot, summary)
    save_report(arguments.output, report)
    print(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False))
    if arguments.area_output:
        area_report = package_area_report(snapshot, analyze_area(snapshot.tables))
        save_report(arguments.area_output, area_report)
        print(json.dumps(area_report, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
