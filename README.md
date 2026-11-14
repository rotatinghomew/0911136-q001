# 不均匀站网下的村域降雨分布与空间支撑

本项目组织领域资料并提供当前摘要，通过本地 HTTP 供数服务连接应用。业务代码在 `rainnet`，字段与单位约定见 `docs/data-contract.md`。示例均为本地合成记录，不代表真实县域记录或业务结论。

## 环境

Python 3.11 及以上，业务运行仅使用标准库。构建工具及其依赖完整固定在 `requirements-build.txt`；先用 `python3 -m pip install -r requirements-build.txt` 安装。源码和业务测试无需访问外部服务，首次安装构建依赖需要访问 Python 包仓库。

## 当前摘要

在项目根目录运行，会临时启动真实本地 HTTP 服务、核验资料代际、计算现有摘要并完整打印 JSON：

```sh
python3 -m rainnet.runner --output /tmp/rainnet-summary.json
```

`point_totals` 就是原有按站摘要，字段与 `scope` 说明保持不变，调用方式也不变。

## 村域面积分析

在同一代资料上另算村域分析，并单独保存为可核验报告；同一次运行里现行摘要仍写到 `--output`：

```sh
python3 -m rainnet.runner --output /tmp/rainnet-summary.json --area-output /tmp/rainnet-area.json
```

输出把点观测与面积口径分开：

- `point_scope`：各位置的原生区间累积，保留现有按站摘要口径；`intervals_consistent`、`valid_intervals_contiguous` 说明区间是否可比。
- `basis`：跨站比较用的"过程相对口径"，按可支撑村域覆盖面积自动选择，并列出其他候选口径的支撑变化。
- `villages`：每村逐过程的 `coverage`（有可用观测的村域面积占比）、`covered_area_mean_mm`（该覆盖面积上的面积加权均值）、`coverage_status` 与 `reasons`。`not_inferable` 的村没有数值；只有单站支撑的村标为 `estimated_single_station`。
- `regional`：只在可估计村之间形成的分布、面积口径极值和独立的点观测极值。
- `dependence`：同过程站间差异给出的空间方差、站间距离相关模型、设计效应，解释为什么站点不能当独立样本。
- `network`：迁移、停测、晚启用位置的记录，以及排除变动站位后结论是否改变。
- `basis_sensitivity`：换用其他累积口径时哪些村还能推断。

缺测、可疑、区间不齐、停测和迁移都不折算成零雨，也不插值或跨位置版本补值。区间估计是空间支撑近似，不随过程数收窄。

## 分开运行供数与应用

```sh
python3 -m rainnet.runner --serve --port 8765
```

另一个终端可运行：

```sh
python3 -m rainnet.runner --source-url http://127.0.0.1:8765 --output /tmp/rainnet-summary.json --area-output /tmp/rainnet-area.json
```

## 测试与检查

```sh
python3 -m unittest discover -s tests -t . -v
python3 -m tests.check_requirements
python3 -m compileall -q rainnet tests
python3 -m build --wheel --no-isolation
```

测试覆盖现行领域摘要、非法引用、真实 HTTP 换代拒绝、表清单内容核验、结果保存，以及面积分析的覆盖状态、缺测/可疑处理、区间不齐、站点迁移、站网敏感性和设计效应。展示现有产品可直接运行上述摘要命令，不需要回放历史日志；新增能力用 `--area-output` 运行，原有按站摘要调用保持不变。

## 容器运行

Dockerfile 固定官方 Python 3.11.14 slim-bookworm 镜像摘要，不依赖本机私有镜像；首次构建需下载基础镜像和固定构建工具。构建与运行实际摘要：

```sh
docker build -t local/rainnet .
docker run --rm local/rainnet
```

生成面积分析报告：

```sh
docker run --rm -v "$PWD/out:/out" local/rainnet python3 -m rainnet.runner \
  --output /out/summary.json --area-output /out/area.json
```

容器内测试：

```sh
docker run --rm --entrypoint python3 local/rainnet -m unittest discover -s tests -t . -v
```
