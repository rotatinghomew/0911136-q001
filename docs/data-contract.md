# 领域资料与供数约定

## 数据含义

坐标为局地等面积平面米，村域为不重叠轴对齐矩形。雨量为 start 含、end 不含区间内累积毫米，不是瞬时强度；不同区间不能直接相加或比较。station_id 位置随 position_version 变化；每个 reading 显式引用对应位置版本。quality 只允许 valid、suspect、missing；缺测 rain_mm 为 null，0 为有效零雨。process_id 表示同一连续降雨过程，不等于独立样本。

## 现行摘要

`rainnet/analysis.py` 生成现行资料摘要。源码和运行结果中的 scope 说明其统计边界，不能把未提供的推断当成已有能力。

## 供数接口

`GET /manifest` 返回 generation、tables（每表 count 与 sha256）及 metadata。`GET /tables/<name>?generation=<generation>` 返回同代 rows。供数换代时旧 generation 返回 HTTP 409，不允许将新代记录接入旧分析；不存在的表返回 404。表内容为稳定键唯一的完整同代记录，不是增量事件流。客户端对清单摘要进行核验。摘要使用 UTF-8、排序对象键、紧凑 JSON 的 SHA256，数组顺序保留供数端顺序。

字段与主键声明位于 `rainnet/schema.json`，合成资料位于 `fixtures/dataset.json`。`metadata.cutoff` 是供数展示截点，不自动为分析选择训练截止点。服务可将整个 dataset 原子替换为新一代，但不在同一 generation 下更改内容。

## 结果保存

`reports.py` 保存现有摘要及来源表摘要，report_id 是内容摘要，不是随机编号。相同内容生成相同编号，临时文件写入完成后原子替换目标文件；历史内容是否保留由调用者选择不同目标路径。服务只监听本机，不接入生产系统。
