# 候选选择说明

本批提供 VQ-VAE Top100 和 DiMA Top100，合并200条仅用于研究选择。
两份50,000候选库相互独立，不将架构融合为一个提交。

## 直接浏览
- [飞书候选选择文档](https://xcne9cbnfv9b.feishu.cn/wiki/N8aJwWSXaiE29Jk6oSdcq4kwnYg)
- [Top200 CSV](../reports/top200-candidates.csv)：可筛选、保留原始排名。
- [Top200 JSON](../reports/top200-candidates.json)：完整精度。
- [两个50k库及Top100 FASTA](../candidates/)。
- [逐指标附表](evaluation-appendix.html)：下载HTML后浏览；包含所有24个RL run。
- [彩色序列浏览稿](selection-review.html)：下载HTML后浏览，40位位置对齐，不是多序列比对。

## 主要观察

| 指标 | VQ-VAE | DiMA |
| --- | ---: | ---: |
| Library总数 / 唯一数 | 50000 / 50000 | 50000 / 50000 |
| 全库训练集精确命中 | 0 | 34 (0.068%) |
| Top100训练集精确命中 | 0 | 0 |
| Top100训练近邻ratio>0.8 | 0 | 2 |
| Top100官方参考ratio>0.8 | 0 | 0 |
| Top100综合代理均值 | 0.84494 | 0.85013 |
| Top100选择性代理均值 | 0.64632 | 0.60129 |
| Top100预测HC50均值 (µM) | 66.426 | 49.206 |
| Top100预测HC50中位数 (µM) | 58.194 | 38.224 |
| Top100平均长度 | 26.83 | 29.01 |

两个library间无精确重叠；Top200也全部唯一。
DiMA活性代理略高，VQ的HC50/选择性代理较高；这不是实验胜负。

## 每条如何读

candidate_id是稳定的路线+原始rank，不是重新混排。
reward_broad / gram_positive / gram_negative / mdr / selectivity / joint
均为固定校准的0–1代理。MDR列没有MDR菌株专属性。
selection_utility_joint_minus_diversity是原始贪心选择时的J−0.15×近邻惩罚。
ania_*为原始log-MIC输出，不能直接拿它判断官方16µM门槛。
hemopi2_hc50_um是预测，不是实测。
NOT_MEASURED / UNAVAILABLE / NOT_ASSESSED表示缺失，不是零或失败。

200条均通过字符、长度、所属library、唯一性和官方参考阈值检查。
没有任何一条可被称为“所有生物性质已满足”。

## 竞赛逐项选择

| 小项 | 当前可参考列 | 必须补的证据 |
| --- | --- | --- |
| 广谱 | reward_broad及EC/PA/SA原始预测 | 20菌株MIC、成功率、MIC90 |
| Gram+ | reward_gram_positive及SA预测 | 5菌株MIC与MIC50 |
| Gram− | reward_gram_negative及EC/PA预测 | 15菌株MIC与MIC50 |
| MDR | 仅作三物种鲁棒性参考 | 8 MDR分离株对应oracle/实验 |
| 最佳选择性 | 同看HC50、活性与选择性代理 | 至少一菌株MIC≤16µM；20菌株MIC50；实测HC50/MIC50 |
| 多样性/新颖性 | near-reference、near-train、CD50/CD95 | 完整官方seqme及数据库筛查 |
| 可合成性/线性自由末端 | 当前仅设计约束 | 合成成功、化学身份、纯度及末端状态确认 |

全库功能预测没有跑满100,000条：每路线只对固定2,048候选池做外部打分。
原始基线/综合RL/混合原始/过滤库/Top100不能混在一个均值里。
CD50/CD95基于2,048样本或Top100；ProtT5基于512或100，均非全50k逐条。
全量统计与抽样统计在文件中分别标明。
