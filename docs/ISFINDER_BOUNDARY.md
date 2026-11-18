# 边界检验：对照 ISfinder 的已知 IS 元件

2026-11-18。独立的检验模块，只读 catalogue，不改任何产物。代码 `scripts/260_isfinder_boundary.py`，测试 `tests/test_isfinder_boundary.py`，元件名单 `config/isfinder_panel.tsv`。未列入 SPEC，未经审查循环。`$W` = `/global/scratch/users/kh36969/fna_ins_discovery`。

## 问题

`86_catalogue.py` 的文件头写明"junction 是否精确到碱基从未检验"。这里检验：catalogue 中的插入序列若就是一个已知的、人工整理过的 IS 元件，我们切出的两端是否正好落在元件的两端。

## 数据

- 参照：ISfinder 序列的公开镜像（GitHub `thanhleviet/ISfinder-sequences`，`IS.fna`，提交 `d5b1102`，2020-10-06；ISfinder 本身不提供整库下载）。存于 `$W/isfinder_bench/ref/IS.fna`，sha256 `47607f9e…d5507`。5,970 条记录，其中 400 条是同名同序列的重复，去重后 5,570 条。
- 被检验者：`$W/catalogue_v3/*_events.tsv` 的全部 109,709 个事件的 `insert_seq`（18 个物种标签）。

## 方法

每个插入序列与全部 IS 比对（minimap2 `-c --eqx -x asm20`），取匹配碱基最多的一条；identity ≥ 0.95 且覆盖 IS ≥ 90% 的才判定。把 IS 的两端投影到插入序列上（比对器在末端没对上的 IS 碱基按共线推算），得到插入序列在元件两侧多出的碱基数 `left_extra`、`right_extra`。

catalogue 的插入序列是空等位两侧 flank 之间的那一段：junction 是正向重复（TSD）时，这一段是"元件 + 一份 TSD"，TSD 可以在任一端或分在两端。catalogue 自己报告了重复长度 `junction_overlap_bp`。所以完整元件的精确切割应满足 `extra = left_extra + right_extra = junction_overlap_bp`。

| 类别 | 含义 |
|---|---|
| `exact` | extra 等于 junction 重复，两端都不切进元件 |
| `exact_shifted` | junction 重复比 TSD 长 1–3 bp：元件端部碱基碰巧与 TSD 之后的碱基相同（每端约 1/4 概率），切口是同一等位的等价写法。不是边界错误 |
| `near` | 其余相差 ≤ 3 bp |
| `over` | 插入序列比"元件 + 重复"多 > 3 bp |
| `under` | 插入序列在元件内部结束 |

`over` 与 `under` 不等于错误：可以是元件加携带序列、两个元件，或真实的截短拷贝。

## 结果（运行 `$W/isfinder_bench/runs/r2/`，登录节点 8 线程 1 分钟）

109,709 个插入中 71,538 个（65.2%）是 ISfinder 中的元件（按上述阈值）。

| 类别 | 数目 | 占判定数 |
|---|---:|---:|
| exact | 37,644 | 52.6% |
| exact_shifted | 21,288 | 29.8% |
| near | 1,930 | 2.7% |
| over | 8,500 | 11.9% |
| under | 2,176 | 3.0% |

**边界正确（exact + exact_shifted）：58,932 个，82.4%。**

独立的佐证：名单中凡填了文献 TSD 长度的元件，实测的众数 extra 都与文献一致——IS1 9、IS2 5、IS3 3、IS4 12、IS5/ISKpn26 4、IS10 9、IS26 8、IS50 9、IS903 9、IS256 8、ISAba1 9、IS1181/IS1251/IS16 8、ISKpn1 3、IS91 0。多数元件该众数占 94–100%（IS1 79%，IS4 51%，这两个元件的 TSD 长度本来就不固定）。catalogue 没有用任何 IS 知识，这些长度是从组装比较中得到的。

`exact_shifted` 的解释经过检验：21,388 个 delta 为 −1 至 −3 的事件中 21,288 个符合"切口在重复允许的范围内平移"；它们的 extra 等于该元件众数 TSD 的比例（89.3%）与 `exact`（90.6%）相同。

按 catalogue 自己的标记分：`S1_structurally_clean = True` 的 67,292 个中 84.3% 正确，`False` 的 4,246 个中 52.6%。`exact` 分解 85.7%，`tolerant_alignment` 75.5%。

按物种（判定数 ≥ 1,000 的）：E. faecalis 91.8%，S. enterica 85.9%，A. baumannii 84.9%，E. coli 84.5%，E. faecium 84.2%，K. pneumoniae 83.8%，S. aureus 78.5%（k12）/ 75.5%（k24），P. aeruginosa 74.9%，S. pneumoniae 67.1%，M. tuberculosis 65.1%。

### 不精确的部分

- **over（8,500）**：3,922 个多出 > 1 kb，1,465 个多出 101–1000 bp，1,202 个 21–100 bp，1,911 个 4–20 bp。1,445 个插入中还有第二个 IS。5,322 个只在一侧多出（另一侧在重复范围内），即元件的一端是对的。已知的生物学例子：ISEc9（ISEcp1 类）801 个中 575 个是 over——这类元件以单端转座带走相邻序列。其余 over 是元件加携带序列还是边界偏移，本模块不能区分。
- **under（2,176）**：1,641 个只短 4–20 bp。截短拷贝与切进元件，本模块不能区分。
- **near（1,930）**：部分来自参照本身。例：IS421 的 755 个没有一个 exact，几乎全部差 2 bp（比对两端各有 1 个 IS 碱基对不上）；IS102 的 1,359 个中 1,100 个 extra 为 8，而它的同型 IS903B 是 9。这是 ISfinder 记录的端点与实际拷贝差 1–2 bp，不是我们的切口错。

## 不精确部分的归因（2026-11-18，一次性分析，未写成脚本）

把判定的 71,538 个逐个归到一类（`per_event.tsv` 与 catalogue 的 `target_seq`、`insertion_point_offset` 连接后，看元件两侧的碱基）：

| 类 | 数目 | 占比 | 是不是边界错误 |
|---|---:|---:|---|
| 正确（exact + exact_shifted） | 58,932 | 82.4% | 否 |
| 元件 + 重复，重复旁有 1–2 bp 的多态（SNP 或同聚物长度差） | 1,642 | 2.3% | 否；切口对，`junction_overlap_bp` 报短了 |
| 元件 + ≤ 25 bp，看不出重复 | 1,659 | 2.3% | 未知 |
| 元件在更大插入的一端 | 2,743 | 3.8% | 一端对；另一端之外还有别的序列 |
| 元件在更大插入的内部 | 2,734 | 3.8% | 否（插入是更大的东西，IS 嵌在里面） |
| ISfinder 记录的末端 ≥ 4 bp 对不上 | 3,062 | 4.3% | 无结论（检验的限制） |
| 插入在元件内部结束 | 766 | 1.1% | 未知（截短拷贝或切进元件） |

可能是真正边界错误的上限是"未知"两类加"一端"类，共 7.2%；确定不是错误的至少 88.5%。物种间的差别主要来自后几类：L. monocytogenes 69% 是"元件在更大插入内部"，M. tuberculosis 与 P. aeruginosa 各 13%，S. pneumoniae 20% 是参照末端对不上。

## 限制

- 只检验了是已知 IS 的插入（65%）；其余 38,171 个没有参照，本检验对它们无结论。
- 最佳命中在近似同型（IS1A / IS1R / IS1X2，IS6110 / IS987，IS629 / IS1203）之间的归属是任意的；按名字读单个元件的数目时要合并同型。名单里没出现在结果中的元件（如 IS6110、IS481、IS150、IS30、IS200）是被同型吸收或没有判定的插入。
- 末端投影假设没对上的 IS 末端与插入序列共线；IS 记录末端与实际拷贝有分歧时（IS102 末 53 bp）投影有 1–2 bp 的不确定。
- `config/isfinder_panel.tsv` 中的文献 TSD 是手工录入的常识值，不是从下载文件读出的。
- 检验的是 catalogue（Layer 2 的产物），不是 Arm A 的家族边界。

## 运行

```
export PATH=$PWD/env/bin:$HOME/.conda/envs/claude-env/bin:$PATH
python scripts/260_isfinder_boundary.py --isfinder $W/isfinder_bench/ref/IS.fna \
  --isfinder-sha256 47607f9e398d200f1954b18b9d91d6f5a612331ee137cb10d63897688f9d5507 \
  --events "$W/catalogue_v3/*_events.tsv" --panel config/isfinder_panel.tsv \
  --out $W/isfinder_bench/runs/<新目录> --threads 8
```

产物：`per_event.tsv`（每个判定的插入一行）、`per_is.tsv`（每个元件一行）、`result.json`（输入 sha256、提交、参数、各类计数）。
