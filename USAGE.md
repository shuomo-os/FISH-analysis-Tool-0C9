# 修复版使用说明

本项目保留原有的探针分析、RNA 探针设计和反向互补三个页面。原作者声明及许可证见 [README.md](README.md)，该文件内容保持原样。

## 安装与启动

建议使用 Python 3.10 或更高版本，下载整个仓库，将 `0C9.py`、`probe_io.py`、`blast_support.py` 和原始 `README.md` 放在同一目录。

```bash
python -m pip install -r requirements.txt
python 0C9.py
```

Python 需要带有 Tkinter；Linux 系统可通过发行版的 Python Tk 软件包安装。使用本地 BLAST 时还需安装 NCBI BLAST+，并准备核酸数据库。

现在可以从其他工作目录使用脚本的绝对路径启动。README 校验按脚本所在目录定位，接受 Git 在 Windows 上使用的 CRLF 换行；仍会校验许可证正文。导入 Python 模块供测试使用时不会打开界面，直接运行脚本时仍执行许可校验。

## 表格输入与输出

| 扩展名 | 读写格式 |
| --- | --- |
| `.csv` | 逗号分隔 |
| `.tsv`、`.txt` | 制表符分隔 |
| `.xlsx` | Excel，使用 openpyxl |
| 其他扩展名 | 默认按 CSV 处理 |

探针分析必须有 `sequence` 列；建议另有 `id` 列。如果没有 `id`，程序按行生成 `probe_1`、`probe_2` 等名称。示例：

```csv
id,sequence
001,CGTTCCAAAGATGTGGGCATGAGCTTAC
sample A,CTGCACTGAATCCGTAGCAT
sample A,ACGTNACGT
```

ID 按文本读取，保留前导零、空格和重复值。Excel 输入中希望保留前导零时，应将该单元格存为文本；已经被 Excel 存成数值的原始格式不能由程序恢复。

缺失序列保持为空，不会转成 `nan`。分析输出新增以下列：

- `input_sequence`：原始序列单元格内容。
- `sequence`：去掉首尾空白并转成大写后的序列。
- `alphabet_valid`：是否只含 DNA 字母 A/T/C/G/N，且非空。
- `valid_sequence`：是否为至少两个碱基、只含 A/T/C/G 的序列。
- `tm_calculable`：所选方法是否成功给出 Tm。
- `sequence_error`：空值、未定碱基、非法字符或无法计算的原因。
- `tm_method`：本次使用的 Tm 方法。

含 N 的序列不计算 Tm/GC，也不进入 BLAST，避免软件删掉 N 后把较短序列的计算值当成原序列结果。含内部空格或其他字符的分析序列会明确标记无效。表格保存先写临时文件，写入成功后替换目标，失败时保留旧文件。

## Tm 模型

| 方法 | 计算模型 | 适用输入 |
| --- | --- | --- |
| `santalucia` | Allawi & SantaLucia 1997，`DNA_NN3` | DNA/DNA |
| `nn` | Biopython 默认最近邻表，目前也是 `DNA_NN3` | DNA/DNA |
| `rna_dna` | Sugimoto 1995，`R_DNA_NN1` | RNA/DNA，仅设计页面 |
| `wallace` | Wallace 经验公式 | DNA 探针 |
| `gc` | Biopython `Tm_GC` 默认经验公式 | DNA 探针 |

原 `santalucia` 选项错误使用了 Breslauer 1986 的 `DNA_NN1`。修复后，相同序列的 Tm 和通过筛选的候选数量可能改变。例如 Biopython 文档中的 `CGTTCCAAAGATGTGGGCATGAGCTTAC`，在默认条件下，`DNA_NN3` 为 60.32°C，原 `DNA_NN1` 为 72.19°C。

默认仍是 DNA/DNA 的 `santalucia`。设计 DNA 探针与 RNA 靶标的杂交时，可以选择 `rna_dna`；程序使用目标 RNA 片段的 5′→3′方向计算该模型，再输出对应 DNA 反向互补探针。

最近邻方法沿用 Biopython 默认计算条件：两条链各 25 nM、Na 50 mM、K/Tris/Mg/dNTPs 为 0、`selfcomp=False`、`saltcorr=5`。界面尚未开放这些条件，也没有加入甲酰胺、DMSO 或具体杂交缓冲液的修正，因此报告的是指定模型及默认条件下的 Tm。

模型与参数依据：[Biopython MeltingTemp 文档](https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html)。

## RNA 探针设计

- 目标序列按 5′→3′输入，支持 A/T/C/G/U，输出 DNA 探针也按 5′→3′显示。
- 加载多个记录的 FASTA 时，程序要求选择一条目标，分别保留记录边界，不会拼接记录。
- “探针间空余碱基数”表示相邻两个已接受探针之间未覆盖的碱基数量。例如探针长度 20、间距 3，对应位置 1–20、24–43。设置 0 时允许紧邻。
- 最小复杂度范围为 0–1。重复过滤能够识别重叠的 6-mer，以及同一完整候选在目标中的多个位置。
- 最后一个完整候选窗口会被检查；目标长度恰好等于探针长度时也可以得到候选。
- 设计完成后点击“运行BLAST检查”才会执行数据库搜索。详情与位置图始终使用本次设计的目标序列快照，后续编辑输入框不会改变已有探针的位置解释。

## BLAST 设置与结果

“程序路径”选择 `blastn.exe`（Windows）或 `blastn`。数据库输入应为包含数据库名称的**前缀**，例如 `E:\Blast_db\refseq_rna`，不能只填目录。浏览按钮可选择 `.nal`、`.nin`、`.nsq` 等索引文件，程序会提取前缀。

两个页面统一使用适合短核酸的参数：

```text
-task blastn-short -word_size 7 -evalue 1000
-dust no -soft_masking false -strand both -max_target_seqs 1000
```

这些设置用于保留短序列的潜在匹配，E-value 阈值本身不代表已确认的生物学特异性。依据：[NCBI 短序列检索说明，表 4.3.1](https://blast.ncbi.nlm.nih.gov/blast/BLAST_guide.pdf)。比对输出采用明确的 14 个字段，包括序列长度、比对长度、查询覆盖率和坐标，见 [NCBI 自定义输出说明](https://www.ncbi.nlm.nih.gov/sites/books/NBK569862/)。

每个查询使用内部 ID（如 `q000001`），按原始行位置对应，用户自己的 ID 不作为关联键。原始 BLAST 文件旁保存 `<BLAST文件名>.queries.csv`，包含内部 ID、原始行号、原 ID 和序列；保存失败会回滚已替换的对应表，避免它和旧 BLAST 文件混用。不同任务不能同时占用相同输入/输出文件。

分析页会在分析结果旁生成 `_with_blast` 文件，保留所选 CSV/TSV/XLSX 格式。设计页会更新表格，导出结果中的 `blast_hits` 使用 JSON 保存详细比对数据。

| 结果字段 | 含义 |
| --- | --- |
| `blast_hits_count` | 已返回结果中的不同 subject ID 数，不是不同基因数 |
| `blast_hsp_count` | 比对片段（HSP）数量，同一 subject 可以有多条 |
| `blast_full_match_count` | 查询全长、100%一致、无错配和无缺口的 subject 数 |
| `blast_hit_limit_reached` | 是否达到当前 1000 条返回上限 |
| `blast_report_limit` | 本次设置的最大返回 subject 数 |
| `blast_query_id` | 本次检索的内部查询 ID |
| `specificity` | 保留旧字段名，内容改为描述性匹配状态 |

界面显示“单一序列命中（待核实）”“多序列命中（待核实）”或“无匹配（当前参数）”。没有指定预期靶标和基因归属时，不自动给出“高特异性”结论。返回数量有上限，搜索也受数据库版本和参数影响，输出计数不能解释为数据库中的全部潜在匹配。

不足 7 nt 的分析序列和无效序列显示为“未检索”，不与“已检索但没有命中”混淆。空文件、缺失查询和截断输出会报错。重新运行时清空旧 BLAST 证据；失败或取消会显示相应状态。

BLAST 在后台线程执行，两个页面都可取消，取消会终止并回收子进程。分析页“暂停”用于逐条探针分析；进入 BLAST 阶段后禁用暂停，保留取消。关闭窗口也会先停止后台任务。

## 回归验证

```bash
python -m unittest discover -s tests -v
```

测试覆盖 Tm 参数与链方向、窗口边界、间距、重叠重复、CSV/TSV/XLSX、ID 映射、多记录 FASTA、结果重跑、进程取消、失败保存回滚及许可路径/换行。安装了 `blastn` 和 `makeblastdb` 时，还会自动创建临时小型数据库，执行真实 20 nt 查询；未安装时仅该项跳过。

本次验证环境为 Python 3.12、Biopython 1.88、pandas 2.2.3、NumPy 2.3.5、openpyxl 3.1.5 和 BLAST+ 2.12.0。界面的线程入口和结果处理已通过无显示回归测试；当前测试环境无法启动显示服务，尚未进行 Windows 图形界面的完整交互验证，也未使用实验室的实际大型数据库。
