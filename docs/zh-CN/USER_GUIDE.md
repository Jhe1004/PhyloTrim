# PhyloTrim 用户手册（简体中文）

[English Version](../USER_GUIDE.md) | [返回中文主页](../../README.zh-CN.md)

---

## 1. 支持的文件格式与压缩

PhyloTrim 自动根据文件内容识别比对格式（自动去除 UTF-8 BOM），并原生支持 `.gz` 压缩读写：

| 格式 | 目录模式自动扫描的后缀 | 说明 |
| :--- | :--- | :--- |
| **FASTA** | `.fa`, `.fas`, `.fasta`, `.fna`, `.ffn`, `.faa`, `.frn`, `.aln` | 支持多行或单行序列；若行长不一致会自动用 `-` 右补齐。 |
| **PHYLIP** | `.phy`, `.phylip` | 同时支持 Sequential 与 Interleaved 两种变体。 |
| **NEXUS** | `.nex`, `.nexus` | 支持标准 `DATA` / `CHARACTERS` 区块、带引号物种名与交错矩阵。 |
| **Gzip** | 上述任意后缀 + `.gz` | 读取时通过 gzip 魔数自动识别；目录模式下保持原 `.gz` 状态，单文件模式下按输出文件名是否以 `.gz` 结尾决定是否压缩。 |

### 单文件模式下的格式转换

单文件模式下，程序会根据输出路径扩展名自动转换格式（未识别扩展名则沿用输入格式）：

```bash
# 直接读取 .fas.gz，按 20% 阈值修剪后输出为 PHYLIP 格式
phylotrim -i gene1.fas.gz -o gene1.trimmed.phy -g 0.2
```

目录模式下，输出目录中**仅包含修剪后的比对文件**（不写入额外统计表或日志文件），方便下游建树脚本直接通配读取。

---

## 2. 阈值与边界规则

设输入比对共有 `n_taxa` 条序列，设定阈值为 `max_gap` $\in [0.0, 1.0]$：

- 某列非标准字符数为 `gap_count`，该列被**保留**当且仅当：
  $$\frac{\text{gap\_count}}{n_{\text{taxa}}} \le \text{max\_gap}$$
- **零容忍（`-g 0.0`）**：仅保留 `gap_count == 0` 的完全无缺失列。
- **全保留（`-g 1.0`）**：保留所有列，自动触发零拷贝直出快径。
- **绝不删行**：即使某条序列修剪后全部为 `-`（甚至所有列都被删完导致长度为 `0`），所有物种名与行顺序依然完整保留。

---

## 3. 密码子模式（`-c / --codon`）

开启 `--codon` 后：

1. 输入必须为核酸（DNA/RNA）比对；若误用于蛋白质比对，程序会直接报错拦截。
2. 序列末尾不足 3 bp 的残余位点（`n_sites % 3`）会自动忽略，保证输出长度严格为 3 的整数倍。
3. 对每组三联体 `(3k, 3k+1, 3k+2)`，只有当 **3 列全部满足** `gap_count <= max_allowed_gaps` 时才整组保留；只要其中任意 1 列超标，整组 3 列一同删除。

---

## 4. 常见建树流程衔接

### 多基因树流程（MAFFT → PhyloTrim → IQ-TREE → ASTRAL）

```bash
# 1. 并行修剪全部比对基因
phylotrim -i 01_mafft_loci/ -o 02_trimmed_loci/ -g 0.2

# 2. 批量使用 IQ-TREE 构建基因树
mkdir -p 03_genetrees
for f in 02_trimmed_loci/*.fas; do
  iqtree2 -s "$f" -m MFP -B 1000 -T AUTO --prefix "03_genetrees/$(basename "$f" .fas)"
done
```

### 编码序列（CDS）流程

```bash
# 严格保持阅读框，适用于下游 PAML / HyPhy 选择压力分析或密码子分区建树
phylotrim -i aligned_cds/ -o trimmed_cds/ -g 0.3 --codon
```

---

## 5. 常见问题

**Q：为什么核酸比对中 `N` 和 `R/Y/S/W/K/M` 等简并碱基也会算作 gap？**  
A：在缺失位点过滤中，不确定碱基与简并符号代表未解析状态。只将标准 `A/T/C/G/U` 计为有效位点可保证严格且可复现的列完整度定义。

**Q：如果 `-i` 和 `-o` 指向同一个路径会怎样？**  
A：程序会在运行前检测规范路径并直接报错退出，防止误覆盖原始比对文件。

**Q：批量处理中途断电或中断会产生半截损坏文件吗？**  
A：不会。每个输出文件都会先写入目标目录下的隐形临时文件，写完后再通过原子 `rename` 替换到位。
