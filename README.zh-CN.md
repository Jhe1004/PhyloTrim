<div align="center">

# PhyloTrim

**极简、极速的多序列比对（MSA）列缺口与 3 bp 密码子修剪工具（Rust 编写）**

[![CI](https://github.com/Jhe1004/PhyloTrim/actions/workflows/ci.yml/badge.svg)](https://github.com/Jhe1004/PhyloTrim/actions/workflows/ci.yml)
[![Version](https://img.shields.io/badge/version-0.2.0-blue.svg)](CHANGELOG.md)
[![Rust](https://img.shields.io/badge/rust-2021%20edition-orange.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#编译安装)

[English](README.md) | **简体中文**

[核心特性](#为什么选择-phylotrim) • [编译安装](#编译安装) • [快速上手](#快速上手) • [参数说明](#命令行参数) • [完整文档](#文档导航) • [性能实测](#性能快照)

</div>

---

## 为什么选择 PhyloTrim？

在系统发育基因组学（Phylogenomics）建树流程中，经常需要一次修剪成千上万个直系同源基因文件，或修剪长达数百万至上亿位点的超大型拼接矩阵。**PhyloTrim**（`phylotrim`）摒弃冗余启发式参数，只专注用**确定性语义与底层硬件并行**做好两件事：

| 核心能力 | 设计保证 |
| :--- | :--- |
| **按列删除高缺失列（默认）** | 只要某一列中非标准字符/缺口比例严格超过 `-g`（`gap_count / n_taxa > max_gap`），即删除该列。 |
| **密码子防移码修剪（`-c`）** | 按连续 3 bp 密码子整组判定：只要 3 列中**任意 1 列**超过阈值 `-g`，整组 3 列一起删除，严格维持阅读框。 |
| **绝不删行** | 仅对列进行修剪，100% 保留输入中的所有物种行及其原始排列顺序。 |
| **自动识别核酸与蛋白质** | 无需手动指定序列类型，自动区分 DNA/RNA 与 Protein，并应用精确的标准字母表判定规则。 |
| **零依赖原生极致速度** | 单一静态 Rust 二进制文件，结合 LPT 文件级调度、L1 缓存分块 SIMD 向量化计数与连续区间零拷贝重构。 |

---

## 编译安装

需要安装 [Rust 工具链](https://rustup.rs/)（`cargo >= 1.70`）：

```bash
git clone https://github.com/Jhe1004/PhyloTrim.git
cd PhyloTrim
cargo build --release
# 生成的可执行文件位于 target/release/phylotrim
```

也可以直接安装到 `$HOME/.cargo/bin`：

```bash
cargo install --git https://github.com/Jhe1004/PhyloTrim.git
```

验证安装：

```bash
phylotrim --version
```

---

## 快速上手

```bash
# 1. 单个基因比对：删除非标准碱基 / gap 比例 > 20% 的列
phylotrim -i alignment.fasta -o alignment.trimmed.fasta -g 0.2

# 2. 密码子模式：3 列中只要有 1 列超 30% 缺失，整组 3 列一起删掉
phylotrim -i cds.fasta -o cds.trimmed.fasta -g 0.3 --codon

# 3. 批量并发处理整个目录的成千上万个基因文件（自动吃满全部 CPU 核心，原生兼容 .gz）
phylotrim -i raw_loci/ -o trimmed_loci/ -g 0.2

# 4. 指定 16 线程处理单个超大型染色体/超矩阵 .gz 文件
phylotrim -i chr1_supermatrix.fas.gz -o chr1_trimmed.fas.gz -g 0.1 -t 16
```

---

## 命令行参数

`phylotrim` 精简至仅 **5 个选项**：

```text
Usage: phylotrim [OPTIONS] --input <PATH> --output <PATH>

Options:
  -i, --input <PATH>     输入比对文件或目录（自动识别 FASTA / PHYLIP / NEXUS，原生支持 .gz）
  -o, --output <PATH>    输出比对文件或目录（支持 .gz）
  -g, --max-gap <FLOAT>  单列允许的最大非标准字符/gap 比例 [0.0..=1.0] [默认: 0.2]
  -c, --codon            按 3 bp 密码子模式修剪（3 列中任一列超阈值则整组删除）
  -t, --threads <INT>    并发线程数（0 = 自动吃满全部 CPU 核心）[默认: 0]
  -h, --help             查看帮助
  -V, --version          查看版本
```

---

## 字符判定规则

程序自动识别输入比对的序列类型：

| 序列类型 | 有效标准字符（大小写不敏感） | 计为 gap / 非标准字符 |
| :--- | :--- | :--- |
| **核酸（DNA / RNA）** | `A, T, C, G, U` | `-`、`.`、`?`、`N` 及全部 IUPAC 简并碱基（`R, Y, S, W, K, M, B, D, H, V` 等） |
| **蛋白质（Protein）** | 20 种标准氨基酸（`A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y`） | `-`、`.`、`?`、`X`、`*`、`B`、`Z`、`J`、`U`、`O` 等 |

> **安全防护**：对蛋白质比对误开启 `--codon` 时会直接报错拦截，避免静默破坏氨基酸列。

---

## 性能快照

在本地 Apple Silicon 工作站上的端到端实测（含完整磁盘读写、解析、修剪与原子落盘，`phylotrim v0.2.0`，默认 `-t 0`）：

| 测试场景 | 数据规模 | 模式 | 端到端耗时 |
| :--- | :--- | :--- | :---: |
| 单个超大比对矩阵（38.1 MB FASTA） | 200 物种 × 200,000 位点 | 按列（`-g 0.2`） | **0.030 s** |
| 单个超大比对矩阵（38.1 MB FASTA） | 200 物种 × 200,000 位点 | 密码子（`-g 0.2 --codon`） | **0.025 s** |
| 多基因目录批处理（总计 ~75 MB） | 1,000 个文件 ×（50 物种 × 1,500 bp） | 按列（`-g 0.2`） | **0.110 s** |

关于底层 L1 缓存分块、无分支 SIMD 累加器与 LPT 调度的详细实现，请参阅 [Algorithm & Performance Architecture](docs/ALGORITHM_AND_PERFORMANCE.md)。

---

## 文档导航

- **[用户手册（简体中文）](docs/zh-CN/USER_GUIDE.md)** — 支持格式详情、单文件与目录批处理行为、密码子边界说明、建树流程衔接示例与常见问题。
- **[User Guide (English)](docs/USER_GUIDE.md)** — Complete English user manual.
- **[Algorithm & Performance Architecture](docs/ALGORITHM_AND_PERFORMANCE.md)** — 底层内存布局、SIMD 向量化、L1 缓存分块与原子 I/O 架构详解。
- **[更新日志 (Changelog)](CHANGELOG.md)** — 版本变更记录。

---

## 开源协议

PhyloTrim 采用 [MIT License](LICENSE) 开源。
