# PhyloTrim

**极简、极速的多序列比对（MSA）列缺口与密码子修剪工具（Rust 编写）**

`phylotrim` 只专注做好两件事，并用底层硬件并行把速度推到极致（无论是一次修剪成千上万个基因文件，还是修剪单个几百万至上亿位点的超大型比对矩阵）：

1. **按列删除高缺失列（默认）**：只要某一列中**非标准字符（gap）**的比例超过设定阈值 `-g`，就删除该列。
2. **按密码子（3 bp）整组删除（`-c / --codon`）**：每 3 列作为一个密码子，只要其中**任意 1 列**超过阈值 `-g`，就将整个密码子（3 列）一起删除，严格保持阅读框。

> **绝不删行**：`phylotrim` 仅修剪列，永远完整保留输入中的所有物种行及其原始顺序。

---

## 字符判定规则（自动识别核酸与蛋白质）

程序自动识别输入序列类型，无需手动指定：

- **核酸（DNA / RNA）**：仅 `A, T, C, G, U`（大小写不敏感）计为有效标准碱基；其余所有字符（`-`, `.`, `?`, `N`, 以及 `R/Y/S/W/K/M` 等简并符号）**一律计为 gap**。
- **蛋白质（Protein）**：仅 **20 种标准氨基酸**（`A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y`，大小写不敏感）计为有效字符；其余所有字符（`-`, `.`, `?`, `X`, `*`, `B`, `Z`, `J`, `U`, `O` 等）**一律计为 gap**。

---

## 编译安装

```bash
git clone https://github.com/Jhe1004/phylotrim.git
cd phylotrim
cargo build --release
# 可执行文件位于 target/release/phylotrim
```

也可以直接安装到 `$HOME/.cargo/bin`：

```bash
cargo install --path .
```

---

## 命令行参数（仅 5 个选项）

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

## 常用示例

```bash
# 1. 单个基因比对：删除非 ATCG 比例 > 20% 的列
phylotrim -i aln.fasta -o aln_trimmed.fasta -g 0.2

# 2. 密码子模式：3 列中只要有 1 列非 ATCG 比例 > 30%，整组 3 列一起删掉
phylotrim -i cds.fasta -o cds_trimmed.fasta -g 0.3 --codon

# 3. 批量并行处理整个目录下的成千上万个基因（自动吃满所有 CPU 核心，兼容 .gz）
phylotrim -i raw_alignments/ -o trimmed_alignments/ -g 0.2

# 4. 指定 16 线程处理单个超大型染色体/全基因组比对矩阵
phylotrim -i chr1_huge.fas.gz -o chr1_trimmed.fas.gz -g 0.1 -t 16
```

---

## 加速设计概览

- **海量基因文件级工作窃取并行**：目录模式下按文件体积降序（LPT）调度至 Rayon 线程池，消除大基因拖尾等待。
- **超大基因内部 L1-Cache 分块 + SIMD 向量化扫描**：将超长序列切分为驻留 L1 缓存的 `8,192` 列块（Tiles）与 `255` 行无溢出 `u8` 微批累加器，触发编译器生成 ARM NEON / x86 AVX2 SIMD 向量指令，单大基因自动多核并行扫描。
- **零拷贝快径与连续区间 `memcpy` 重构**：平铺式零二次分配 FASTA 解析；全部列达标时零拷贝直出；部分列删除时自动合并连续保留列区间（Run-Length Coalescing）并行整块 `copy_from_slice`。

## License

MIT
