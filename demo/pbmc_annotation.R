# 使用外部 GPT、DeepSeek 或 Kimi 注释 Seurat PBMC 3K 数据。
#
# 数据下载链接：
# https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
#
# 请在 DeepCellSeek.Rproj 中运行本脚本，数据会下载到 demo/inputs/。
#
# Kimi：Sys.setenv(KIMI_API_KEY = "你的 Kimi API Key")
# DeepSeek：Sys.setenv(DEEPSEEK_API_KEY = "你的 DeepSeek API Key")
# 外部 GPT 中转站（OpenAI Responses API）：
# 与 Codex auth.json 相同：Sys.setenv(OPENAI_API_KEY = "你的 API Key")
# 也支持：Sys.setenv(DEEPCELLSEEK_EXTERNAL_API_KEY = "你的 API Key")
# Sys.setenv(DEEPCELLSEEK_EXTERNAL_BASE_URL = "https://sub2.hongliantina.xyz")
# 若中转站要求 /v1/responses，可设置：
# Sys.setenv(DEEPCELLSEEK_EXTERNAL_ENDPOINT_PATH = "/v1/responses")
# 可选：Sys.setenv(DEEPCELLSEEK_REASONING_EFFORT = "xhigh")
# 在下方的 model 变量中选择使用外部 GPT、Kimi 或 DeepSeek。

if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("请先安装 Seurat：install.packages('Seurat')")
}

library(DeepCellSeek)

project_root <- getwd()
input_dir <- file.path(project_root, "demo", "inputs")
output_dir <- file.path(project_root, "demo", "outputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 使用 Seurat 官方 PBMC 3K 教程中的 10x Genomics 计数矩阵。
dataset_url <- paste0(
  "https://cf.10xgenomics.com/samples/cell/pbmc3k/",
  "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)
archive_file <- file.path(input_dir, "pbmc3k_filtered_gene_bc_matrices.tar.gz")
matrix_dir <- file.path(input_dir, "filtered_gene_bc_matrices", "hg19")
markers_file <- file.path(input_dir, "pbmc3k_official_markers.rds")
seurat_file <- file.path(input_dir, "pbmc3k_official_seurat.rds")
allowed_cell_types_file <- file.path(input_dir, "PeripheralBlood_celltype.rds")

# 第一次运行下载并解压；后续运行直接复用下载的数据。
if (!dir.exists(matrix_dir)) {
  if (!file.exists(archive_file)) {
    message("正在下载 PBMC 3K 数据：", dataset_url)
    utils::download.file(dataset_url, archive_file, mode = "wb", quiet = FALSE)
  }
  message("正在解压数据到：", input_dir)
  utils::untar(archive_file, exdir = input_dir)
}

if (!dir.exists(matrix_dir)) {
  stop("PBMC 3K 数据解压后未找到计数矩阵：", matrix_dir)
}

# 严格按 Seurat PBMC 3K 教程预处理。该教程的数据只有一个样本，不包含批次校正。
if (file.exists(seurat_file)) {
  message("读取已缓存的 Seurat 对象：", seurat_file)
  pbmc <- readRDS(seurat_file)
} else {
  message("正在运行 Seurat 官方 PBMC 3K 预处理、聚类和 UMAP 流程。")
  pbmc <- Seurat::CreateSeuratObject(
    counts = Seurat::Read10X(data.dir = matrix_dir),
    project = "pbmc3k",
    min.cells = 3,
    min.features = 200
  )
  pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  pbmc <- Seurat::NormalizeData(pbmc)
  pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- Seurat::ScaleData(pbmc, features = all.genes)
  pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))
  pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
  pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
  pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
  saveRDS(pbmc, seurat_file)
  message("已缓存 Seurat 对象：", seurat_file)
}

# 兼容不完整的官方流程缓存。
if (!"umap" %in% names(pbmc@reductions)) {
  pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
  saveRDS(pbmc, seurat_file)
}

# 缓存 Seurat 的 FindAllMarkers 结果，避免每次调用 API 前重复计算。
if (file.exists(markers_file)) {
  message("读取已缓存的 Seurat marker：", markers_file)
  markers_df <- readRDS(markers_file)
} else {
  markers_df <- Seurat::FindAllMarkers(
    pbmc,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  saveRDS(markers_df, markers_file)
  message("已缓存 Seurat marker：", markers_file)
}

if (!file.exists(allowed_cell_types_file)) {
  stop("未找到允许的细胞类型列表：", allowed_cell_types_file)
}
allowed_cell_types <- readRDS(allowed_cell_types_file)

# 模型选择（三选一；只保留所选模型的 model 赋值）：
# model <- "kimi-k2.6"                # Kimi
# Sys.setenv(KIMI_API_KEY = "你的 API Key")
#
# model <- "deepseek-v4-flash"          # DeepSeek
# Sys.setenv(DEEPSEEK_API_KEY = "")

# 默认使用外部 GPT-5.6-sol；建议使用 xhigh 推理强度以获得更充分的分析。
model <- "gpt-5.6-sol"
Sys.setenv(DEEPCELLSEEK_REASONING_EFFORT = "xhigh")

# 如果不使用 OPENAI_API_KEY 环境变量，可取消下一行注释并填入真实密钥：
# Sys.setenv(OPENAI_API_KEY = "sk-xxx")
Sys.setenv(DEEPCELLSEEK_EXTERNAL_BASE_URL = "https://hk1.r7z.net")
# 若中转站要求 /v1/responses，可设置：
# Sys.setenv(DEEPCELLSEEK_EXTERNAL_ENDPOINT_PATH = "/v1/responses")


api_key_envs <- switch(
  model,
  "kimi-k2.6" = "KIMI_API_KEY",
  "deepseek-v4-flash" = "DEEPSEEK_API_KEY",
  "gpt-5.6-sol" = c("OPENAI_API_KEY", "DEEPCELLSEEK_EXTERNAL_API_KEY")
)
api_key_values <- Sys.getenv(api_key_envs, unset = "")
if (!any(nzchar(api_key_values))) {
  stop("请先设置以下任一变量：", paste(api_key_envs, collapse = " 或 "), "，再运行此 demo。")
}

# 对于基因数量较多或经过多批次整合的数据集，建议每个 cluster 提供前 50 个
# 差异 marker 基因，以保留足够的注释信息；小型或低质量数据可酌情减少该数值。
# allowed_cell_types 限制模型只能使用 RDS 中的候选类型；混合群可返回
# "候选类型 A + 候选类型 B"，其中每个组成部分均须来自候选列表。
annotations <- llm_celltype(
  input = markers_df,
  tissuename = "PBMC",
  species = "Human",
  model = model,
  topgenenumber = 50,
  wait_indefinitely = TRUE,
  allowed_cell_types = allowed_cell_types
)

print(data.frame(
  cluster = names(annotations),
  annotation = unname(annotations),
  row.names = NULL
))

# 始终使用不变的 seurat_clusters，支持在同一会话重复比较不同模型。
cluster_ids <- as.character(pbmc$seurat_clusters)
pbmc$LLM_Annotation <- unname(annotations[cluster_ids])
seurat_reference <- c(
  "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
  "FCGR3A+ Mono", "NK", "DC", "Platelet"
)
names(seurat_reference) <- as.character(0:8)
pbmc$Seurat_Reference <- unname(seurat_reference[cluster_ids])

if (anyNA(pbmc$LLM_Annotation)) {
  stop("模型注释未覆盖全部 Seurat cluster，无法绘制对照图。")
}

# 左图为 Seurat 教程参考注释，右图为本次模型的注释。
model_label <- switch(
  model,
  "kimi-k2.6" = "Kimi",
  "deepseek-v4-flash" = "DeepSeek",
  "gpt-5.6-sol" = "External GPT",
  model
)
seurat_umap <- Seurat::DimPlot(
  pbmc,
  reduction = "umap",
  group.by = "Seurat_Reference",
  label = TRUE,
  pt.size = 0.5
) + Seurat::NoLegend() + ggplot2::ggtitle("Seurat PBMC 3K reference")
llm_umap <- Seurat::DimPlot(
  pbmc,
  reduction = "umap",
  group.by = "LLM_Annotation",
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  pt.size = 0.5
) +
  Seurat::NoLegend() +
  ggplot2::ggtitle(paste0(model_label, " annotation (", model, ")")) +
  ggplot2::theme(plot.margin = ggplot2::margin(12, 28, 12, 28))
umap_comparison <- seurat_umap + llm_umap
output_file <- file.path(output_dir, paste0("pbmc3k_umap_", model, ".pdf"))
ggplot2::ggsave(
  filename = output_file,
  plot = umap_comparison,
  width = 14,
  height = 7,
  units = "in"
)
message("UMAP 对比图已保存到：", output_file)
