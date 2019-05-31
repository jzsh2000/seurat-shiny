# 简介

这是一个 [rstudio shiny](https://shiny.rstudio.com/) 应用，借助
[Seurat](http://satijalab.org/seurat/) 包对 single-cell RNA-seq 数据进行可视化
分析。本应用最初为 [10X](https://www.10xgenomics.com/) 数据设计，从原理上来说
其它平台的 single-cell RNA-seq 数据经 Seurat 处理后也能适用，但尚未测试。

本应用目前实现了以下功能：

* 对于一套 single-cell RNA-seq 数据，显示 t-SNE 投影图以及不同参数下的聚类结果。
* 根据输入的基因名称，在 t-SNE 投影图中以颜色深浅显示图中不同区域的细胞的基因表
  达量的高低，同时使用 violin plot 表示不同小群的基因表达情况。
* 根据输入的两个基因，显示这两个基因在某一个或某几个小群中的表达相关关系，并提
  供 pearson 相关系数。
* 通过选中某一个小群，找出这一小群细胞相比于其它某个小群或者剩余所有细胞的特征
  性高表达基因。
* 对于选定数据集的聚类结果，选择小群重新进行聚类和 t-SNE 分析。

# 环境准备

赖于的 R 包显示如下：(除了`limma`需要使用[Bioconductor](http://bioconductor.org/)
的`biocLite()`命令安装之外，其它包均可通过`install.packages()`命令安装):

* [shiny](https://github.com/rstudio/shiny)
* [shinyjs](https://github.com/daattali/shinyjs)
* [tidyverse](https://github.com/tidyverse/tidyverse)
* [magrittr](https://github.com/tidyverse/magrittr)
* [DT](https://github.com/rstudio/DT)
* [cowplot](https://github.com/wilkelab/cowplot)
* [Seurat](https://github.com/satijalab/seurat)
* [limma](https://github.com/cran/limma)

# 应用部署

1. 更改配置文件

    进入 `shiny/` 目录，将 `config.txt.example` 复制为 `config.txt`，编辑
    `config.txt` 以修改应用的全局变量

2. 建立数据列表

    进入 `shiny/data/` 目录，将 `resource_list.json.example` 复制为
    `resource_list.json`。这一 JSON 文件记录了应用关联的数据的详细信息。

    以`resource_list.json.example` 文件中的内容说明各键值对的含义：
    ```json
    [
      {
        "label": "bl1",
        "description": "[BL1] mouse blood",
        "species": "mouse",
        "dims.use": 12,
        "resolution": 0.3,
        "subsets": [
          {
            "label": "bl1_dc",
            "description": "[BL1_DC] mouse blood DC"
          }
        ]
      },
      {
        "label": "sp2",
        "description": "[SP2] mouse spleen,mouse",
        "species": "mouse",
        "dims.use": 12,
        "resolution": 0.6
      }
    ]
    ```

    * __label__ - 样本简称，这一项不应该包含空格等空白字符。同时这也是样本数据所在的目录（相对于 `shiny/data/`  文件夹），如果该项的值为`foobar`，则样本数据位于`shiny/data/foobar`文件夹下
    * __description__ - 样本描述，即关于该样品的详细信息
    * __species__ - 样本来源的物种名，对于人和小鼠的样本应该分别为`human`和`mouse`
    * __dims.use__ - 样本处理过程中 `Seurat` 包中 `RunTSNE()` 和 `RunUMAP()` 使用的 `dims` 参数，表示主成分（PC）数量，若该项未设置则使用 `shiny/config.txt` 中 `dims_default` 的值。
    * __resolution__ - 样本处理过程中 `Seurat` 包中 `FindClusters()` 使用的`resolution` 参数，若该项未设置则使用 `shiny/config.txt` 中 `res_default` 的值。设计这一项的意义在于：如果曾经使用不同的`resolution` 参数多次执行 `FindClusters()` 函数，这一项给出了数据提供者认为的最适合的参数值。
    * __subsets__ - 预先计算的样本的子集，其下同样可以设定 `label`、`description` 等键。

3. 数据准备

    对于每一个样本（即上述 `resource_list.json` 文件中的每一条记录），需要在相应
    的数据文件夹中放置两类数据：

    1. _XXX_.rds （这里的 _XXX_ 应该与前述的 __label__ 保持一致，为使用 `Seurat`
       处理得到的 R 对象文件，使用 `saveRDS()` 或者 `readr::write_rds()` 函数保
       存到本地）
    2. 0至若干个 _YYY_.rds，表示数据的子集处理得到的 R 对象文件。

    为了得到 _XXX_.rds 或 _YYY_.rds 文件，需要用 `Seurat` 对 single-cell RNA-seq 数据进行一系列处理，包括数据过滤，降维处理，聚类分析等，实际的处理流程可以参考[analysis.Rmd](shiny/example/analysis.Rmd) 文件。


# 其它说明

* Cell Ranger 输出的 t-SNE 坐标可以通过 `Seurat` 包中的 `SetDimReduction()`
  函数整合到 R 对象文件中，参考脚本 `utils/update-seurat-obj.R`。Cell Ranger
  输出的 t-SNE 坐标文件为 projection.csv。对于 10X 数据，该文件可以在 10X pipeline 输出结果的路径 `./outs/analysis/tsne/2_components/projection.csv`
  找到)

* 对基因别名的支持：对于人(human) 和小鼠(mouse) 的数据，基因名的输入框可以使用
  基因别名。例如输入的 CD56 会被自动转换为基因的标准名 NCAM1。若输入的基因名有
  效，则基因名输入框的边框会变为绿色，否则会变为红色。

* 可以直接在网址中指定使用的数据集：若应用部署在 localhost:4533 这个页面，则在
  浏览器地址栏中输入 localhost:4533/?dataset=bl1 将直接打开名为 `bl1` 的数据集
  （这里的 `bl1` 在前述的 `resource_list.csv` 文件中体现为一条记录的 __label__
  列）
