# FCA_Fetal_Skin_priv
Private repo for the preparation of the fetal skin manuscript

See [description of contents](toc.md)

## Data Science Pipeline

Click on the name of the file of interest to access it. 

```mermaid
flowchart LR;
    classDef green color:#022e1f,fill:#00f500;
    classDef red color:#022e1f,fill:#DF4F2F;
    classDef blue color:#022e1f,fill:#6EC7F7;
    classDef white color:#022e1f,fill:#fff;
    classDef black color:#fff,fill:#000;
    subgraph Raw Data pre-Processing
    A[(Raw *raw.h5 Database)]:::green-->B>0_collectRaw.ipynb]:::red;
    AA[(*cell_source.csv Database)]:::blue-.->B;
    click B "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/0_collectRaw.ipynb"
    AAA[(fetal_skin_samples.txt)]:::blue-.->B;
    AAAA[(fetal_annotation_rachel-201904.txt)]:::blue-.->B;
    click AAA "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/data/fetal_skin_samples.txt"
    click AAAA "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/data/fetal_annotation_rachel-201904.txt"
    B-->C[(fetal_skin_raw.20190926.h5ad)]:::green;
    BB[(hierarchy1.txt)]:::blue-.->C
    click BB "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/data/hierarchy1.txt"
    end
    subgraph Preliminary Data Analysis
    C-->D>1_broadCategory.ipynb]:::red;
    CC[(hierarchy1.txt)]:::blue;
    CC-.->D;
    click CC "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/data/hierarchy1.txt"
    click D "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/1_broadCategory.ipynb"
    D-->E[(fetal_skin_processed.hierarchy1.h5ad)]:::green;
    D-->F[(fetal_skin_bbknn.h5ad)]:::green;
    D-->G[(fetal_skin.hierarchy1_stroma.h5ad)]:::green;
    D-->H[(fetal_skin.hierarchy1_myloid.h5ad)]:::green;
    D-->I[(fetal_skin.hierarchy1_T-cells.h5ad)]:::green;
    D-->J[(fetal_skin.hierarchy1_B-cells.h5ad)]:::green;
    D-->K[(fetal_skin.hierarchy1_endothelium.h5ad)]:::green;
    D-->L[(fetal_skin.hierarchy1_mast-cells.h5ad)]:::green;
    D-->M[(fetal_skin.hierarchy1_keratinocytes.h5ad)]:::green;
    D-->N[(fetal_skin.hierarchy1_melanocytes.h5ad)]:::green;
    D-->O[(fetal_skin.hierarchy1_erythroid.h5ad)]:::green;
    end
    subgraph Cell-Type-Specific Data Analysis
    G & C-->P>2.1_Stroma_v2.ipynb]:::red;
    H & C-->Q>2.2_Myloid_cells_v2.ipynb]:::red;
    II[(G1-S_genes.list)]:::blue;
    III[(G2-M_genes.list)]:::blue;
    IIII[(JP_cycle_genes.list)]:::blue;
    I & C-->R>2.3_T_cells_v2.ipynb]:::red;
    AAAA & BB & II & III & IIII-.->R & S & V & W
    click II "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/data/G1-S_genes.list"
    click III "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/data/G2-M_genes.list"
    click IIII "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/data/JP_cycle_genes.list"
    J & C-->S>2.4_B_cells_v2.ipynb]:::red;
    K & C-->T>2.5_Endothelium_v2.ipynb]:::red;
    L & C-->U>2.6_Mast_cells_v2.ipynb]:::red;
    M & C-->V>2.7_Keratinocytes_v2.ipynb]:::red;
    N & C-->W>2.8_Melanocytes_v2.ipynb]:::red;
    O & C-->X>2.9_Erythroid_v2.ipynb]:::red;
    click P "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/2.1_Stroma_v2.ipynb"
    click Q "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/2.2_Myloid_cells_v2.ipynb"
    click R "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/2.3_T_cells_v2.ipynb"
    click S "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/2.4_B_cells_v2.ipynb"
    click T "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/2.5_Endothelium_v2.ipynb"
    click U "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/2.6_Mast_cells_v2.ipynb"
    click V "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/2.7_Keratinocytes_v2.ipynb"
    click W "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/2.8_Melanocytes_v2.ipynb"
    click X "https://github.com/haniffalab/FCA_Fetal_Skin_priv/blob/main/code/notebooks/analysis/2.9_Erythroid_v2.ipynb"
    P-->PP[(fetal_skin_hierarch1_stroma_processed.h5ad)]:::green;
    Q-->QQ[(fetal_skin_hierarch1_myloid_processed.h5ad)]:::green;
    R-->RR[(fetal_skin_hierarch1_T-cells_processed.h5ad)]:::green;
    S-->SS[(fetal_skin_hierarch1_B-cells_processed.h5ad)]:::green;
    T-->TT[(fetal_skin_hierarch1_endothelium_processed.h5ad)]:::green;
    U-->UU[(fetal_skin_hierarch1_mast-cells_processed.h5ad)]:::green;
    V-->VV[(fetal_skin_hierarch1_keratinocytes_processed.h5ad)]:::green;
    W-->WW[(fetal_skin_hierarch1_melanocytes_processed.h5ad)]:::green;
    X-->XX[(fetal_skin_hierarch1_erythroid_processed.h5ad)]:::green;
    PP-->P
    QQ-->Q
    RR-->R
    SS-->S
    TT-->T
    UU-->U
    VV-->V
    WW-->W
    XX-->X
    end
```
