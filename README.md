# ScType operator

##### Description

The `ScType operator` : Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic or flow cytometry data.

##### Usage

Input projection|.
---|---
`y-axis`        | numeric, value (count or MFI)
`row`           | represents the variables (e.g. genes, channels, markers)
`column`        | represents the observations (e.g. cells, samples, individuals)
`colors`        | Cluster id 
`labels`        | documentID of the file containing the cell_type-marker table.

Input parameters|.
---|---
`RNA_check_name`        | boolean, check the genes name for RNA seq data 
`tissue`        | character, select the tissue: Immune system(default), Liver, Pancreas, Kidney, Eye, Brain or other
`confidence threshold`        | numeric, select the confidence threshold to be used

Output relations|.
---|---
`Population`        | list of cellular population
`solo score`        | return the score for each observation
`cluster score`        | return the score for each of the cluster id
`cluster`        | return the cluster id set as input in colors
##### Details

ScType is a computational method for automated selection of marker based merely on scRNA-seq data.
Custom database can be inputed (as labels), those manual DB file should contain four columns:

| tissueType      | cellName | geneSymbolmore1   | geneSymbolmore2   |
| :---        |    :----:   |    :----:   |          ---: | 
| Immune system      | Naive B cells       | CD19,IgD,CD38,CD24,... |          |
| Immune system   | Natural killer  cells        | CD56,CD2,CD16,CD94,...      |      CD3,CD4,CD8     |

with geneSymbolmore1 as positive markers and geneSymbolmore2 as markers not expected to be expressed by a cell type.

Low-confident scores (sctype cluster score less than the number of observation of this cluster divided by the confidence threshold) are set to "unknown".

##### See Also

[Sctype git](https://github.com/IanevskiAleksandr/sc-type)

