# ScType operator

##### Description

The `ScType operator` : Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data.

##### Usage

Input projection|.
---|---
`y-axis`        | numeric, value (count or MFI)
`row`           | type, genes or channels
`column`        | type, row_Id 
`colors`        | type, Cluster id 
`labels`        | type, documentID of the file containing the cell_type-marker table.

Input parameters|.
---|---
`RNA_check_name`        | Check the genes name for RNA seq data
`tissue`        | Select the tissue: Immune system(default), Liver, Pancreas, Kidney, Eye, Brain or other
`confidence threshold`        | Select the confidence threshold to be used

Output relations|.
---|---
`Population`        | output relation
`solo score`        | view of the Shiny application
`cluster score`        | view of the Shiny application

##### Details

ScType a computational method for automated selection of marker genes based merely on scRNA-seq data.

##### See Also

[Sctype git](https://github.com/IanevskiAleksandr/sc-type)

