![scTransferLearning](https://raw.githubusercontent.com/dosorio/scTransferLearning/main/Figures/F1.png?token=GHSAT0AAAAAABQZDAMLAJDQ3VZRRSKWIA5EYSXCFLQ)
# Interested in mapping your cells?
Please install [symphony](https://github.com/immunogenomics/symphony) and download the following files: `Code > umapBRCA` and `Results > refBRCA.RData`, then you can use the following commands to map your cells to the reference:

```R
library(symphony)
load('refBRCA.RData')
Q <- mapQuery(yourData, metadata_query = yourMetaData, ref_obj = clBRCA)
Q <- knnPredict(Q, clBRCA, clBRCA$meta_data$cellLine, k = 5)
```

The object Q will contain the learned cell line for each cell used as input.
