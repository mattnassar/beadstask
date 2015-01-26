function entropy=computeEntropy(PDF)
    entropy=nansum(PDF.*log2(1./PDF), 2);