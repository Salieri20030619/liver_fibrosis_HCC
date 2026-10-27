function (cds_subset, cth, curr_node, frequency_thresh) 
{
    next_nodes <- c()
    for (child in V(cth@classificationTree)[suppressWarnings(.nei(curr_node, 
        mode = "out"))]) {
        child_cell_class_func <- V(cth@classificationTree)[child]$classify_func[[1]]
        type_res <- child_cell_class_func(exprs(cds_subset))
        type_res <- unlist(type_res)
        names(type_res) <- row.names(pData(cds_subset))
        cell_type_name <- V(cth@classificationTree)[child]$name
        if (length(frequency_thresh) > 1) 
            required_thresh <- frequency_thresh[cell_type_name]
        else required_thresh <- frequency_thresh
        if ((sum(type_res)/length(type_res)) > frequency_thresh) {
            next_nodes <- c(next_nodes, cell_type_name)
        }
    }
    if (length(next_nodes) == 1) {
        CellType <- cth_classifier_cds(cds_subset, cth, next_nodes[1], 
            frequency_thresh)
    }
    else if (length(next_nodes) == 0) {
        if (curr_node == "root") 
            CellType = "Unknown"
        else CellType = curr_node
    }
    else if (length(next_nodes) > 1) {
        CellType = "Ambiguous"
    }
    else {
        CellType = "Unknown"
    }
    return(CellType)
}
