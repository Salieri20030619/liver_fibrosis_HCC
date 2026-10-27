function (cell_name, cth, curr_node, gate_res) 
{
    next_nodes <- c()
    for (child in V(cth@classificationTree)[suppressWarnings(.nei(curr_node, 
        mode = "out"))]) {
        type_res <- gate_res[[V(cth@classificationTree)[child]$name]]
        if (type_res[cell_name] == TRUE) 
            next_nodes <- c(next_nodes, V(cth@classificationTree)[child]$name)
    }
    if (length(next_nodes) == 1) {
        CellType <- cth_classifier_cell(cell_name, cth, next_nodes[1], 
            gate_res)
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
