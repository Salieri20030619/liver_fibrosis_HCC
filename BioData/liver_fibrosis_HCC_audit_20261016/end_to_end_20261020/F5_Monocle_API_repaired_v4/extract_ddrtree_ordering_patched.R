function (cds, root_cell, verbose = T) 
{
    dp <- cellPairwiseDistances(cds)
    dp_mst <- minSpanningTree(cds)
    curr_state <- 1
    res <- list(subtree = dp_mst, root = root_cell)
    states = rep(1, ncol(dp))
    names(states) <- V(dp_mst)$name
    pseudotimes = rep(0, ncol(dp))
    names(pseudotimes) <- V(dp_mst)$name
    parents = rep(NA, ncol(dp))
    names(parents) <- V(dp_mst)$name
    mst_traversal <- igraph::dfs(dp_mst, root = root_cell, mode = "all", 
        unreachable = FALSE, parent = TRUE)
    mst_traversal$parent <- as.numeric(mst_traversal$parent)
    curr_state <- 1
    for (i in 1:length(mst_traversal$order)) {
        curr_node = mst_traversal$order[i]
        curr_node_name = V(dp_mst)[curr_node]$name
        if (is.na(mst_traversal$parent[curr_node]) == FALSE) {
            parent_node = mst_traversal$parent[curr_node]
            parent_node_name = V(dp_mst)[parent_node]$name
            parent_node_pseudotime = pseudotimes[parent_node_name]
            parent_node_state = states[parent_node_name]
            curr_node_pseudotime = parent_node_pseudotime + dp[curr_node_name, 
                parent_node_name]
            if (degree(dp_mst, v = parent_node_name) > 2) {
                curr_state <- curr_state + 1
            }
        }
        else {
            parent_node = NA
            parent_node_name = NA
            curr_node_pseudotime = 0
        }
        curr_node_state = curr_state
        pseudotimes[curr_node_name] <- curr_node_pseudotime
        states[curr_node_name] <- curr_node_state
        parents[curr_node_name] <- parent_node_name
    }
    ordering_df <- data.frame(sample_name = names(states), cell_state = factor(states), 
        pseudo_time = as.vector(pseudotimes), parent = parents)
    row.names(ordering_df) <- ordering_df$sample_name
    return(ordering_df)
}
