function (pq_tree) 
{
    type <- NA
    canonical_pq <- pq_tree
    V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, 
        mode = "out") == 2]$color = "black"
    V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, 
        mode = "out") == 2]$type = "Q"
    single_child_p <- V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, 
        mode = "out") == 1]
    V(canonical_pq)[type == "P" & igraph::degree(canonical_pq, 
        mode = "out") == 1]$color = "blue"
    for (p_node in single_child_p) {
        child_of_p_node <- V(canonical_pq)[suppressWarnings(nei(p_node, 
            mode = "out"))]
        parent_of_p_node <- V(canonical_pq)[suppressWarnings(nei(p_node, 
            mode = "in"))]
        for (child_of_p in child_of_p_node) {
            canonical_pq[parent_of_p_node, child_of_p] <- TRUE
        }
    }
    canonical_pq <- delete.vertices(canonical_pq, V(canonical_pq)[type == 
        "P" & igraph::degree(canonical_pq, mode = "out") == 1])
    return(canonical_pq)
}
