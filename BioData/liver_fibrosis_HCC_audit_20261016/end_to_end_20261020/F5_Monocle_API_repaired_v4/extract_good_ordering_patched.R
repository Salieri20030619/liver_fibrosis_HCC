function (pq_tree, curr_node, dist_matrix) 
{
    if (V(pq_tree)[curr_node]$type == "leaf") {
        return(V(pq_tree)[curr_node]$name)
    }
    else if (V(pq_tree)[curr_node]$type == "P") {
        p_level <- list()
        for (child in V(pq_tree)[suppressWarnings(.nei(curr_node, 
            mode = "out"))]) {
            p_level[[length(p_level) + 1]] <- extract_good_ordering(pq_tree, 
                child, dist_matrix)
        }
        p_level <- order_p_node(p_level, dist_matrix)
        p_level <- unlist(p_level)
        return(p_level)
    }
    else if (V(pq_tree)[curr_node]$type == "Q") {
        q_level <- list()
        for (child in V(pq_tree)[suppressWarnings(.nei(curr_node, 
            mode = "out"))]) {
            q_level[[length(q_level) + 1]] <- extract_good_ordering(pq_tree, 
                child, dist_matrix)
        }
        q_level <- order_q_node(q_level, dist_matrix)
        q_level <- q_level$ql
        q_level <- unlist(q_level)
        return(q_level)
    }
}
