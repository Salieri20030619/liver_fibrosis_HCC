function (pq_tree, curr_node, path_lengths) 
{
    if (V(pq_tree)[curr_node]$type != "Q") {
        path_lengths[curr_node] = 0
        return(path_lengths)
    }
    else {
        children_count = 0
        for (child in V(pq_tree)[suppressWarnings(.nei(curr_node, 
            mode = "out"))]) {
            children_counts <- count_leaf_descendents(pq_tree, 
                child, children_counts)
            if (V(pq_tree)[child]$type == "leaf") {
                children_count <- children_count + 1
            }
            else {
                children_count <- children_count + children_counts[child]
            }
        }
        path_lengths[curr_node] = children_count
        return(children_counts)
    }
}
