function (pq_tree, curr_node, children_counts) 
{
    if (V(pq_tree)[curr_node]$type == "leaf") {
        children_counts[curr_node] = 0
        return(children_counts)
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
        children_counts[curr_node] = children_count
        return(children_counts)
    }
}
