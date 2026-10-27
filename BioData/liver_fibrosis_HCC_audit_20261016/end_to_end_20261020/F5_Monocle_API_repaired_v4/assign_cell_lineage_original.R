function (pq_tree, curr_node, assigned_state, node_states) 
{
    if (V(pq_tree)[curr_node]$type == "leaf") {
        node_states[V(pq_tree)[curr_node]$name] = assigned_state
        return(node_states)
    }
    else {
        for (child in V(pq_tree)[suppressWarnings(nei(curr_node, 
            mode = "out"))]) {
            node_states <- assign_cell_lineage(pq_tree, child, 
                assigned_state, node_states)
        }
        return(node_states)
    }
}
