function (orig_pq_tree, curr_node, dist_matrix, num_branches, 
    reverse_main_path = FALSE) 
{
    requireNamespace("plyr")
    nei <- NULL
    type <- NA
    pseudo_time <- NA
    pq_tree <- orig_pq_tree
    branch_node_counts <- V(pq_tree)[type == "Q"]$diam_path_len
    names(branch_node_counts) <- V(pq_tree)[type == "Q"]$name
    if (length(names(branch_node_counts)) < num_branches) 
        stop("Number of branches attempted is larger than the branches constructed from pq_tree algorithm")
    branch_node_counts <- sort(branch_node_counts, decreasing = TRUE)
    cell_states <- rep(NA, length(as.vector(V(pq_tree)[type == 
        "leaf"])))
    names(cell_states) <- V(pq_tree)[type == "leaf"]$name
    cell_states <- assign_cell_lineage(pq_tree, curr_node, 1, 
        cell_states)
    branch_point_roots <- list()
    branch_tree <- graph.empty()
    for (i in 1:num_branches) {
        branch_point_roots[[length(branch_point_roots) + 1]] <- names(branch_node_counts)[i]
        branch_id <- names(branch_node_counts)[i]
        branch_tree <- branch_tree + vertex(branch_id)
        parents <- V(pq_tree)[suppressWarnings(nei(names(branch_node_counts)[i], 
            mode = "in"))]
        if (length(parents) > 0 && parents$type == "P") {
            p_node_parent <- V(pq_tree)[suppressWarnings(nei(names(branch_node_counts)[i], 
                mode = "in"))]
            parent_branch_id <- V(pq_tree)[suppressWarnings(nei(p_node_parent, 
                mode = "in"))]$name
            branch_tree <- branch_tree + edge(parent_branch_id, 
                branch_id)
        }
        pq_tree[V(pq_tree)[suppressWarnings(nei(names(branch_node_counts)[i], 
            mode = "in"))], names(branch_node_counts)[i]] <- FALSE
    }
    branch_pseudotimes <- list()
    for (i in 1:length(branch_point_roots)) {
        branch_ordering <- extract_good_ordering(pq_tree, branch_point_roots[[i]], 
            dist_matrix)
        branch_ordering_time <- weight_of_ordering(branch_ordering, 
            dist_matrix)
        names(branch_ordering_time) <- branch_ordering
        branch_pseudotimes[[length(branch_pseudotimes) + 1]] = branch_ordering_time
        names(branch_pseudotimes)[length(branch_pseudotimes)] = branch_point_roots[[i]]
    }
    cell_ordering_tree <- graph.empty()
    curr_branch <- "Q_1"
    extract_branched_ordering_helper <- function(branch_tree, 
        curr_branch, cell_ordering_tree, branch_pseudotimes, 
        dist_matrix, reverse_ordering = FALSE) {
        nei <- NULL
        curr_branch_pseudotimes <- branch_pseudotimes[[curr_branch]]
        curr_branch_root_cell <- NA
        for (i in 1:length(curr_branch_pseudotimes)) {
            cell_ordering_tree <- cell_ordering_tree + vertex(names(curr_branch_pseudotimes)[i])
            if (i > 1) {
                if (reverse_ordering == FALSE) {
                  cell_ordering_tree <- cell_ordering_tree + 
                    edge(names(curr_branch_pseudotimes)[i - 1], 
                      names(curr_branch_pseudotimes)[i])
                }
                else {
                  cell_ordering_tree <- cell_ordering_tree + 
                    edge(names(curr_branch_pseudotimes)[i], names(curr_branch_pseudotimes)[i - 
                      1])
                }
            }
        }
        if (reverse_ordering == FALSE) {
            curr_branch_root_cell <- names(curr_branch_pseudotimes)[1]
        }
        else {
            curr_branch_root_cell <- names(curr_branch_pseudotimes)[length(curr_branch_pseudotimes)]
        }
        for (child in V(branch_tree)[suppressWarnings(nei(curr_branch, 
            mode = "out"))]) {
            child_cell_ordering_subtree <- graph.empty()
            child_head <- names(branch_pseudotimes[[child]])[1]
            child_tail <- names(branch_pseudotimes[[child]])[length(branch_pseudotimes[[child]])]
            curr_branch_cell_names <- names(branch_pseudotimes[[curr_branch]])
            head_dist_to_curr <- dist_matrix[child_head, curr_branch_cell_names]
            closest_to_head <- names(head_dist_to_curr)[which(head_dist_to_curr == 
                min(head_dist_to_curr))]
            head_dist_to_anchored_branch = NA
            branch_index_for_head <- NA
            head_dist_to_anchored_branch <- dist_matrix[closest_to_head, 
                child_head]
            tail_dist_to_curr <- dist_matrix[child_tail, curr_branch_cell_names]
            closest_to_tail <- names(tail_dist_to_curr)[which(tail_dist_to_curr == 
                min(tail_dist_to_curr))]
            tail_dist_to_anchored_branch = NA
            branch_index_for_tail <- NA
            tail_dist_to_anchored_branch <- dist_matrix[closest_to_tail, 
                child_tail]
            if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch) {
                reverse_child <- TRUE
            }
            else {
                reverse_child <- FALSE
            }
            res <- extract_branched_ordering_helper(branch_tree, 
                child, child_cell_ordering_subtree, branch_pseudotimes, 
                dist_matrix, reverse_child)
            child_cell_ordering_subtree <- res$subtree
            child_subtree_root <- res$root
            for (v in V(child_cell_ordering_subtree)) {
                cell_ordering_tree <- cell_ordering_tree + vertex(V(child_cell_ordering_subtree)[v]$name)
            }
            edge_list <- get.edgelist(child_cell_ordering_subtree)
            for (i in 1:nrow(edge_list)) {
                cell_ordering_tree <- cell_ordering_tree + edge(V(cell_ordering_tree)[edge_list[i, 
                  1]]$name, V(cell_ordering_tree)[edge_list[i, 
                  2]]$name)
            }
            if (tail_dist_to_anchored_branch < head_dist_to_anchored_branch) {
                cell_ordering_tree <- cell_ordering_tree + edge(closest_to_tail, 
                  child_subtree_root)
            }
            else {
                cell_ordering_tree <- cell_ordering_tree + edge(closest_to_head, 
                  child_subtree_root)
            }
        }
        return(list(subtree = cell_ordering_tree, root = curr_branch_root_cell, 
            last_cell_state = 1, last_cell_pseudotime = 0))
    }
    res <- extract_branched_ordering_helper(branch_tree, curr_branch, 
        cell_ordering_tree, branch_pseudotimes, dist_matrix, 
        reverse_main_path)
    cell_ordering_tree <- res$subtree
    curr_state <- 1
    assign_cell_state_helper <- function(ordering_tree_res, curr_cell) {
        nei <- NULL
        cell_tree <- ordering_tree_res$subtree
        V(cell_tree)[curr_cell]$cell_state = curr_state
        children <- V(cell_tree)[suppressWarnings(nei(curr_cell, 
            mode = "out"))]
        ordering_tree_res$subtree <- cell_tree
        if (length(children) == 1) {
            ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, 
                V(cell_tree)[children]$name)
        }
        else {
            for (child in children) {
                curr_state <<- curr_state + 1
                ordering_tree_res <- assign_cell_state_helper(ordering_tree_res, 
                  V(cell_tree)[child]$name)
            }
        }
        return(ordering_tree_res)
    }
    res <- assign_cell_state_helper(res, res$root)
    assign_pseudotime_helper <- function(ordering_tree_res, dist_matrix, 
        last_pseudotime, curr_cell) {
        nei <- NULL
        cell_tree <- ordering_tree_res$subtree
        curr_cell_pseudotime <- last_pseudotime
        V(cell_tree)[curr_cell]$pseudotime = curr_cell_pseudotime
        V(cell_tree)[curr_cell]$parent = V(cell_tree)[suppressWarnings(nei(curr_cell, 
            mode = "in"))]$name
        ordering_tree_res$subtree <- cell_tree
        children <- V(cell_tree)[suppressWarnings(nei(curr_cell, 
            mode = "out"))]
        for (child in children) {
            next_node <- V(cell_tree)[child]$name
            delta_pseudotime <- dist_matrix[curr_cell, next_node]
            ordering_tree_res <- assign_pseudotime_helper(ordering_tree_res, 
                dist_matrix, last_pseudotime + delta_pseudotime, 
                next_node)
        }
        return(ordering_tree_res)
    }
    res <- assign_pseudotime_helper(res, dist_matrix, 0, res$root)
    cell_names <- V(res$subtree)$name
    cell_states <- V(res$subtree)$cell_state
    cell_pseudotime <- V(res$subtree)$pseudotime
    cell_parents <- V(res$subtree)$parent
    ordering_df <- data.frame(sample_name = cell_names, cell_state = factor(cell_states), 
        pseudo_time = cell_pseudotime, parent = cell_parents)
    ordering_df <- plyr::arrange(ordering_df, pseudo_time)
    return(list(ordering_df = ordering_df, cell_ordering_tree = cell_ordering_tree))
}
