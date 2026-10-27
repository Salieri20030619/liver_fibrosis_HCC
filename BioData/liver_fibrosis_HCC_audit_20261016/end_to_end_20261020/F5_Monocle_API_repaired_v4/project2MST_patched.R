function (cds, Projection_Method) 
{
    dp_mst <- minSpanningTree(cds)
    Z <- reducedDimS(cds)
    Y <- reducedDimK(cds)
    cds <- findNearestPointOnMST(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex_names <- colnames(Y)[closest_vertex]
    closest_vertex_df <- as.matrix(closest_vertex)
    row.names(closest_vertex_df) <- row.names(closest_vertex)
    tip_leaves <- names(which(degree(dp_mst) == 1))
    if (!is.function(Projection_Method)) {
        P <- Y[, closest_vertex]
    }
    else {
        P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
        for (i in 1:length(closest_vertex)) {
            neighbors <- names(V(dp_mst)[suppressWarnings(.nei(closest_vertex_names[i], 
                mode = "all"))])
            projection <- NULL
            distance <- NULL
            Z_i <- Z[, i]
            for (neighbor in neighbors) {
                if (closest_vertex_names[i] %in% tip_leaves) {
                  tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], 
                    neighbor)])
                }
                else {
                  tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], 
                    neighbor)])
                }
                projection <- rbind(projection, tmp)
                distance <- c(distance, dist(rbind(Z_i, tmp)))
            }
            if (!is(projection, "matrix")) 
                projection <- as.matrix(projection)
            P[, i] <- projection[which(distance == min(distance))[1], 
                ]
        }
    }
    colnames(P) <- colnames(Z)
    dp <- as.matrix(dist(t(P)))
    min_dist = min(dp[dp != 0])
    dp <- dp + min_dist
    diag(dp) <- 0
    cellPairwiseDistances(cds) <- dp
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
    cds
}
