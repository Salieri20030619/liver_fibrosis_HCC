function (object, newdata, comps = object$computed_nt, type = c("link", 
    "response", "terms", "scores", "class", "probs"), se.fit = FALSE, 
    weights, dispersion = NULL, methodNA = "adaptative", verbose = TRUE, 
    ...) 
{
    if (!inherits(object, "plsRglmmodel")) 
        stop("Primary argument much be a plsRglmmodel object")
    if (missing(type)) {
        type = "link"
    }
    type <- match.arg(type)
    if (!(type %in% c("link", "response", "terms", "scores", 
        "class", "probs"))) 
        stop("Invalid type specification")
    if (comps > object$computed_nt) {
        stop("Cannot predict using more components than extracted.")
    }
    is_pls <- isTRUE(identical(as.character(object$call$modele), 
        "pls")) || (is.null(object$FinalModel) && !is.null(object$CoeffC))
    predict_pls_from_scores <- function(scores_mat, type) {
        link_pred <- as.matrix(scores_mat) %*% object$CoeffC
        if (type == "link") {
            return(link_pred)
        }
        if (type == "response") {
            return(attr(object$RepY, "scaled:center") + attr(object$RepY, 
                "scaled:scale") * link_pred)
        }
        return(NULL)
    }
    predict_glm_from_scores <- function(scores_mat, type) {
        if (isTRUE(se.fit)) {
            ttpred <- data.frame(tt = scores_mat)
            return(predict(object$FinalModel, newdata = ttpred, 
                type = type, se.fit = se.fit, dispersion = dispersion, 
                ...))
        }
        pls_glm_predict_score_model(object$FinalModel, scores_mat, 
            type = type)
    }
    if (missing(newdata) || is.null(newdata)) {
        nrtt <- nrow(object$tt)
        ttpredY <- data.frame(cbind(object$tt[, 1:comps], matrix(0, 
            nrow = nrtt, ncol = object$computed_nt - comps)))
        colnames(ttpredY) <- NULL
        if (is_pls && (type %in% c("link", "response"))) {
            return(predict_pls_from_scores(ttpredY, type))
        }
        if (type == "link") {
            if ("glm" %in% class(object$FinalModel)) {
                return(predict_glm_from_scores(ttpredY, "link"))
            }
        }
        if (type == "response") {
            if ("glm" %in% class(object$FinalModel)) {
                return(predict_glm_from_scores(ttpredY, "response"))
            }
        }
        if (type == "terms") {
            if ("glm" %in% class(object$FinalModel)) {
                return(predict_glm_from_scores(ttpredY, "terms"))
            }
        }
        if (type == "class") {
            if ("polr" %in% class(object$FinalModel)) {
                ttpred <- data.frame(tt = ttpredY)
                return(predict(object$FinalModel, newdata = ttpred, 
                  type = "class", ...))
            }
        }
        if (type == "probs") {
            if ("polr" %in% class(object$FinalModel)) {
                ttpred <- data.frame(tt = ttpredY)
                return(predict(object$FinalModel, newdata = ttpred, 
                  type = "probs", ...))
            }
        }
        if (type == "scores") {
            return(object$tt[, 1:comps])
        }
    }
    else {
        nrnd <- nrow(newdata)
        if (any(apply(is.na(newdata), MARGIN = 1, "all"))) {
            return(vector("list", 0))
            cat("One of the rows of newdata is completely filled with missing data\n")
            stop()
        }
        if (any(is.na(newdata))) {
            na.miss.newdata <- TRUE
        }
        else {
            na.miss.newdata <- FALSE
        }
        if (!is.null(object$call$formula)) {
            mf <- match.call(expand.dots = FALSE)
            m <- match(c("subset", "weights"), names(mf), 0L)
            mf <- mf[c(1L, m)]
            mf$data <- newdata
            mf$formula <- object$call$formula[-2]
            mf$drop.unused.levels <- TRUE
            mf$na.action <- na.pass
            mf[[1L]] <- as.name("model.frame")
            mf <- eval(mf, parent.frame())
            mt <- attr(mf, "terms")
            newdata.frame <- if (!is.empty.model(mt)) 
                model.matrix(mt, mf, contrasts)[, -1]
            else matrix(, nrnd, 0L)
            weights <- as.vector(model.weights(mf))
        }
        else {
            newdata.frame <- newdata
        }
        newdata.scaled <- sweep(sweep(newdata.frame, 2, attr(object$ExpliX, 
            "scaled:center")), 2, attr(object$ExpliX, "scaled:scale"), 
            "/")
        newdataNA <- !is.na(newdata)
        newdata.scaledwotNA <- as.matrix(newdata.scaled)
        newdata.scaledwotNA[!newdataNA] <- 0
        ttpredY <- NULL
        if (methodNA == "adaptative") {
            for (ii in 1:nrnd) {
                if (all(newdataNA[ii, ])) {
                  ttpredY <- rbind(ttpredY, c(newdata.scaledwotNA[ii, 
                    ] %*% object$wwetoile[, 1:comps], rep(0, 
                    object$computed_nt - comps)))
                }
                else {
                  if (verbose) {
                    cat("Missing value in row ", ii, ".\n")
                  }
                  ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii, 
                    ], , drop = FALSE]) %*% object$pp[newdataNA[ii, 
                    ], , drop = FALSE]) %*% t(object$pp[newdataNA[ii, 
                    ], , drop = FALSE]) %*% (newdata.scaledwotNA[ii, 
                    ])[newdataNA[ii, ]])[1:comps], rep(0, object$computed_nt - 
                    comps)))
                }
            }
        }
        if (methodNA == "missingdata") {
            if (verbose) {
                cat("Prediction as if missing values in every row.\n")
            }
            for (ii in 1:nrnd) {
                ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii, 
                  ], , drop = FALSE]) %*% object$pp[newdataNA[ii, 
                  ], , drop = FALSE]) %*% t(object$pp[newdataNA[ii, 
                  ], , drop = FALSE]) %*% (newdata.scaledwotNA[ii, 
                  ])[newdataNA[ii, ]])[1:comps], rep(0, object$computed_nt - 
                  comps)))
            }
        }
        colnames(ttpredY) <- NULL
        ttpred <- data.frame(tt = ttpredY)
        if (is_pls && (type %in% c("link", "response"))) {
            return(predict_pls_from_scores(ttpredY, type))
        }
        if ("glm" %in% class(object$FinalModel)) {
            if (type == "link") {
                return(predict_glm_from_scores(ttpredY, "link"))
            }
            if (type == "response") {
                return(predict_glm_from_scores(ttpredY, "response"))
            }
            if (type == "terms") {
                return(predict_glm_from_scores(ttpredY, "terms"))
            }
        }
        if ("polr" %in% class(object$FinalModel)) {
            if (type == "class") {
                return(predict(object$FinalModel, newdata = ttpred, 
                  type = "class", ...))
            }
            if (type == "probs") {
                return(predict(object$FinalModel, newdata = ttpred, 
                  type = "probs", ...))
            }
        }
        if (type == "scores") {
            colnames(ttpredY) <- paste("Comp_", 1:object$computed_nt, 
                sep = "")
            return(ttpredY[, 1:comps, drop = FALSE])
        }
    }
}
