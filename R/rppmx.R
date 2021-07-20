#Miscellanea of useful/maybenot functions

#' postquant
#'
# @export
#'
postquant <- function(y, output, data, lab, plot){#, minbinder = F){
  cls <- as.matrix(output$label)
  psm <- mcclust::comp.psm(cls)
  #if(minbinder == T){
  #  mc <- minbinder.ext(psm)
  #} else {
    mc <- minVI(psm)
  #  }
  yhat <- output$pred
  vec <- (c(y)-c(yhat))
  mmse <- mean((y-yhat)^2)
  ari <- adjustedRandIndex(mc$cl, data$label)
  ess <- effectiveSize(output$nclu)

  mypostquant <- list("nclupost" = mean(output$nclu), "MSE" = mmse,
                      "lpml" = output$lpml, "ARI" = ari, "ESS" = ess)
                      #"comp_time" = time[1],
  if(lab==T){
    mypostquant <- list("nclupost" = mean(output$nclu), "MSE" = mmse,
                        "lpml" = output$lpml, "ARI" = ari,
                        #"comp_time" = time[1],
                        "ESS" = ess, "lab" = mc$cl)
  }
  if(plot == T){
    par(mfrow=c(1,2))
    plot(output$nclu, type="l")
    acf(output$nclu)
  }
  return(mypostquant)
}

#' postquant Multinomial Dirichlet
#'
# @export
#'
postquant_dm <- function(y, yp, output, data, plot, minbinder = F){
  cls <- as.matrix(output$label)
  psm <- comp.psm(cls)
  #if(minbinder == T){
    mc_b <- minbinder.ext(psm)
  #} else {
    mc_vi <- minVI(psm)
  #}
  ari_b <- adjustedRandIndex(mc_b$cl, data$labeltrain)
  ari_vi <- adjustedRandIndex(mc_vi$cl, data$labeltrain)
  ess <- effectiveSize(output$nclu)

  #predictive AUC multiclass
  probs <- apply(output$pipred, c(1, 2), mean)
  colnames(probs) <- colnames(Ytest) <- c(1:4)

  catvec <- c()
  for(i in 1:nrow(Ytest)){
    catvec <- c(catvec, which(yp[i,] == 1))
  }

  auc <- multiclass.roc(catvec, probs)$auc[1]
  mypostquant <- list("nclupost" = mean(output$nclu),
                      "ARI_B" = ari_b, "ARI_VI" = ari_vi, "ESS" = ess, "labb" = mc_b$cl, "labvi" = mc_vi$cl,
                      "waic" = output$WAIC, "lpml" = output$lpml, "auc" = auc)
  if(plot == T){
    par(mfrow=c(1,2))
    plot(output$nclu, type="l")
    acf(output$nclu)
  }
  return(mypostquant)
}

postquant_dm_ct <- function(y, yp, output, data, plot){
  ari_b <- ari_vi <- ess <- c()
  label <- output$label
  orig_names <- names(label)
  label <- label[order(names(label))]
  truelabel <- output$asstreat
  num_treat <- output$num_treat
  ntreat <- data$ntreat
  for(t in 1:ntreat){
    currtreat <- as.character(data$treatnames[[t]])
    mat <- label[[currtreat]]
    mat <- mat[1:num_treat[[currtreat]],]
    cls <- as.matrix(t(mat))
    psm <- comp.psm(cls)
    mc_b <- minbinder.ext(psm)
    mc_vi <- minVI(psm)
    ari_b[t] <- adjustedRandIndex(mc_b$cl, data$labeltrain[which(data$treat == currtreat)])
    ari_vi[t] <- adjustedRandIndex(mc_vi$cl, data$labeltrain[which(data$treat == currtreat)])
    ess[t] <- effectiveSize(output$nclu[t,])
  }

  #
  A0 <- output$pipred
  A0 <- array(unlist(A0), dim = c(dim(A0[[1]]), length(A0)))
  A0 <- apply(A0, c(1, 2, 3), mean)
  Al <- list()
  for(t in 1:dim(A0)[3]){
    Al[[t]] <- A0[,,t]
  }
  cat("og: ", orig_names, "\n")
  names(Al) <- orig_names
  Al <- Al[order(names(Al))]

  tpred <- c()
  for(i in 1:dim(A0)[1]){
    myvec <- c(A0[i,,])
    cat("myvec", myvec, "\n")
    tpred[i] <- switch(which(myvec == max(myvec)), "A1", "A2", "B1", "B2")
    cat("sel", tpred[i], "\n")
  }

  ##predictive AUC multiclass
  #probs <- apply(output$pipred, c(1, 2), mean)
  #colnames(probs) <- colnames(Ytest) <- c(1:4)

  #catvec <- c()
  #for(i in 1:nrow(Ytest)){
  #  catvec <- c(catvec, which(yp[i,] == 1))
  #}

  #auc <- multiclass.roc(catvec, probs)$auc[1]
  mypostquant <- list("nclupost" = apply(output$nclu, 1, mean),
                      "ARI_B" = ari_b, "ARI_VI" = ari_vi, "ESS" = ess, "labb" = mc_b$cl, "labvi" = mc_vi$cl,
                      "waic" = output$WAIC, "lpml" = output$lpml, "pipredm" = A0,
                      "pipredl" = Al, "mypred" = tpred)#, "auc" = auc)
  #if(plot == T){
  #  par(mfrow=c(1,2))
  #  plot(output$nclu, type="l")
  #  acf(output$nclu)
  #}
  return(mypostquant)
}

#' Multiclass AUC
#'
# @export
#'
plot_auc <- function(output_ppm, output_ppmx){
  probs_ppm <- apply(output_ppm$pipred, c(1, 2), mean)
  probs_ppmx <- apply(output_ppmx$pipred, c(1, 2), mean)

  colnames(Ytest) <- c("a_true", "b_true", "c_true", "d_true")
  colnames(probs_ppm) <- c("a_pred_ppm", "b_pred_ppm", "c_pred_ppm", "d_pred_ppm")
  colnames(probs_ppmx) <- c("a_pred_ppmx", "b_pred_ppmx", "c_pred_ppmx", "d_pred_ppmx")
  final_df <- data.frame(cbind(Ytest, probs_ppm, probs_ppmx))

  roc_res <- multi_roc(final_df, force_diag = T)
  plot_roc_df <- plot_roc_data(roc_res)

  ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
    geom_path(aes(color = Group, linetype=Method), size=0.5)
}

#' Multiclass AUC noplot
#'
# @export
#'
compute_auc <- function(output_ppm, output_ppmx){
  probs_ppm <- apply(output_ppm$pipred, c(1, 2), mean)
  probs_ppmx <- apply(output_ppmx$pipred, c(1, 2), mean)

  colnames(Ytest) <- c("a_true", "b_true", "c_true", "d_true")
  colnames(probs_ppm) <- c("a_pred_ppm", "b_pred_ppm", "c_pred_ppm", "d_pred_ppm")
  colnames(probs_ppmx) <- c("a_pred_ppmx", "b_pred_ppmx", "c_pred_ppmx", "d_pred_ppmx")
  final_df <- data.frame(cbind(Ytest, probs_ppm, probs_ppmx))

  roc_res <- multi_roc(final_df, force_diag = T)
}

#' post processing sketch as described in Richardson & Greene (1997)
#'
# @export
#'
pp_cs_rg <- function(output, post, nout, dim, refdim = 1){
  C_binder <- max(post$lab)
  median_C <- median(output$nclu)
  nc <- C_binder#median_C#

  eta_work <- matrix(0, nc, dim,)
  eta_pp <- matrix(0, nc, dim)
  noutC <- 0
  for(l in 1:nout){
    if(output$nclu[l]==nc){
      eta_work <- output$eta[c(1:nc),,l]
      eta_pp = eta_pp + eta_work[order(eta_work[,refdim]),]
      noutC <- noutC + 1
    }
  }
  eta_pp <- eta_pp/noutC
  return(eta_pp)
}

#' avg AUC/ROC
#' I am sure there is a better way
#'
# @export
#'
avg_auc <- function(roc_res, KK){

  roc_res_avg <- roc_res[[1]]
  for(kk in 2:KK){
    ##PPM
    roc_res_avg$Specificity$ppm$a <- roc_res_avg$Specificity$ppm$a +
      roc_res[[kk]]$Specificity$ppm$a
    roc_res_avg$Specificity$ppm$b <- roc_res_avg$Specificity$ppm$b +
      roc_res[[kk]]$Specificity$ppm$b
    roc_res_avg$Specificity$ppm$c <- roc_res_avg$Specificity$ppm$c +
      roc_res[[kk]]$Specificity$ppm$c
    roc_res_avg$Specificity$ppm$d <- roc_res_avg$Specificity$ppm$d +
      roc_res[[kk]]$Specificity$ppm$d

    roc_res_avg$Sensitivity$ppm$a <- roc_res_avg$Sensitivity$ppm$a +
      roc_res[[kk]]$Sensitivity$ppm$a
    roc_res_avg$Sensitivity$ppm$b <- roc_res_avg$Sensitivity$ppm$b +
      roc_res[[kk]]$Sensitivity$ppm$b
    roc_res_avg$Sensitivity$ppm$c <- roc_res_avg$Sensitivity$ppm$c +
      roc_res[[kk]]$Sensitivity$ppm$c
    roc_res_avg$Sensitivity$ppm$d <- roc_res_avg$Sensitivity$ppm$d +
      roc_res[[kk]]$Sensitivity$ppm$d

    roc_res_avg$AUC$ppm$a <- roc_res_avg$AUC$ppm$a +
      roc_res[[kk]]$AUC$ppm$a
    roc_res_avg$AUC$ppm$b <- roc_res_avg$AUC$ppm$b +
      roc_res[[kk]]$AUC$ppm$b
    roc_res_avg$AUC$ppm$c <- roc_res_avg$AUC$ppm$c +
      roc_res[[kk]]$AUC$ppm$c
    roc_res_avg$AUC$ppm$d <- roc_res_avg$AUC$ppm$d +
      roc_res[[kk]]$AUC$ppm$d

    ##PPMx
    roc_res_avg$Specificity$ppmx$a <- roc_res_avg$Specificity$ppmx$a +
      roc_res[[kk]]$Specificity$ppmx$a
    roc_res_avg$Specificity$ppmx$b <- roc_res_avg$Specificity$ppmx$b +
      roc_res[[kk]]$Specificity$ppmx$b
    roc_res_avg$Specificity$ppmx$c <- roc_res_avg$Specificity$ppmx$c +
      roc_res[[kk]]$Specificity$ppmx$c
    roc_res_avg$Specificity$ppmx$d <- roc_res_avg$Specificity$ppmx$d +
      roc_res[[kk]]$Specificity$ppmx$d

    roc_res_avg$Sensitivity$ppmx$a <- roc_res_avg$Sensitivity$ppmx$a +
      roc_res[[kk]]$Sensitivity$ppmx$a
    roc_res_avg$Sensitivity$ppmx$b <- roc_res_avg$Sensitivity$ppmx$b +
      roc_res[[kk]]$Sensitivity$ppmx$b
    roc_res_avg$Sensitivity$ppmx$c <- roc_res_avg$Sensitivity$ppmx$c +
      roc_res[[kk]]$Sensitivity$ppmx$c
    roc_res_avg$Sensitivity$ppmx$d <- roc_res_avg$Sensitivity$ppmx$d +
      roc_res[[kk]]$Sensitivity$ppmx$d

    roc_res_avg$AUC$ppmx$a <- roc_res_avg$AUC$ppmx$a +
      roc_res[[kk]]$AUC$ppmx$a
    roc_res_avg$AUC$ppmx$b <- roc_res_avg$AUC$ppmx$b +
      roc_res[[kk]]$AUC$ppmx$b
    roc_res_avg$AUC$ppmx$c <- roc_res_avg$AUC$ppmx$c +
      roc_res[[kk]]$AUC$ppmx$c
    roc_res_avg$AUC$ppmx$d <- roc_res_avg$AUC$ppmx$d +
      roc_res[[kk]]$AUC$ppmx$d
  }

  roc_res_avg$Specificity$ppm$a <- roc_res_avg$Specificity$ppm$a/KK
  roc_res_avg$Specificity$ppm$b <- roc_res_avg$Specificity$ppm$b/KK
  roc_res_avg$Specificity$ppm$c <- roc_res_avg$Specificity$ppm$c/KK
  roc_res_avg$Specificity$ppm$d <- roc_res_avg$Specificity$ppm$d/KK

  roc_res_avg$Sensitivity$ppm$a <- roc_res_avg$Sensitivity$ppm$a/KK
  roc_res_avg$Sensitivity$ppm$b <- roc_res_avg$Sensitivity$ppm$b/KK
  roc_res_avg$Sensitivity$ppm$c <- roc_res_avg$Sensitivity$ppm$c/KK
  roc_res_avg$Sensitivity$ppm$d <- roc_res_avg$Sensitivity$ppm$d/KK

  roc_res_avg$AUC$ppm$a <- roc_res_avg$AUC$ppm$a/KK
  roc_res_avg$AUC$ppm$b <- roc_res_avg$AUC$ppm$b/KK
  roc_res_avg$AUC$ppm$c <- roc_res_avg$AUC$ppm$c/KK
  roc_res_avg$AUC$ppm$d <- roc_res_avg$AUC$ppm$d/KK

  roc_res_avg$Specificity$ppmx$a <- roc_res_avg$Specificity$ppmx$a/KK
  roc_res_avg$Specificity$ppmx$b <- roc_res_avg$Specificity$ppmx$b/KK
  roc_res_avg$Specificity$ppmx$c <- roc_res_avg$Specificity$ppmx$c/KK
  roc_res_avg$Specificity$ppmx$d <- roc_res_avg$Specificity$ppmx$d/KK

  roc_res_avg$Sensitivity$ppmx$a <- roc_res_avg$Sensitivity$ppmx$a/KK
  roc_res_avg$Sensitivity$ppmx$b <- roc_res_avg$Sensitivity$ppmx$b/KK
  roc_res_avg$Sensitivity$ppmx$c <- roc_res_avg$Sensitivity$ppmx$c/KK
  roc_res_avg$Sensitivity$ppmx$d <- roc_res_avg$Sensitivity$ppmx$d/KK

  roc_res_avg$AUC$ppmx$a <- roc_res_avg$AUC$ppmx$a/KK
  roc_res_avg$AUC$ppmx$b <- roc_res_avg$AUC$ppmx$b/KK
  roc_res_avg$AUC$ppmx$c <- roc_res_avg$AUC$ppmx$c/KK
  roc_res_avg$AUC$ppmx$d <- roc_res_avg$AUC$ppmx$d/KK

  roc_res_avg$Specificity$ppm <- within(roc_res_avg$Specificity$ppm, rm(macro, micro))
  roc_res_avg$Specificity$ppmx <- within(roc_res_avg$Specificity$ppmx, rm(macro, micro))
  roc_res_avg$Sensitivity$ppm <- within(roc_res_avg$Sensitivity$ppm, rm(macro, micro))
  roc_res_avg$Sensitivity$ppmx <- within(roc_res_avg$Sensitivity$ppmx, rm(macro, micro))
  roc_res_avg$AUC$ppm <- within(roc_res_avg$AUC$ppm, rm(macro, micro))
  roc_res_avg$AUC$ppmx <- within(roc_res_avg$AUC$ppmx, rm(macro, micro))

  return(roc_res_avg)
}

my_plot_multi_roc_data <- function (roc_res){
  n_method <- length(unique(roc_res$Methods))
  n_group <- length(unique(roc_res$Groups))
  roc_res_df <- data.frame(Specificity = numeric(0), Sensitivity = numeric(0),
                           Group = character(0), AUC = numeric(0), Method = character(0))
  for (i in 1:n_method) {
    for (j in 1:n_group) {
      temp_data_1 <- data.frame(Specificity = roc_res$Specificity[[i]][j],
                                Sensitivity = roc_res$Sensitivity[[i]][j], Group = unique(roc_res$Groups)[j],
                                AUC = roc_res$AUC[[i]][j], Method = unique(roc_res$Methods)[i])
      colnames(temp_data_1) <- c("Specificity", "Sensitivity",
                                 "Group", "AUC", "Method")
      roc_res_df <- rbind(roc_res_df, temp_data_1)
    }
  }
  return(roc_res_df)
}
