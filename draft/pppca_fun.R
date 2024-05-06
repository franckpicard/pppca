library(parallel)
library(RSpectra)
pppca = function(PP,Jmax,mc.cores=24){
    grid              = sort(unique(Reduce('c',PP)))
    dgrid             = diff(grid)
    nb_occ_by_process = sapply(PP, length)
    nb_occ_total      = length(grid)
    nb_processes      = length(PP)
    
    Lambda  = 0:(nb_occ_total-1) / nb_processes			
    Pi      = t(mcmapply(PP,FUN=function(x){
        c(0,cumsum(hist(x,breaks = grid, plot=F)$counts))
    },mc.cores=mc.cores))
    K_Delta = crossprod(Pi)/nb_processes - tcrossprod(Lambda)
    M_Delta = tcrossprod(sqrt(dgrid))*K_Delta[2:nb_occ_total,2:nb_occ_total]
    PP_PC   = eigs_sym(M_Delta,Jmax)
    
    eigenval   = PP_PC$values[1:Jmax]				
    
    baseline_pmf     = apply(Pi,2,mean)[-nb_occ_total]
    Pi_centered_norm = t(apply(Pi, 1, FUN=function(x){
        (x[-nb_occ_total]-baseline_pmf)*sqrt(dgrid)
    }))
    scores     = Pi_centered_norm%*%PP_PC$vectors[,1:Jmax]
    scores     = t(apply(scores,1,FUN=function(x){x/sqrt(eigenval)}))
    
    change_sign = sign(PP_PC$vectors[nb_occ_total/2,1])
    
    scores[,1] = change_sign*scores[,1]
    scores     = as.data.frame(scores)
    colnames(scores) = paste0("axis",c(1:Jmax))
    
    eigenfun        = PP_PC$vectors[,1:Jmax]/sqrt(dgrid)
    eigenfun[,1]    = change_sign*eigenfun[,1]
    eigenfun        = data.frame(eigenfun)
    names(eigenfun) = paste0("fun",c(1:Jmax))
    
    return(
        list(
            grid       = grid,
            eigenval   = eigenval,
            percentvar = 100* eigenval / sum(eigenval),
            eigenfun   = eigenfun, 
            scores     = scores
        )
    )
}

w_pred = function(i,baseline_pmf, eigenval, scores, eigenfun){
    J = dim(eigenfun)[2]
    w = baseline_pmf
    for (j in 1:J){
        w = w + sqrt(eigenval[j])*scores[i,j]*eigenfun[,j]
    }
    return(w)
}


get_pred_pmf = function(PP, eigenval, scores, eigenfun, mc.cores=24){
    n     = length(PP)
    grid  = sort(unique(Reduce('c',PP)))
    PPbin = mclapply(PP, grid, FUN=function(x,gg){	
        y = rep(0,length(gg))
        y[which(gg %in% x)] = 1
        return(y)
    },mc.cores=mc.cores)
    baseline_pmf = cumsum(Reduce('+',PPbin))[-1]#/n
    pred_pmf     = mcmapply(1:n,FUN=function(i){w_pred(i,baseline_pmf, eigenval, scores, eigenfun)},mc.cores=mc.cores)
    return(pred_pmf)
}

get_pred_int = function(PP, eigenval, scores, eigenfun, mc.cores=24){
    n     = length(PP)
    grid  = sort(unique(Reduce('c',PP)))
    PPbin = mclapply(PP, grid, FUN=function(x,gg){	
        y = rep(0,length(gg))
        y[which(gg %in% x)] = 1
        return(y)
    },mc.cores=mc.cores)
    baseline_pmf = cumsum(Reduce('+',PPbin))[-1]#/n
    pred_pmf     = mcmapply(1:n,FUN=function(i){w_pred(i,baseline_pmf, eigenval, scores, eigenfun)},mc.cores=mc.cores)
    pred_int     = apply(pred_pmf,2,FUN=function(x){diff(x)})
    return(pred_int)
}


get_pmf = function(PP, mc.cores=24){
    n     = length(PP)
    grid  = sort(unique(Reduce('c',PP)))
    PPbin = mclapply(PP, grid, FUN=function(x,gg){	
        y = rep(0,length(gg))
        y[which(gg %in% x)] = 1
        return(y)
    },mc.cores=mc.cores)
    pmf = mclapply(PPbin, FUN=function(x){return(cumsum(x)[-1])},mc.cores=mc.cores)
    return(pmf)
}

get_int = function(PP, mc.cores=24){
    n     = length(PP)
    grid  = sort(unique(Reduce('c',PP)))
    PPbin = mclapply(PP, grid, FUN=function(x,gg){	
        y = rep(0,length(gg))
        y[which(gg %in% x)] = 1
        return(y)
    },mc.cores=mc.cores)
    pmf = mclapply(PPbin, FUN=function(x){return(cumsum(x)[-1])},mc.cores=mc.cores)
    int = mclapply(pmf,FUN=function(x){diff(x)},mc.cores=mc.cores)
    return(int)
}

