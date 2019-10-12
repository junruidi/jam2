#' @title Multilelvel Principal Component Analysis
#' @description Conduct multilevel (functional) principal component analysis.
#' A mixed model framework is used to estimate scores and obtain variance estimates.
#'
#' @param Y \code{matrix} with columns as variables/grid points, and rows as subject-visit, must be supplied.
#' @param id Must be supplied, a vector containing the id information used to identify clusters.
#' @param twoway logical, indicating whether to carry out twoway ANOVA and calculate
#' visit-specific means. Defaults to \code{FALSE}.
#' @param cov.method covariance estimation approaches.
#' @param pve proportion of variance explained: used to choose the number of principal components.
#' @param tlength For functional data,the length of the interval that the functions are
#' evaluated at. Defualt to be number of grid points which is exactly number of variables
#' for multivariate data.
#' @param smoothing logical, indicating whether smoothing should be carried out or not, only use this
#' for functional data.
#' @param smooth.method Method of smoothing, "sf" indicates directly smoothing the functions using gam.
#' "sc" indicates smoothing the covariance using penalized smoothing splines.
#' @param nk number of knots for the smoothing
#'
#' @return An object of class \code{mpca} containing:
#' \item{Y.df}{The original data.}
#' \item{mu}{estimated mean function.}
#' \item{eta}{the estimated visit specific shifts from overall mean.}
#' \item{npc}{number of PCs.}
#' \item{pctvar}{percentage of variation explained by the kept pcs within each level.}
#' \item{varofTot}{percentage of variation of the total variability explained by level 1 and level 2.}
#' \item{evectors}{level 1/2 eigen vectors/functions.}
#' \item{evalues}{level 1/2 eigen values.}
#' \item{scores}{level 1/2 pc scores.}
#'
#' @importFrom mgcv gam
#' @importFrom dplyr %>% group_by slice
#' @importFrom ICSNP pair.diff
#' @importFrom SemiPar spm predict.spm
#' @importFrom  MASS ginv
#' @importFrom stats ave
#'
#' @references Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
#' Multilevel functional principal component analysis. \emph{Annals of Applied
#' Statistics}, 3, 458--488.
#'
#' Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2014).
#' Multilevel sparse functional principal component analysis. \emph{Stat}, 3, 126--143.
#'
#' Goldsmith, J., Greven, S., and Crainiceanu, C. (2013). Corrected confidence
#' bands for functional data using principal components. \emph{Biometrics},
#' 69(1), 41--51.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mpca.y = mfpca.sc(Y = Y, id = id, twoway = T, cov.method = "m1")
#' }


multilevel_pca = function(Y = NULL, id = NULL, twoway = TRUE,
                          cov.method = c("m1","m2"), pve = 0.9,
                          tlength = ncol(Y),
                          smoothing = FALSE, smooth.method = c("sf","sc") ,nk = 15){

  x1 = x2 = myknots = NULL
  rm(list = c("x1","x2","myknots"))
  # 1. Preliminary check on all possible mistakes ---------------------------
  stopifnot((!is.null(Y) && !is.null(id)))
  if(length(id) != nrow(Y)){
    stop("number of id visits should be the same as number of rows of Y")
  }
  cov.method = match.arg(cov.method)
  if(smoothing == TRUE){
    tlength = 1
  }

  if(smoothing  == FALSE){
    smooth.method = "nothing"
    nk = NULL
  }

  if(smoothing == TRUE & is.null(smooth.method)){
    stop("choosing the smoothing method")
  }
  if(smoothing == TRUE){
    cov.method = match.arg(cov.method)
  }

  if(smoothing == TRUE & smooth.method == "sf"){
    Y_sm = matrix(NA,ncol=ncol(Y),nrow=nrow(Y))
    tY = 1:ncol(Y_sm)

    for(i in 1:nrow(Y_sm)){
      Y_sm[i,] = gam(Y[i,] ~ s(tY,k = nk))$fitted.values
    }
    Y = Y_sm
  }


  #########################################################################################

  # 2. Check the dimension and preliminary ----------------------------------
  uni.id = unique(id)
  visit = ave(id, id, FUN = seq_along)
  Y.df = data.frame(id = id, visit = visit)
  Y.df$Y = Y



  J = length(unique(Y.df$visit)) # largest number of visit
  M = length(unique(Y.df$id))    # number of subjects
  N = NCOL(Y.df$Y)               # number of grid points/variables
  I = NROW(Y.df$Y)               # number of subject-visits

  num_visit = c(t(table(Y.df$id))) # number of visit for each subjects
  cum_visit_number = cumsum(num_visit)

  kx = tlength/N


  #########################################################################################


  # 3. Estimate fixed effects \mu and eta_j ---------------------------------

  # estimate grand mean mu
  mu = apply(Y, 2,mean)

  # estimate eta_j
  mueta = matrix(0, J, N)
  eta = matrix(0, J, N)
  if (twoway) {
    Y.split = as.list(rep(NA, J))
    for (j in 1:J) {
      Y.split[[j]] = subset(Y.df, visit == j)
      mueta[j, ] = apply(Y.split[[j]]$Y, 2, mean)
      eta[j,] = mueta[j,] - mu
      Y.split[[j]]$Y.tilde = Y.split[[j]]$Y - matrix(mueta[j,], NROW(Y.split[[j]]), N, byrow = TRUE)
    }
    Y.df.new = Reduce(function(...) merge(..., by = c("id", "visit", "Y", "Y.tilde"), all = TRUE, sort = FALSE), Y.split)
    Y.df.new = Y.df.new[order(Y.df.new$id, Y.df.new$visit),  ]
    Y.tilde = Y.df.new$Y.tilde
  }else{
    Y.df.new = Y.df[order(Y.df$id, Y.df$visit),  ]
    Y.df.new$Y.tilde = Y.df.new$Y - matrix(mu, I, N, byrow = TRUE)
    Y.tilde = Y.df.new$Y.tilde
  }

  #########################################################################################


  # 4. Estimate the covariates matrix/surface Gb and Gw ---------------------
  ###     Estimate the three covariance functions: overall covariance G,
  ###     between covariance Gb and within covariance Gw

  # covariance estimation method 1: all J_i terms inside sum i term, see Wang 2018
  if(cov.method == "m1"){
    Gt = matrix(0, N, N)
    for(m in 1:M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])$Y.tilde
      Gt = Gt + t(Y.tilde.m) %*% Y.tilde.m/nrow(Y.tilde.m)
    }
    Gt = Gt/M

    Gw.temp =  matrix(0, N, N)
    for (m in 1: M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])$Y.tilde
      Ji = nrow(Y.tilde.m)
      if(Ji > 1){
        resid.pair.temp = Y.tilde.m[c(1:Ji), ]
        resid.pair.diff = pair.diff(resid.pair.temp)
        Gw.temp = Gw.temp + crossprod(resid.pair.diff)/nrow(resid.pair.diff)
      }
    }

    Gw = Gw.temp/(2*sum(num_visit > 1))


    # Gb = matrix(0, N, N)
    # for(m in 1:M){
    #   Y.tilde.m = subset(Y.df.new, id == uni.id[m])$Y.tilde
    #   Ji = nrow(Y.tilde.m)
    #   if(Ji > 1){
    #     gb.temp = matrix(0,N,N)
    #     for(j1 in 1:Ji){
    #       for(j2 in setdiff(c(1:Ji),j1)){
    #         gb.temp = gb.temp + Y.tilde.m[j1,] %*% t(Y.tilde.m[j2,])
    #       }
    #     }
    #     Gb = Gb + gb.temp/(Ji * (Ji-1))
    #   }
    # }
    # M_gt1j = sum(num_visit > 1)
    # Gb = Gb/M_gt1j


    Gb = Gt - Gw
  }

  #  covariance estimation method 2: based on MoM proposed in Koch 1968, Shou 2013.
  if(cov.method == "m2"){
    num_visit_gt1 = num_visit[which(num_visit > 1)]
    Gw = matrix(0, N, N)
    for(m in 1:M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])$Y.tilde
      Ji = nrow(Y.tilde.m)
      if(Ji > 1){
        # for(j1 in 1:Ji){
        #   for(j2 in setdiff(c(1:Ji),j1)){
        #     diff.j1j2 = Y.tilde.m[j1,] - Y.tilde.m[j2,]
        #     Gw = Gw + diff.j1j2 %*% t(diff.j1j2)
        #   }
        # }
        resid.pair.temp = Y.tilde.m[c(1:Ji), ]
        resid.pair.diff = pair.diff(resid.pair.temp)
        Gw = Gw + crossprod(resid.pair.diff) * 2
      }
    }

    denom_1 = sum(num_visit_gt1^2) - sum(num_visit_gt1)
    Gw = 0.5/denom_1 * Gw

    h1 = matrix(0,N,N)
    for(m in 1:M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])$Y.tilde
      Ji = nrow(Y.tilde.m)
      temp = matrix(0,N,N)
      mult = I - Ji
      for(j in 1:Ji){
        temp = temp + Y.tilde.m[j,] %*% t(Y.tilde.m[j,])
      }
      h1 = h1 + mult * temp
    }

    g1 = colSums(Y.df.new$Y.tilde) %*% t(colSums(Y.df.new$Y.tilde))
    g2 = matrix(0,N,N)
    for(m in 1:M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])$Y.tilde
      y.ij.sum = colSums(Y.tilde.m)
      g2 = g2 + y.ij.sum %*% t(y.ij.sum)
    }
    gm = g1 - g2
    denom_2 = I^2 - sum(num_visit^2)
    Hz = 2 * h1/denom_2 - 2 * gm/denom_2

    Gb = Hz/2 - Gw
    Gt = Hz/2
  }


  #########################################################################################


  # 5. Smoothing the covariance surface Gb and Gw for functional dat --------

  ###     If "smoothing=TRUE" (which is recommended when functions are measured with noise),
  ###     this option will smooth Gw(s,t) and Gb(s,t) before carrying out principal component
  ###     analysis. This step will reduce the bias in estimating eigenvalues and generate
  ###     smooth eigenfunctions.
  if(smoothing == TRUE & smooth.method == "sc") {
    gw.temp = Gw
    diag(gw.temp) = rep(NA, N)

    data.gb = data.frame(gb = as.vector(Gb),gw = as.vector(gw.temp),
                         x1 = rep(seq(0,1,length = N), N), x2 = rep(seq(0,1,length = N), each = N))
    attach(data.gb)

    myknots = data.frame(x1 = rep(seq(0,1,length = nk),each = nk), x2 = rep(seq(0,1,length = nk),nk))




    # attach(data.gb)
    fit1 = spm(gb ~ f(x1, x2,knots = myknots))
    fit =  spm(gw ~ f(x1, x2,knots = myknots),omit.missing=T)
    newdata <- data.frame(x1=x1,x2=x2)
    pred1 = predict.spm(fit1,newdata)
    pred = predict.spm(fit,newdata)




    var.noise = mean( diag(Gw) - diag(matrix(pred,N,N)) )
    s.gw = matrix(pred,N,N)
    Gw = (s.gw + t(s.gw) )/2
    s.gb = matrix(pred1, N,N)
    Gb = (s.gb + t(s.gb))/2
  }


  #########################################################################################

  # 6. Create level 1/2 eigenvectors/functions and eigen values -------------
  ###     Estimate eigen values and eigen functions at two levels
  e1 = eigen(Gb)
  e2 = eigen(Gw)

  ###     Estimate amount of variability explained by X and W
  pw = sum(diag(Gw))/sum(diag(Gt))
  px = 1 - pw



  ###    get eigen values
  fpca1.value = e1$values * kx
  fpca2.value = e2$values * kx



  ###     Keep only non-negative eigenvalues
  fpca1.value = ifelse(fpca1.value>=0, fpca1.value, 0)
  fpca2.value = ifelse(fpca2.value>=0, fpca2.value, 0)

  ###     Calculate the percentage of variance that are explained by the components
  percent1 = (fpca1.value)/sum(fpca1.value)
  percent2 = (fpca2.value)/sum(fpca2.value)

  ###     Decide the number of components that are kept at level 1 and 2. The general
  ###     rule is to stop at the component where the cumulative percentage of variance
  ###     explained is greater than 90% and the variance explained by any single component
  ###     after is less than 1/N.
  K1 = max( which(cumsum(percent1) < pve | percent1 > 1/N ) + 1)
  K2 = max( which(cumsum(percent2) < pve | percent2 > 1/N ) + 1)

  ###     estimate eigen vectors
  fpca1.vectors = e1$vectors[, 1:K1] * sqrt(1/kx)
  fpca2.vectors = e2$vectors[, 1:K2] * sqrt(1/kx)


  ###     The eigen vectosrs are unique only up to a change of signs.
  ###     Select the signs of eigenfunctions so that the integration over the domain
  ###     is non-negative
  for(i in 1:K1) {
    v2 = fpca1.vectors[,i]
    tempsign = sum(v2)
    fpca1.vectors[,i] = ifelse(tempsign<0, -1,1) * v2
  }
  for(i in 1:K2) {
    v2 = fpca2.vectors[,i]
    tempsign = sum(v2)
    fpca2.vectors[,i] = ifelse(tempsign<0, -1,1) * v2
  }



  #########################################################################################

  # 7. Get level 1/2 PC scores based on the projection methods by so --------

  id.visit.dat = Y.df.new[,c(1:2)]
  id.visit.dat = as.data.frame(cbind(id.visit.dat,Y.df.new$Y.tilde))

  id.visit.all = data.frame(id = rep(uni.id, each = J), visit = rep(c(1:J), M))

  dat = merge(x = id.visit.all, y = id.visit.dat, all.x = T)
  row.index.add = which(rowSums(is.na(dat)) == N)
  dat[row.index.add,-c(1:2)] = 0
  resid.temp2 = as.matrix(dat[,-c(1:2)])


  cross.integral = t(fpca1.vectors)%*%fpca2.vectors * kx # Phi_w * t(Phik)
  int1 = matrix(0, M*J, K1)
  int2 = matrix(0, M*J, K2)
  for(i in 1:(M*J))   {
    for(j in 1:K1) int1[ i ,j] = sum( resid.temp2[i,] * fpca1.vectors[,j] ) * kx
    for(j in 1:K2) int2[ i ,j] = sum( resid.temp2[i,] * fpca2.vectors[,j] ) * kx
  }

  s1 = matrix(NA, M*J, K1)
  s2 = matrix(NA, M*J, K2)

  design.xi = ginv( diag(rep(1,K1)) - cross.integral %*% t(cross.integral) )
  for(m in 1:M) {
    resd = rep(0, K1)
    for(j in 1:J) {
      index = (m-1) * J + j
      resd = resd + ( int1[index,] - drop(cross.integral %*% int2[index,]) )/num_visit[m]
    }
    index.m = ( (m-1) * J + 1 ) : (m*J)
    xi.temp = design.xi %*% resd
    s1[index.m,] = matrix(rep(xi.temp, each=J), nrow=J)
    s2[index.m,] = t(t(int2[index.m,]) - drop( t(cross.integral) %*% xi.temp ))
  }

  s1 = as.data.frame(cbind(dat[,1:2],s1)) %>% group_by(id) %>% slice(1)
  s1$visit = NULL
  names(s1)[2:ncol(s1)] = paste0("lv1_",names(s1)[2:ncol(s1)])

  s2 = as.data.frame(cbind(dat[,1:2],s2))
  s2 = merge(x = Y.df.new[,c(1:2)], y = s2, all.x = T)
  names(s2)[3:ncol(s2)] = paste0("lv2_",names(s2)[3:ncol(s2)])
  # s2$id = s2$visit = NULL





  #########################################################################################
  # 8. Create the result list as a class ------------------------------------


  npc = list(K1,K2)
  evalues = list(fpca1.value[1:K1], fpca2.value[1:K2])
  evectors = list(fpca1.vectors[,1:K1], fpca2.vectors[,1:K2])
  scores = list(s1, s2)
  pctvar = list(percent1[1:K1], percent2[1:K2])
  varofTot = list(px, pw)
  Y.df = Y.df.new[,-4]

  names(evectors) = names(evalues) = names(npc) = names(scores) = names(pctvar) = names(varofTot) = c("level1", "level2")
  ret.objects = c("Y.df", "mu", "eta","npc", "pctvar","varofTot","evectors", "evalues", "scores")
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  class(ret) = "mpca"
  return(ret)


}

