#' @title calcCONS
#' @description A function for calculating conservative LR
#' @details The function calculates LR based on MCMC simulations
#' @param mlefitHp A fitted object returned from contLikMLE (under Hp)
#' @param mlefitHd A fitted object returned from contLikMLE (under Hd)
#' @param niter Number of samples in the MCMC-sampling.
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Should be higher to obtain lower acception rate.
#' @param quantile The quantile used to report conservative LR 
#' @param seed The user can set seed if wanted
#' @param verbose Whether progress should be printed
#' @param accRateGolden aimed acceptance rate
#' @param diffTol Difference tolerance regarding acceptance rate
#' @export

calcCONS = function(mlefitHp, mlefitHd,niter=2000,delta=2,quantile=0.05,seed=NULL, verbose=TRUE, accRateGolden=0.25,diffTol=0.1) {
  #PERFORM CALIBRATING OF DELTA BEFORE RUNNING ALL SAMPLE
  #Tweak delta to find  suitable acceptance rate:
  if(verbose) print("Calibrating MCMC simulator...")
  while(TRUE) { 
    if(verbose) print(paste0("Check with delta=",delta))
    hpmcmc <- euroformix::contLikMCMC(mlefitHp,niter=200,delta=delta)
    acc0 = hpmcmc$accrat #obtain acceptance rate
    if(verbose) print(paste0("Acceptance rate=",acc0))
    if( abs(acc0-accRateGolden)<diffTol) break
    
    if(acc0==0) {
      scaleAcc = 0.5 #reduce sampling variation by 1/2 if none is accepted
    } else {
      scaleAcc = acc0/accRateGolden #obtain scale between accepted and Golden
    }
    delta = delta*scaleAcc #update delta
  }
  
  if(verbose) print("Sampling under Hp...")
  hpmcmc <- euroformix::contLikMCMC(mlefitHp,niter=niter,delta=delta,seed=seed)

  if(verbose) print("Sampling under Hd...")
  seed2 = seed+999 #same seed diff as in EFM
  if(length(seed2)==0) seed2 = NULL
  hdmcmc <- euroformix::contLikMCMC(mlefitHd,niter=niter,delta=delta,seed=seed2) 
  
  #Post-evaluation:
  log10LRdistr <- (hpmcmc$postlogL - hdmcmc$postlogL)/log(10) #calculate log10LR
  consLR  <- quantile(log10LRdistr,quantile)
  if(is.null(hpmcmc$logmargL)) { #old version of EFM (v3.1.0 and older) was used. MUST USE margL argument
    hpmcmc$logmargL = log(hpmcmc$margL)
    hdmcmc$logmargL = log(hdmcmc$margL)
  } 
  bayesLR <- (hpmcmc$logmargL-hdmcmc$logmargL)/log(10) #convert log to log10
  if(is.infinite(bayesLR)) bayesLR <- NaN #invalid number due to zero Hd

  mcmcList = list(bayesLR=bayesLR, consLR=consLR,
                  niter=niter,delta=delta,quantile=quantile,seed=seed,
                  accRateGolden=accRateGolden,diffTol=diffTol)
 return(mcmcList)
}