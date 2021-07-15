LHC_MR_betX_v9_LD <- function(theta,betX,pi1,sig1,weights,m0,nX,bn=2^8,bins=15){
  # return ( logL, dlogL, nmiss)
  
  M = 1e7;
  #M = 6372305;
  piX = abs(theta[1]);
  h2X = abs(theta[2]);
  iX = abs(theta[3]);

  sigX = sqrt(h2X/(piX*M))
  # Number of genotyped SNPs
  m = length(betX)
  
  Rp = iX/nX
  bX = array(betX, c(1,m)) # reshape(betXY(:,1),[1,1,m]);
  
  # Define grid for FFT
  minX = mean(bX)-(5*sd(bX));
  maxX = mean(bX)+(5*sd(bX));
  dX = (maxX-minX)/(bn-1);
  minX = minX-dX/2;
  maxX = maxX+dX/2;

  bXi = ceiling((bX-minX)/dX);

  #### change!! ####
  bXi[bXi<1] = 1;
  bXi[bXi>bn] = bn;
  #su = which(bXi[1,]>=1 & bXi[1,]<=bn)
  #bXi = bXi[su];
  # Define grid for rho
  #pi1 = pi1[su];
  #sig1 = sig1[su];
  
  ## update weights vector if some SNPs are excluded
  #my_weights = weights[su]
  #nmiss = length(su)
  
  #h2X = M*piX*sigX^2;
  #h2Y = M*piY*sigY^2;
  
  if(piX > 0.2 || piX < 1e-6 || h2X > 1 || h2X < 1e-6 || iX > 1.5 || iX < 0.5 ){
  #if( piX>0.2 || h2X>1 || iX<=0.8 ){
    #|| min(c(iX,iY))<=0.8 || max(c(iX,iY))>=1.5
    logL = 1e10
  }else{
    min_pi1 = min(pi1)-1e-10;
    max_pi1 = max(pi1)+1e-10;
    dp = (max_pi1-min_pi1)/bins;
    pc = min_pi1 + (dp * matrix(seq(0.5,(bins-0.5),1),ncol=1)) # min_pi1+dp*(.5:1:(bins-.5))';
    pix = ceiling((pi1-min_pi1)/dp);
    min_sig1 = min(sig1)-1e-10;
    max_sig1 = max(sig1)+1e-10;
    ds = (max_sig1-min_sig1)/bins;
    sc = min_sig1 + (ds * matrix(seq(0.5,(bins-0.5),1),ncol=1)) #min_sig1+ds*(.5:1:(bins-.5))';
    six = ceiling((sig1-min_sig1)/ds);
    cix = pix + bins*(six-1);
    uni_cix = sort(unique(cix))
    #[~,ucix,ixMap] = unique(cix);
    ucix = match(uni_cix, cix)
    ixMap = match(cix, uni_cix)
    Sig1 = sc[six[ucix]];
    Pi1 = pc[pix[ucix]];
    mm = length(Sig1);
    
    #sigK = aperm(array(rep(Sig1, bn*mm), c(mm,bn,bn)), c(2,3,1)) # reshape(Sig1,[1,1,mm]);
    #piK = aperm(array(rep(Pi1, bn*mm), c(mm,bn,bn)), c(2,3,1)) # reshape(Pi1,[1,1,mm]);
    Ax = aperm( array(rep((1/sigX) / Sig1, bn), c(mm,bn)), c(2,1))
    Qx = aperm( array(rep(Pi1 * piX, bn), dim = c(mm, bn)),c(2,1))
    

    j = 0:(bn-1) #';
    vi = 2*pi*(j-bn/2)/(maxX-minX-dX);
    #Rx = array(rep(vi, bn*mm), dim = c(bn, bn, mm)) # reshape(vi,[bn,1,1]);

    #coX = (Rx+alp*Ry)*sigK;
    Rx = array(rep(vi, mm), c(bn,mm))
    
    Lx = -m0 * ( 1 - 1 / sqrt(1 + Rx^2/Ax^2))*Qx
    #Lx = -m0*piX*( 1 - 1/sqrt(1+sigX^2*coX^2) )*piK;
    Le = -(1/2)*(Rp*Rx^2)
    #Le = -(1/2)*( (Rp[1,1]*Rx^2+Rp[2,2]*Ry^2) + 2*Rp[1,2]*(Rx*Ry) );
  
    # /!\ complex numbers here!
    mf_init = -2 * log(as.complex(-1)) * ( (minX+dX/2) / (maxX-minX-dX) )*j
    mf = array(rep(mf_init, mm), dim=c(bn, mm))
    
    phi = exp(Lx+Le+mf);
    
    # in R, not possible to chose dimension for FFT, so we need a loop to do it for all rho bins
    FFT=array(NA, dim=c(bn, mm))
    for(l in 1:mm){
      FFT[,l] = fft(phi[,l])
    }
    
    # FFTmod = (1/((maxX-minX-dX)*(maxY-minY-dY)))*(-1).^(bn*((minX+dX/2)/(maxX-minX-dX) + ...
    #           (minY+dY/2)/(maxY-minY-dY)) + (j+j'));
    #FFTmod_init = (1/((maxX-minX-dX)*(maxY-minY-dY)))*(as.complex(-1))^(bn*((minX+dX/2)/(maxX-minX-dX) +
    #                                                                          (minY+dY/2)/(maxY-minY-dY)) + (j%+%t(j)))
    FFTmod_init = (1/(maxX-minX-dX))*(as.complex(-1))^(bn*((minX+dX/2)/(maxX-minX-dX)) + j)
    FFTmod = array(rep(FFTmod_init, mm), dim=c(bn,mm))
    
    FFT0 = Re(FFT*FFTmod)
    pfE = FFT0[cbind(t(bXi), ixMap)];
    length(which(pfE<0))
    # if some, remove them & update weights to keep the correct set of SNPs
    my_weights = weights[pfE>0]
    pfE=pfE[pfE>0]

    # we use m * mean(...) to account for SNPs that may have been excluded before
    logL = -m * mean(log(pfE*my_weights))
  }
  return(logL)
}
