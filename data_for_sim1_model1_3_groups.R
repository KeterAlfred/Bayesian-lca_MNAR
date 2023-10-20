#===========================================================================
#
#
# Simultaneous imputation of the missing bacteriological test results in the 
# analysis model under the assumption that the missing bacteriological test 
# results are missing not at random (MNAR)
#===========================================================================



#Design Matrix

K = 5   #Number of tests
D5 = matrix(NA, nrow = 2^K, ncol = K+1)

D5[,1] = 1
D5[,2] = rep(c(1,0), 2^K/2)
D5[,3] = rep(c(rep(1,2), rep(0,2)),2^K/4)
D5[,4] = rep(c(rep(1,4), rep(0,4)),2^K/8)
D5[,5] = rep(c(rep(1,8), rep(0,8)),2^K/16)
D5[,6] = rep(c(rep(1,16), rep(0,16)),2^K/32)


#Questions
#What do you think is the prevalence of PTB in those not tested with Xpert Ultra & culture?
#What do you think is the sensitivity of Xpert Ultra and culture among those not tested using Xpert Ultra and culture?
#Should we assume perfect specificity for Xpert Ulta and culture? among the unverified?






sim1_dat <- foreach(m=c(1:M)+j, .inorder = FALSE) %do% {




#Data for the model


  dat.tested = as.data.frame(sim_dat_Final1[,,m])
  #First create a variable with 1s
  ones = rep(1, length.out = dim(dat.tested)[1])
  
  var.nm  <- c("Y1", "Y2", "Y4",  "Y50", "Y60","groups")
  
  
  
  dat.tested[,"Y50"] <- ifelse( dat.tested[,"Y50"] %in% 9, NA, dat.tested[,"Y50"])
  dat.tested[,"Y60"] <- ifelse( dat.tested[,"Y60"] %in% 9, NA, dat.tested[,"Y60"])
  
  dat.tested$ones = ones
  dat.tested$g23 = (dat.tested$groups == c(2,3))*1
  # dat.tested$g3 = (dat.tested$groups == 3)*1
  dat.tested$g4 = (dat.tested$groups == 4)*1
  jags.data = dat.tested[,c("ones",var.nm,c("g23","g4"))]
  
  jags.data$g1 = (jags.data$groups %in% c(2,3,4))*1
  
  
  #Prepare the data
  jags.sim.data_3_groups = list( N = nrow(jags.data), 
                                 y = jags.data[,c(1:(dim(jags.data)[2]))],
                                 J = dim(jags.data)[1],
                                 K = 5,
                                 G = 3,
                                 nv = c(nrow(jags.data[jags.data$groups==1,]),
                                        nrow(jags.data[jags.data$groups %in% c(2,3),]),
                                        nrow(jags.data[jags.data$groups==4,])),
                                 D = D5,
                                 dim_D = dim(D5)[1]
  )
  #Assign names to the objects in the list
  names(jags.sim.data_3_groups) <- c("N","y","J","K","G","nv","D","dim_D")
  
  
  
  
  
  return( jags.sim.data_3_groups )
}
