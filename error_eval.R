##########################
#                        #
#  Error evaluation      #
#                        #
##########################
source('scaling.R')



## for each index, compare approximation against ground truth
## (taking gcm output to be ground truth)
for (ndx in index_names)
{
  ## name of index
  ndx_objs <- ls(pattern = paste0('^',ndx,'.*scaled','.*RCP[0-9]{2}'))
  
  ## load value
  load(grep(ndx,ndxs,value = TRUE))
  
  
  for(index in ndx_objs)
  {
    
    scaled = eval(parse(text = index))
    
    
    ## get correct matrix for ground truth
    if (length(grep('_2_',index)) == 1) truth_name <- paste(ndx,'2pt0degC',sep = '.')
    else truth_name <- paste(ndx,'1pt5degC',sep = '.')
    
    
    ## evaluate string as object
    truth <- eval(parse(text = truth_name))
    print(paste('comparing',index,'against',truth_name))
    
    ## get century change in ground truth, because this is the approximation we made
    first_mean_truth = apply(truth[,,1:20,],c(1,2,4),mean, na.rm=TRUE)
    second_mean_truth = apply(truth[,,(dim(truth)[3]-20):dim(truth)[3],],c(1,2,4),mean,na.rm = TRUE)
    truth = second_mean_truth - first_mean_truth
    
    
    ## compute each component of error
    error = array(dim = c(288, 192, 4))
    scaled_mean <- apply(scaled,c(1,2),mean, na.rm=TRUE)
    scaled_var <- apply(scaled,c(1,2),var, na.rm=TRUE)
    truth_mean <- apply(truth,c(1,2),mean, na.rm=TRUE)
    truth_var <- apply(truth,c(1,2),var, na.rm=TRUE)
    meandiff = (truth_mean - scaled_mean)^2
    
    ## problems with model output, set NA
    meandiff[,1:32] <- NA
    scaled_mean[,1:32] <- NA
    scaled_var[,1:32] <- NA
    
    ## Replace outlier with 0.975th percentile
    meandiff_cutoff <- as.numeric(quantile(meandiff, 0.975,na.rm = TRUE))
    scaled_var_cutoff <- as.numeric(quantile(scaled_var, 0.975,na.rm = TRUE))
    truth_var_cutoff <- as.numeric(quantile(truth_var, 0.975,na.rm = TRUE))
    meandiff[which(meandiff > meandiff_cutoff,arr.ind = TRUE)] <- meandiff_cutoff
    scaled_var[which(scaled_var > scaled_var_cutoff,arr.ind = TRUE)] <- scaled_var_cutoff
    truth_var[which(truth_var > truth_var_cutoff,arr.ind = TRUE)] <- truth_var_cutoff
    
    
    ind = which(truth_var == 0,arr.ind = TRUE)
    
    error[,,1] <- meandiff
    error[,,2] <- scaled_var
    error[,,3] <- truth_var
    error[,,4] <- (scaled_var + meandiff) / truth_var
    
    
    ## save error object
    errname = paste(index,'error',sep = '_')
    assign(errname,error)
    
    #png(paste0(errname,'.png'))
    #image.plot(error)
    #dev.off()
    
    
  }
  to_rm = ls(pattern = paste0('^',ndx,'\\.'))
  rm(list = to_rm)
}



###########################################
#          Error Aggregation              #
###########################################


load('error_objects')

for (n in ls(pattern = '_error$'))
{
  e = eval(parse(text = n))
  ind = which(e[,,3] == 0,arr.ind = TRUE)
  
  if(nrow(ind) > 0)
  {
    print(n)
    e[ind[,1],ind[,2],c(1,2,3,4)] <- rep(NA, dim(ind)[1] * dim(ind)[2])
  }
  
  assign(n,e)
  
  
}


