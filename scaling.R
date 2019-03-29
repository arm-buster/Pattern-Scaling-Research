## pattern scaling
## run timeshift.r
library(stringi)
library(dplyr)



setwd('~/Desktop/PatternScaling/xtreme_indices')
ndxs <- list.files(pattern = 'RData_allscenarios')
index_names = c('cdd','fd','gsl','r10mm','r95ptot','rx5day',
                'sdii','tnn','txx','wsdi')


## scenarios
scenarios = c('1pt5degC','2pt0degC','RCP45','RCP85')
highscenarios = scenarios[3:4]
lowscenarios = scenarios[1:2]

source('timeshift.R')

## all global temp objs
glbtemps <- ls(pattern = "^annglbtas*")




## take average global temperature for begining and end of century
temp_avgs <- c()
for(g in scenarios){
  
  scen = paste('annglbtas',g,sep = '.')
  
  glbtemp = eval(parse(text = scen))
  change <- apply(tail(glbtemp,20),2,mean) - apply(glbtemp[1:20,],2,mean)
  temp_avgs <- rbind(temp_avgs,c(g,change))
  
}

########################################################
#                                                      #
#          pattern scaling loop                        #
#          each iteration we approximate               #
#          1.5C and 2.0 values for one index           #
#                                                      #
########################################################
i = 1
for (ndx in ndxs)
{ 
  
  print(paste('pattern scaling:',ndx))
  load(ndx)
  
  for (scen in highscenarios)
  {
    print(paste(index_names[i],scen,sep = '.'))
    x = paste(index_names[i],scen,sep = '.')
    
    
    mx = eval(parse(text = x))
    
    
    ## mean at each gridpoint of first 20 years
    first_mean = apply(mx[,,1:20,],c(1,2,4),mean, na.rm=TRUE)
    
    ## mean at each gridpoint of last 20 years
    second_mean = apply(mx[,,(dim(mx)[3]-20):dim(mx)[3],],c(1,2,4),mean,na.rm = TRUE)
    
    ## century change
    change = second_mean - first_mean
    
    ## gat
    gat = as.numeric(temp_avgs[temp_avgs[,1] == scen,][2:11])
    
    ## grid (pattern) of local change per degree of warming
    pattern = aperm(change,c(3,1,2)) / rep(gat,288*192)
    
    
    ## get GAT for the scenarios
    gat1.5 = as.numeric(temp_avgs[temp_avgs[,1] == "1pt5degC",][2:11])
    gat2 = as.numeric(temp_avgs[temp_avgs[,1] == "2pt0degC",][2:11])
    
    
    ## now scale the pattern to the scenario we want to approximate
    scaled1.5 = pattern * rep(gat1.5,288*192)
    scaled2 = pattern * rep(gat2,288*192)
    
    scaled1.5 = round(aperm(scaled1.5,c(2,3,1)),digits = 2)
    scaled2 = round(aperm(scaled2,c(2,3,1)),digits = 2)
    
    cutoffs1.5 = as.numeric(filter(tshift_cutoffs, hscen == scen, lscen == '1pt5degC')[3:4])
    cutoffs2 = as.numeric(filter(tshift_cutoffs, hscen == scen, lscen == '2pt0degC')[3:4])
    
    
    scaled1.5_tshift = apply(mx[,,cutoffs1.5[1]:cutoffs1.5[2],],c(1,2,4),mean,na.rm = TRUE) - apply(mx[,,1:20,],c(1,2,4),mean,na.rm = TRUE)
    scaled2_tshift = apply(mx[,,cutoffs2[1]:cutoffs2[2],],c(1,2,4),mean,na.rm = TRUE) - apply(mx[,,1:20,],c(1,2,4),mean,na.rm = TRUE)
    
    
    
    
    assign(paste(index_names[i],'scaled_1.5','from',scen,sep = '_'),scaled1.5)
    assign(paste(index_names[i],'scaled_2','from',scen,sep = '_'),scaled2)
    
    assign(paste(index_names[i],'timeshift_scaled_1.5','from',scen,sep = '_'),scaled1.5_tshift)
    assign(paste(index_names[i],'timeshift_scaled_2','from',scen,sep = '_'),scaled2_tshift)
    
    
    
    
  }
  
  to_rm = ls(pattern = paste0('^',index_names[i],'\\.'))
  rm(list = to_rm)
  
  i = i + 1
}
