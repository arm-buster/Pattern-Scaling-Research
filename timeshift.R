## timeshift method


load('RData_glbtas')
load('RData_coo_gat')
annglbtas.RCP85 = allgat[76:181,,'RCP85']



timeshift_cutoff <- c()
count = 1



for(highscen in highscenarios)
{
  for(lowscen in lowscenarios)
  {
    
    
    i = 1
    
    hscen = eval(parse(text = paste('annglbtas',highscen,sep = '.')))
    lscen = eval(parse(text = paste('annglbtas',lowscen,sep = '.')))
    
    cutoff = mean(tail(lscen,20))
    
    while(mean(hscen[i:(i+20),]) < cutoff) i = i + 1
    
    timeshift_cutoff[[i]] <- c(highscen,lowscen,i, i + 20)
    
    
    count = count + 1
    
  }
}



tshift_cutoffs <- as.data.frame(do.call(rbind,timeshift_cutoff))
names(tshift_cutoffs) <- c('hscen','lscen','start_yr','end_yr')
tshift_cutoffs$start_yr <- as.numeric(as.character(tshift_cutoffs$start_yr))
tshift_cutoffs$end_yr <- as.numeric(as.character(tshift_cutoffs$end_yr))










