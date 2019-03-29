library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(ggrepel)

##############################
#### barplot scaled error ####
##############################

#final = final[!grepl('wsdi', final$scaling),]



plots = letters[1:8]
i=1
for (scale_method in c('Linear Rescaling','Timeshift'))
{
  for(hscen in c('RCP4.5','RCP8.5'))
  {
    for(lscen in c('1.5C', '2.0C'))
    {
      k <- filter(final, starting_scen == hscen & target_scen == lscen & method == scale_method)
      
      ## should this be k$meandiff^2 / k$truth_var ??
      ## no, see line 157, this is (scaledmean - truthmean)^2
      
      k$emulator_error = k$meandiff / k$truth_var
      k$scaled_scaled_var = k$scaled_var / k$truth_var
      k = select(k, scaling,emulator_error,scaled_scaled_var)
      k$true_var <- 1
      
      k_long <- melt(k, id.vars = 'scaling')
      
      
      k_long$index = as.character(as.matrix(stri_match(k_long$scaling, regex = '^[^_]+(?=_)')))
      names(k_long)[4] <- 'index'
      
      total_error = summarise(group_by(k_long,index), tot_error = sum(value))
      indx_orders = total_error$index[order(total_error$tot_error)]
      
      k_long$index = factor(k_long$index, levels = indx_orders)
      #levels(k_long$index) = indx_orders
      #levels(k_long$index) = filter(k_long, variable == 'emulator_error')$index[indx_orders]
      
      t = paste(scale_method,'of',hscen,'to',lscen)
      
      print(t)
      
      png(paste0(t,'.png'))
      
      #y_scale = 0:ceiling(max(rowSums(k[2:4])))
      
      
      
      if(i != 1){assign(plots[i],ggplot(k_long) + geom_col(aes(x = index, y = value, fill = variable)) +
                          scale_y_continuous(limits = c(0,3.2)) + guides(fill=FALSE) +
                          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                          geom_hline(yintercept = 2,color = 'orange',linetype = 2, size = 1.4) + 
                          scale_fill_manual(labels = c('Squared Mean Difference \n (Systematic Error)','Emulation Variance','Truth Variance'),values = c('red','green','blue')) +
                          labs(y = 'Total Error', title = t))}
      
      else {assign(plots[i],ggplot(k_long) + geom_col(aes(x = index, y = value, fill = variable)) +
                     scale_y_continuous(limits = c(0,3.2)) +
                     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                     geom_hline(yintercept = 2,color = 'orange',linetype = 2, size = 1.4) +
                     scale_fill_manual(labels = c('Squared Mean Difference \n (Systematic Error)','Emulation Variance','Truth Variance'),values = c('red','green','blue')) +
                     labs(y = 'Total Error', title = t))}
      
      
      
      
      print(eval(parse(text = plots[i])))
      
      dev.off()
      i = i + 1
    }
  }
}


grid.arrange(e,b,c,d,nrow = 2,top = 'Error Components Normalized by Ground Truth Variance')

###################################
## error against GAT change #######
###################################
load('RData_glbtas')
load('RData_coo_gat')
annglbtas.RCP85 = allgat[76:181,,'RCP85']
scenarios = c('1pt5degC','2pt0degC','RCP45','RCP85')





temp_avgs <- c()
for(g in scenarios){
  
  scen = paste('annglbtas',g,sep = '.')
  
  glbtemp = eval(parse(text = scen))
  change <- apply(tail(glbtemp,20),2,mean)
  temp_avgs <- rbind(temp_avgs,c(g,change))
  
}



temps <- apply(apply(temp_avgs[,2:11], c(1,2), as.character), c(1,2), as.numeric)
gt_change <- as.data.frame(cbind(temp_avgs[,1],rowMeans(temps)))

names(gt_change) <- c('scen','gat')
gt_change$gat <- as.numeric(as.character(gt_change$gat))
gt_change$scen = c('1.5C','2.0C','RCP4.5','RCP8.5')




v <- c()
i = 1
for(hscen in c('RCP4.5','RCP8.5'))
{
  for(lscen in c('1.5C','2.0C'))
  {
    v[[i]] <- c(hscen,lscen, gt_change[gt_change[,1] == hscen,]$gat - gt_change[gt_change[,1] == lscen,]$gat)
    i = i + 1
  }
  
}

tchanges <- as.data.frame(do.call(rbind, v))
names(tchanges) <- c('starting_scen','target_scen','distance')

final <- left_join(final, tchanges, by = c('starting_scen','target_scen'))
final$distance <- as.numeric(as.character(final$distance))


x_breaks = unique(final$distance)[c(3,1,4,2)]
x_labs = c('0.32 \n (RCP4.5 to 2.0C)',
           '0.85 \n (RCP4.5 to 1.5C)',
           '2.54 \n (RCP8.5 to 2.0C)',
           '3.07 \n (RCP8.5 to 1.5C)')


outlier_lab <- final[grepl('wsdi',final$scaling),]
outlier_lab <- outlier_lab[grepl('Linear Rescaling',outlier_lab$method),]
outlier_lab$label = 'WSDI'


a <- ggplot(final) + geom_point(aes(distance, error, color = method)) +
  scale_x_continuous(breaks = unique(final$distance))+
  geom_hline(yintercept = 1,color = 'orange',linetype = 2, size = 1.4) +
  scale_x_continuous(breaks = x_breaks, labels = x_labs)+
  geom_text_repel(data = outlier_lab, aes(x = distance,y = error,label = label)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = NULL)





nonwsdi <- final[!grepl('wsdi',final$scaling),]



errs <- summarise(group_by(nonwsdi, target_scen, starting_scen), error = mean(error))
nonwsdi_ = select(nonwsdi, 'target_scen','starting_scen','distance')
nonwsdi_ = nonwsdi_[!duplicated(nonwsdi_),]

errs = left_join(errs,nonwsdi_, by = c('target_scen','starting_scen'))



tshift <- nonwsdi[grepl('Timeshift',nonwsdi$method),]
linear <- nonwsdi[grepl('Linear Rescaling',nonwsdi$method),]







## get rid of legend on right side, tilt labels 45 degrees

b <- ggplot(nonwsdi) + geom_point(aes(distance, error, color = method)) +
  scale_x_continuous(breaks = unique(final$distance)) + guides(color=FALSE) +
  geom_hline(yintercept = 1,color = 'orange',linetype = 2, size = 1.4) +
  scale_x_continuous(breaks = x_breaks, labels = x_labs)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = 'WSDI removed',
       x = 'Difference Between End-Of-Century GAT Average \n of Starting and Target Scenarios')





grid.arrange(a,b,ncol = 1,top = 'Error Against Difference Between Starting and Target Scenarios')

