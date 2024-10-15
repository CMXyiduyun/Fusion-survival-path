


rm(list=ls())


# import libraries
library(readxl)
library(dplyr)  
library(survminer)
library(survival)
library(Hmisc)
library(SurvivalPath)
library(survivalROC)
library(party)
library(randomForestSRC)
library(dbscan)



# set parameters 

merge <- 1
if(merge){
  p.cutoff <- 0.05
}else{
  p.cutoff <- 1
}

max_ts <- 2
min_samples <- 3
this_eps <- 1
min_median <- 0
  


# define function

##  function: split parent node to child node(s)
result <- function(ID_name, ts){
  
  # ratio 5%~95% AND significant in uni-variate analysis
  significant_vars <- c()
  for (var in variables) {
    sum1.ratio <- sum(ID_name[[var]])/dim(ID_name)[1]*100
    
    if( (sum1.ratio>5) & (sum1.ratio<95) ){
      cox_fit <- coxph(Surv(`OS time (Month)`, `Status of death`) ~ .,
                       data = ID_name[, c('OS time (Month)', 'Status of death', var)])
      if(summary(cox_fit)$coefficients[5] < 0.05){
        if(summary(cox_fit)$coefficients[5] > 0){
          significant_vars <- c(significant_vars, var)
        }
      }
    }
  }
  selected.varibels <- significant_vars
  
  
  if(length(selected.varibels)==0){
    ID_new_name <- ID_name[, c('ID')]
    ID_new_name[[paste('ts_', ts, '_var', sep='')]] <- 'None'
    ID_new_name[[paste('ts_', ts, '_value', sep='')]] <- 'all'
    ID_new_name[[paste('Class_', ts, sep='')]] <- 1
  }else{
    # multi-variate analysis
    set.seed(20231211) 
    rfsrc_fit <- rfsrc(Surv(`OS time (Month)`, `Status of death`) ~ ., 
                       data = ID_name[, c('OS time (Month)', 'Status of death', selected.varibels)],
                       ntree = 500,  
                       nodedepth = 10,
                       importance=TRUE,
                       seed=20231211)
    
    rfsrc.FI <- data.frame(rfsrc_fit$importance)
    rfsrc.FI$Features <- selected.varibels 
    rfsrc.FI <- rfsrc.FI[order(rfsrc.FI$rfsrc_fit.importance, decreasing=T), ]
    VIP.most <- rfsrc.FI$Features[1]
    
    ID_name[[paste('ts_', ts, '_var', sep='')]] <- VIP.most
    ID_name[[paste('ts_', ts, '_value', sep='')]] <- ID_name[[VIP.most]]
    ID_name[[paste('Class_', ts, sep='')]] <- ID_name[[VIP.most]] + 1
    
    ID_new_name <- ID_name[, c('ID', 
                               paste('ts_', ts, '_var', sep=''),
                               paste('ts_', ts, '_value', sep=''),
                               paste('Class_', ts, sep=''))]
  }
  
  
  ID_new_name
}


## function: to determine whether differing in OS 
KM.p <- function(TS.res.bind, nodes){
  
  tmp <- TS.res.bind[TS.res.bind$tmp.node %in% nodes, ]
  
  fit <- survfit(Surv(`OS time (Month)`, `Status of death`) ~ tmp.node, data=tmp)
  ms <- surv_median(fit)$median
  
  msl <- surv_median(fit)$lower
  msu <- surv_median(fit)$upper
  
  na_pos <- which(is.na(ms))
  ms[na_pos] <- msl[na_pos]
  
  na_pos <- which(is.na(ms))
  ms[na_pos] <- msu[na_pos]
  
  na_pos <- which(is.na(ms))
  ms[na_pos] <- mean(c(ms[na_pos-1], ms[na_pos+1]), na.rm = TRUE) 
  
  sdiff <- survdiff(Surv(`OS time (Month)`, `Status of death`) ~ tmp.node, data=tmp)
  p.val <- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  
  p.val
  
}



## function: to determine whether differing in BCLC distribution
define.group <- function(TS.res.bind, nodes, ts){
  
  cnt <- 1
  for(n in nodes){
    tmp <- TS.res.bind[TS.res.bind$tmp.node==n, ]
    a <- tmp[tmp[, c(paste('BCLC stage', ts, sep=''))]=='0/A', ]
    b <- tmp[tmp[, c(paste('BCLC stage', ts, sep=''))]=='B', ]
    c <- tmp[tmp[, c(paste('BCLC stage', ts, sep=''))]=='C', ]
    abc <- c(n, tmp[[paste('Median_ts_', ts, sep='')]][1],
             dim(a)[1] / dim(tmp)[1] * 100, 
             dim(b)[1] / dim(tmp)[1] * 100, 
             dim(c)[1] / dim(tmp)[1] * 100)
    
    if(cnt==1){ABC <- abc}else{ABC <- rbind(ABC, abc)}
    cnt <- cnt + 1
  }
  
  colnames(ABC) <- c('node', 'Median', 'A', 'B', 'C')
  rownames(ABC) <- nodes
  ABC <- data.frame(ABC)
  
  res <- dbscan(ABC[, c('A', 'B','C')], eps=this_eps, minPts=1)
  ABC$group <- res[["cluster"]]
  
  
  ABC
}


# set root
setwd('../Data/')



# load data
dataset <- read_excel('Mock data.xlsx')
dataset <- dataset[dataset$`time slice` <= max_ts, ]

variables <- c(
  'Viable lesion',	
  'New lesion',	
  'Distant metastasis',	
  'Local lymph node metastasis',	
  'Varicosity',	
  'Vascular invasion',	
  'Number of hepatic lesions > 1',
  'Number of hepatic lesions > 3',	
  'Diameter of main lesion > 50 mm',	
  'Diameter of main lesion > 70 mm',
  'ALBI grade = 2/3',	
  'ALBI grade = 3',	
  'Child-Pugh score = Class B/C',	
  'Child-Pugh score = Class C',	
  'Serum AFP level > 200 IU/ml',
  'Serum AFP level > 400 IU/ml', 
  'Serum AFP level > 800 IU/ml'	
)


dataset <- dataset[, c('ID', 'Status of death', 'OS time (Month)',  'time slice', 'BCLC stage', variables)] 

for (i in 1:max_ts) {
  ts.BCLC <- dataset[dataset$`time slice` ==i, ][,c('ID', 'BCLC stage')]
  colnames(ts.BCLC) <- c('ID', paste('BCLC stage', i, sep=""))
  dataset <- left_join(dataset, ts.BCLC, by='ID')
}






# train
Node.number <- 1

for(ts in 1:max_ts){
  
  cat('\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
  print(paste('Time_slice=', ts, sep=''))
  
  
  # parent 
  if(ts==1){
    parent.nodes <- c('0/A', 'B', 'C', 'D')
    dataset$PN <- dataset$`BCLC stage1`
  }else{
    dataset$PN <- dataset[[paste('Class_', (ts-1), '_merge', sep='')]]
  }
  
  
  # Split: pdata.0/A, pdata.B...
  for (i in parent.nodes){
    
    ID_name <-  dataset[ (dataset$`time slice`==ts) & (!is.na(dataset$PN)) & (dataset$PN==i), ]
    pdata_name <- paste('pdata.', i, sep='')
    
    if(dim(ID_name)[1]==0){
      assign(pdata_name, ID_name)
    }else{
      if(dim(ID_name)[1]>min_samples){
        ID_new_name <- result(ID_name, ts)
        assign(pdata_name, dplyr::left_join(ID_name, ID_new_name, by='ID'))
      }else{
        pseudo.pdata <- ID_name
        pseudo.pdata[[paste('ts_', ts, '_var', sep='')]] <- 'None'
        pseudo.pdata[[paste('ts_', ts, '_value', sep='')]] <- 'all'
        pseudo.pdata[[paste('Class_', ts, sep='')]] <- 1
        assign(pdata_name, pseudo.pdata)
      }
    }
    
  }
  
  
  # obtain info of child node(s): TS1.res
  TS.res.name <- paste('TS', ts, '.res', sep='')
  assign(TS.res.name, list())
  all.m <- c()
  for(i in parent.nodes){
    pdata.name <- paste('pdata.', i, sep='')
    pdata <- get(pdata.name)
    if(dim(pdata)[1]==0){
      next
    }else{
      nClass <- max(pdata[[paste('Class_', ts, sep='')]])
      for(cl in 1:nClass){
        lst <- pdata[pdata[, c(paste('Class_', ts, sep=''))]==cl, ]
        
        fit <- survfit(Surv(`OS time (Month)`, `Status of death`) ~ 1, data=lst)
        m = surv_median(fit)$median[1]
        if(is.na(m)){m = surv_median(fit)$lower}
        if(is.na(m)){m = surv_median(fit)$upper}
        lst[, c(paste('Median_ts_', ts, sep=''))] = m
        all.m <- append(all.m, c(m))
        
        assign(TS.res.name, append(get(TS.res.name), list(lst)))
      }
    }
    
  }

  
  
  # place info of child node(s) in order by estimated OS: TS1.res.order
  TS.res.order.name <- paste('TS', ts, '.res.order', sep='')
  assign(TS.res.order.name, list())
  median.order <- order(all.m, decreasing = 1)     
  TS.res <- get(TS.res.name)
  for(i in median.order){
    assign(TS.res.order.name, 
           append(get(TS.res.order.name), list(TS.res[[i]]))
    )
  }

  
  
  # bind info of child node(s) together: TS1.res.bind
  TS.res.bind.name <- paste('TS', ts, '.res.bind', sep='')
  TS.res.order <- get(TS.res.order.name)
  for(i in 1:length(TS.res.order)){
    to.bind <- TS.res.order[[i]][, c('ID', 
                                     'OS time (Month)', 'Status of death', 
                                     paste('BCLC stage', ts, sep=''),
                                     paste('ts_', ts, '_var', sep=''),
                                     paste('ts_', ts, '_value', sep=''),
                                     paste('Class_', ts, sep=''),
                                     paste('Median_ts_', ts, sep=''))]
    to.bind$tmp.node <- i
    
    if(i==1){
      assign(TS.res.bind.name, to.bind)
    }else{
      assign(TS.res.bind.name, 
             rbind(get(TS.res.bind.name), to.bind)
      )
    }
    
  }
  
  TS.res.bind.original.name <- paste('TS', ts, '.res.bind.original', sep='')
  assign(TS.res.bind.original.name, get(TS.res.bind.name))
  
  
  # determine: whether merge child nodes: renew TS1.res.bind
  TS.res.bind <- get(TS.res.bind.name)
  child.nodes <- max(TS.res.bind$tmp.node)
  ABCs = define.group(TS.res.bind, c(1:child.nodes), ts)

  all.possible.group <- c(1:max(ABCs$group))
  for(g in all.possible.group){
    group <- ABCs$node[ABCs$group==g]
    if(length(group)>1){
      
      
      for(i in 1:(length(group)-1)){
        for(j in length(group):(i+1)){
          
          if( ( group[i] %in% c(unique(TS.res.bind$tmp.node)) ) &  
              ( group[j] %in% c(unique(TS.res.bind$tmp.node)) )
            ){
            
              nodes <- group[c(i:j)]
              judge1 <- KM.p(TS.res.bind, nodes)
              
              
              if(judge1 > p.cutoff){
                TS.res.bind$tmp.node[TS.res.bind$tmp.node %in% nodes] <- paste(nodes, collapse = " - ")
              }
            
          }else{next}
        }
      }
      
      
    }else{
      next
    }
    
    
  }
  assign(TS.res.bind.name, TS.res.bind)
  
  print(ABCs)
  print(table(TS.res.bind$tmp.node))
  
  
  
  # obtain info of merged child node(s):list：TS1.res.merge
  TS.res.merge.name <- paste('TS', ts, '.res.merge', sep='')
  assign(TS.res.merge.name, list())
  all.m.merge <- c()
  for(i in unique(TS.res.bind$tmp.node)){
    lst <- TS.res.bind[TS.res.bind$tmp.node==i, ]
    fit <- survfit(Surv(`OS time (Month)`, `Status of death`) ~ 1, data=lst)
    m = surv_median(fit)$median
    if(is.na(m)){m = surv_median(fit)$lower}
    if(is.na(m)){m = surv_median(fit)$upper}
    lst[, c(paste('Median_ts_', ts, sep=''))] = m
    all.m.merge <- append(all.m.merge, c(m))
    assign(TS.res.merge.name, append(get(TS.res.merge.name), list(lst)))
  }


  # place info of merged child node(s) in order by estimated OS: TS1.res.merge.order
  TS.res.merge.order.name <- paste('TS', ts, '.res.merge.order', sep='')
  assign(TS.res.merge.order.name, list())
  median.merge.order <- order(all.m.merge, decreasing = 1)
  TS.res.merge <- get(TS.res.merge.name)
  for(i in median.merge.order){
    assign(TS.res.merge.order.name, 
           append(get(TS.res.merge.order.name), list(TS.res.merge[[i]]))
    )
  }

  
  # bind info of merged child node(s) together：TS3.Merged
  TS.Merge.name <- paste('TS', ts, '.Merged', sep='')
  TS.res.merge.order <- get(TS.res.merge.order.name)
  for(i in 1:length(TS.res.merge.order)){
    to.bind <- TS.res.merge.order[[i]][, c('ID',
                                           paste('ts_', ts, '_var', sep=''),
                                           paste('ts_', ts, '_value', sep=''),
                                           paste('Median_ts_', ts, sep=''))]
    to.bind[[paste('Class_',ts,'_merge',sep='')]] <- paste('Node', Node.number, sep='')
    
    if(i==1){
      assign(TS.Merge.name, to.bind)
    }else{
      assign(TS.Merge.name, 
             rbind(get(TS.Merge.name), to.bind)
      )
    }
    
    Node.number <- Node.number + 1
  }


  parent.nodes <- unique(get(TS.Merge.name)[[paste('Class_',ts,'_merge',sep='')]])
  print(parent.nodes)
  
  dataset <- dplyr::left_join(dataset, get(TS.Merge.name), by='ID')

  
}
  

View(dataset)






