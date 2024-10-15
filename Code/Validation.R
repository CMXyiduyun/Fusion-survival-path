

if(0){
  # Patients' data must be formatted to match the structure of our mock data. 
  # Once formatted, the data can be run through the validation part of our code, 
  # which will return specific nodes and an estimated overall survival length. 
}


rm(list=ls())


# import libraries
library(readxl)
library(dplyr)  


# set parameters
max_ts <- 8


# set root
setwd('../Data/')


# load data
Fusion_SP <- read_excel('Fusion-SP result.xlsx')
validation <- read_excel('Mock data.xlsx')
validation <- validation[validation$`time slice` <= max_ts, ]


# assign node
for(ts in 1:max_ts){
  
  # locate time slice-specific data
  ts.dev <- Fusion_SP[Fusion_SP$`time slice`==ts, ]
  ts.val <- validation[validation$`time slice`==ts, ]
  ts.val[[paste('TS ', ts, ' merge', sep='')]] <- NA
  ts.val[[paste('TS ', ts, ' est OS', sep='')]] <- NA
  
  
  # identify parent node
  if(ts==1){
    parent.node <- 'BCLC stage'
  }else{
    parent.node <- paste('TS ', (ts-1), ' merge', sep='')
  }
  
  
  # assign child node and corresponding estimated overall survival length
  for(i in 1:nrow(ts.dev)){
    
    pn <- ts.dev[['Parent node']][i]
    var <- ts.dev[['Variable to split']][i]
    cn0 <- ts.dev[['Child node upper']][i]
    estOS0 <- round(ts.dev[['Child node upper median OS']][i], 1)
    cn1 <- ts.dev[['Child node lower']][i]
    
    if(var=='None'){
      ts.val[ (!is.na(ts.val[[parent.node]])) & 
              (ts.val[[parent.node]]==pn), 
      ][[paste('TS ', ts, ' merge', sep='')]]  <- cn0
      
      ts.val[ (!is.na(ts.val[[parent.node]])) & 
                (ts.val[[parent.node]]==pn), 
      ][[paste('TS ', ts, ' est OS', sep='')]]  <- estOS0
      
    }else{
      # child-upper
      estOS1 <- round(as.numeric(ts.dev[['Child node lower median OS']][i]), 1)
      
      ts.val[ (!is.na(ts.val[[parent.node]])) & 
              (ts.val[[parent.node]]==pn) & (ts.val[[var]]==0), 
            ][[paste('TS ', ts, ' merge', sep='')]]  <- cn0
      
      ts.val[ (!is.na(ts.val[[parent.node]])) & 
                (ts.val[[parent.node]]==pn) & (ts.val[[var]]==0), 
      ][[paste('TS ', ts, ' est OS', sep='')]]  <- estOS0
    
      # child-lower
      ts.val[ (!is.na(ts.val[[parent.node]])) & 
                (ts.val[[parent.node]]==pn) & (ts.val[[var]]==1), 
      ][[paste('TS ', ts, ' merge', sep='')]]  <- cn1
      
      ts.val[ (!is.na(ts.val[[parent.node]])) & 
                (ts.val[[parent.node]]==pn) & (ts.val[[var]]==1), 
      ][[paste('TS ', ts, ' est OS', sep='')]]  <- estOS1
    }
  
  }
  
  # merge to original validation data
  validation <- left_join(validation, 
                          ts.val[, c('ID',
                                     paste('TS ', ts, ' merge', sep=''),
                                     paste('TS ', ts, ' est OS', sep=''))], 
                          by='ID')
  
}

View(validation)


