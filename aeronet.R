################################################
## extract DRAGON data from website#############
################################################
setwd('~/projects/aerosol/')
source('src/funcs.R')
library(gdata)
library(plyr)
info = read.xls('src/MISR_INFO.xls')
location = 'Baltimore'
id = info$Location == location
info = info[id,c('Dates','Time')]
info_conv = t(apply(info,1,info_str))
colnames(info_conv) = c('Date','Time')
filenames = list.files(path = paste0('aeronet/raw/',location),full.names = T,pattern = '\\.lev20$')


for (i in 1:nrow(info_conv)){
  raw = get_raw(filenames,info_conv[i,1])
  stdtime = strptime(paste(info_conv[i,1],info_conv[i,2]),"%d:%m:%Y %H:%M:%S", tz='UTC')
  lower = -30
  upper = 30
  if (!is.na(stdtime)){
    data = get_data(raw,stdtime,lower,upper)
    if (nrow(data)!=0){
      
      output = ddply(data,.(Site_Name,lon,lat),summarize,AOT_446 = mean(AOT_446),AOT_558 = mean(AOT_558), AOT_672 = mean(AOT_672), AOT_867 = mean(AOT_867))
    
      out_file = paste(info$Dates[i],'aeronet.csv',sep = '_')
      out_path = paste('aeronet/processed2',location,sep='/')
      if (!dir.exists(out_path))
        dir.create(out_path,recursive = T)
      write.table(output,file = paste(out_path,out_file,sep='/'), sep=',', eol = "\n",row.names=F)
      cat(as.character(stdtime),'is done!\n')
    }
    else{
      cat('\n no data available on',as.character(stdtime),'\n')
    }
  }
  else{
    cat(stdtime,'time stamp error!\n')
  }  
}




