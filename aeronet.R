######################################
## extract AERONET data from Yang ####
######################################
raw1 = read.csv(file = 'DisAQ_AERONET_lev02_06012011_08232011_10AM.csv',header = T,sep=',')
raw2 = read.csv(file = 'DisAQ_AERONET_lev02_06012011_08242011_11AM.csv',header = T,sep=',')
raw = rbind(raw1,raw2)
get_data = function(raw,stdtime,doy,lower,upper){
  repo = data.frame()
  for (site in unique(raw$site_name)){
    new_data = raw[raw$site_name == site,]
    Hour = new_data$hour
    Minute = new_data$minute
    Time = paste(Hour,Minute)
    epoch = strptime(Time, "%H %M", tz = 'UTC')
    repo = rbind(repo,new_data[new_data$doy == doy & as.numeric(difftime(epoch,stdtime,units='mins')) < upper & as.numeric(difftime(epoch,stdtime,units='mins')) > lower,])
  }
  return (repo)
}

time_str = "15:57"
doy = 203
lower = -30
upper = 30
stdtime = strptime(time_str,"%H:%M", tz='UTC')
data = get_data(raw,stdtime,doy,lower,upper)
output = aggregate(data$AOT_555, list(site_name = data$site_name, lon = data$lon, lat = data$lat), mean)
colnames(output)[4] = 'AOT_555'

write.table(output,file='2011.07.22_aeronet.csv',sep=',',eol = "\n",row.names=F)

#sites = read.table(file = 'dragon_2011_locations_2011_lev20.txt', header = T, sep = ',', skip = 1)
#colnames(sites)[2:3] = c('Lon','Lat')
#output = merge(x = output, y = sites[,1:3], by = 'Site_Name', all.x = T)
#output = output[,c(1,3,4,2)]
################################################
## extract DRAGON data from website#############
################################################
rm(list=ls())
setwd('~/Graduate_Research/Aerosol_Retrieval/AERONET_DATA')
get_raw = function(filenames,date_str){
  repo = data.frame()
  for (filename in filenames){
    site = readLines(filename,4)
    raw = read.table(file = filename, header = T, sep=',', skip = 4,colClasses = "character")
    site_name = gsub(',','',gsub('Location=','',regmatches(site[3], regexpr('Location=(.+?),',site[3]))))
    lon = gsub(',','',gsub('long=','',regmatches(site[3], regexpr('long=(.+?),',site[3]))))
    lat = gsub(',','',gsub('lat=','',regmatches(site[3], regexpr('lat=(.+?),',site[3]))))
    id = raw$Date.dd.mm.yy. == date_str
    if (sum(id & raw$AOT_440!='N/A') > 0){
      new_data = raw[id & raw$AOT_440!='N/A', c('Date.dd.mm.yy.','Time.hh.mm.ss.','AOT_440','X440.675Angstrom')]
      colnames(new_data) = c('Date','Time','AOT_440','AE')
      new_data$Date = gsub(':','/',new_data$Date)
      new_data$AE = as.numeric(new_data$AE)
      new_data$AOT_440 = as.numeric(new_data$AOT_440)
      new_data$AOT_555 = new_data$AOT_440 * exp( - log(558/440) * new_data$AE )    
      new_data$Site_Name = site_name
      new_data$lon = as.numeric(lon)
      new_data$lat = as.numeric(lat)
      new_data = new_data[,c(6,7,8,1,2,3,4,5)]
      repo = rbind(repo,new_data)
    }
  } 
  return (repo)
}

get_data = function(raw,stdtime,lower,upper){
  repo = data.frame()
  for (site in unique(raw$Site_Name)){
    new_data = raw[raw$Site_Name == site,]
    epoch = strptime(paste(new_data$Date,new_data$Time), "%d/%m/%Y %H:%M:%S", tz = 'UTC')
    idx = as.numeric(difftime(epoch,stdtime,units = 'mins')) < upper & as.numeric(difftime(epoch,stdtime,units = 'mins')) > lower
    if (sum(idx)>0)
      repo = rbind(repo,new_data[idx,])
  }
  return (repo)
}

info_str = function(info){
  Date = info[1]
  Time = info[2]
  Year = regmatches(Date,regexpr('^[0-9]+', Date))
  Month = gsub('\\.','',regmatches(Date,regexpr('\\.[0-9]+\\.',Date)))
  Day = regmatches(Date,regexpr('[0-9]+$', Date))
  date_str = paste(Day,Month,Year,sep=':')
  time_str = gsub("\\'",'',Time)
  #time_str = Time
  info_str = c(date_str,time_str)
  return(info_str)
}

library(gdata)
info = read.xls('MISR_INFO.xls')
location = 'Beijing'
#id = info$Location == location
id = 62:75
info = info[id,c('Dates','Time')]
info_conv = t(apply(info,1,info_str))
colnames(info_conv) = c('Date','Time')
filenames = list.files(path = paste('./RAW',location,sep='/'),full.names = T,pattern = '\\.lev20$')

for (i in 1:nrow(info_conv)){
  raw = get_raw(filenames,info_conv[i,1])
  stdtime = strptime(paste(info_conv[i,1],info_conv[i,2]),"%d:%m:%Y %H:%M:%S", tz='UTC')
  lower = -30
  upper = 30
  if (!is.na(stdtime)){
    data = get_data(raw,stdtime,lower,upper)
    if (nrow(data)!=0){
      output = aggregate(data$AOT_555, list(Site_Name = data$Site_Name,Lon = data$lon, Lat = data$lat), mean)
      colnames(output)[4] = 'AOT_555'
      out_file = paste(info$Dates[i],'aeronet.csv',sep = '_')
      out_path = paste('./PROCESSED',location,sep='/')
      dir.create(out_path,showWarnings = F)
      write.table(output,file = paste(out_path,out_file,sep='/'), sep=',', eol = "\n",row.names=F)
    }
    else{
      cat('\n no data available on',as.character(stdtime))
    }
  }
  else{
    cat(stdtime,'time stamp error!\n')
  }  
}

