get_raw = function(filenames,date_str){
  repo = data.frame()
  cid = c('Date.dd.mm.yy.','Time.hh.mm.ss.','AOT_440','AOT_675','AOT_870','X440.675Angstrom','X500.870Angstrom','X440.870Angstrom')
  for (filename in filenames){
    site = readLines(filename,4)
    raw = read.table(file = filename, header = T, sep=',', skip = 4,colClasses = "character")
    site_name = gsub(',','',gsub('Location=','',regmatches(site[3], regexpr('Location=(.+?),',site[3]))))
    lon = gsub(',','',gsub('long=','',regmatches(site[3], regexpr('long=(.+?),',site[3]))))
    lat = gsub(',','',gsub('lat=','',regmatches(site[3], regexpr('lat=(.+?),',site[3]))))
    rid = raw$Date.dd.mm.yy. == date_str & (raw$AOT_440!='N/A' | raw$AOT_500!='N/A')
    num_valid = sum(rid)
    
    if ( num_valid > 0){
      new_data = raw[rid,cid]
      colnames(new_data) = c('Date','Time','AOT_440','AOT_675','AOT_870','AE_440_675','AE_500_870','AE_440_870')
      
      new_data$Date = gsub(':','/',new_data$Date)
      new_data$AE_440_675 = as.numeric(new_data$AE_440_675)
      new_data$AOT_440 = as.numeric(new_data$AOT_440)
      new_data$AOT_675 = as.numeric(new_data$AOT_675)
      
      new_data$AOT_446 = new_data$AOT_440 * exp( - log(446/440) * new_data$AE_440_675 )
      new_data$AOT_558 = 0.5 * ( new_data$AOT_440 * exp( - log(558/440) * new_data$AE_440_675 ) + new_data$AOT_675 * exp( - log(558/675) * new_data$AE_440_675 ) )
      new_data$AOT_672 = new_data$AOT_675 * exp( - log(672/675) * new_data$AE_440_675 )

      if (new_data$AOT_870[1] == 'N/A'){
        new_data$AE_440_870 = as.numeric(new_data$AE_440_870)
        new_data$AOT_867 = new_data$AOT_675 * exp( - log(867/675) * new_data$AE_440_870 )
        new_data$AE_440_870 = NULL
        new_data$AE_500_870 = NULL
        new_data$AOT_870 = NULL
      }else{
        new_data$AE_500_870 = as.numeric(new_data$AE_500_870)
        new_data$AOT_870 = as.numeric(new_data$AOT_870)
        new_data$AOT_867 = new_data$AOT_870 * exp( - log(867/870) * new_data$AE_500_870 )
        new_data$AE_440_870 = NULL
        new_data$AE_500_870 = NULL
        new_data$AOT_870 = NULL
      }
      
      
      tmp = data.frame(Site_Name=rep(site_name,num_valid),lon=rep(as.numeric(lon),num_valid),lat=rep(as.numeric(lat),num_valid))
      
      new_data = cbind(tmp,new_data)
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
  info_str = c(date_str,time_str)
  return(info_str)
}