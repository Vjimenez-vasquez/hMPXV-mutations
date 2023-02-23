#generate a subset metadata with the following information#
# column 1 : ID #
# column 2 : epicodes #
# column 3 : province #
# column 4 : country #
# column 5 : lineage #
# column 6 : dates (Collection.date) #
# column 7 : $substitutions #
# column 8 : $aaSubstitutions #

library(lubridate)
library(plotly)
library(dplyr)
library(epical)

#install.packages("remotes")#
#remotes::install_github("chrismerkord/epical")#
#install.packages("plotly")#
#install.packages("dplyr")#
#install.packages("lubridate")#

#examples#
#mutations_time_nuc(data=e,xmin="2022-06-15",xmax="2022-09-15",ymin=0,ymax=8,freq_sup=11,freq_inf=2,title="Monkey_Peru",gene="nuc")#
#mutations_time_aa(data=e,xmin="2022-06-15",xmax="2022-09-15",ymin=0,ymax=8,freq_sup=11,freq_inf=2,title="Monkey : aasubstitutions",gene="prot")#

##################################
## CODE 1 : mutations_time_nuc ##
#################################

mutations_time_nuc <- function(data,xmin,xmax,ymin,ymax,freq_sup,freq_inf,title,gene){
  data = data 
  xmin = xmin 
  xmax = xmax
  ymin = ymin
  ymax = ymax
  freq_sup = freq_sup
  freq_inf = freq_inf
  title = title
  gene = gene
  region = paste0(gene,"_")
  
  #add new column dates#
  names(data) <- c("ID","epicodes","province","country","lineage","Collection.date","substitutions","aaSubstitutions")
  data <- add_epi_week(data,"Collection.date", system="cdc")
  data$week_date <- epi_week_date(data$epi_week,data$epi_year,system="cdc")
  data <- as.data.frame(data)
  data$dec_date <- decimal_date(ymd(data$week_date))
  data.frame(names(data))
  head(data)
  data.frame(names(data))
  
  #extract unique mutation lists#
  #separate by commas#
  elements <- unlist(strsplit(data$substitutions, ","))
  a <- sort(unique(elements))
  #unique(gsub(":.*","",a))#
  
  #generate table#
  
  c <- 0
  d <- 0 
  e <- 0 
  f <- 0 
  g <- 0 
  h <- 0 
  j <- 0
  for (i in 1:length(a)){
    c <- append(c,data[grep(a[i],data$substitutions),1])
    d <- append(d,data[grep(a[i],data$substitutions),2])
    e <- append(e,data[grep(a[i],data$substitutions),3])
    f <- append(f,data[grep(a[i],data$substitutions),4])
    g <- append(g,data[grep(a[i],data$substitutions),5])
    h <- append(h,data[grep(a[i],data$substitutions),12])
    j <- append(j,rep(a[i],length(grep(a[i],data$substitutions))))
  }
  
  k <- data.frame(c,d,e,f,g,h,j)
  
  l <- k%>%group_by(h,j)%>%summarise(n=n())
  head(l)
  #plot#
  l <- as.data.frame(l)
  l$date <- date_decimal(l$h)
  names(l) <- c("decimal_date","mutation","freq","week_date")
  head(l)
  
  #freq <- 1#
  mut_freqs <- data.frame(table(k$j))
  mut_freqs$perc <- (mut_freqs$Freq/nrow(data))*100
  mut_selec <- mut_freqs[mut_freqs$perc >= freq_inf | mut_freqs$perc <= freq_sup ,1]
  mut_freqs2 <- mut_freqs[order(-mut_freqs$perc),]
  mut_freqs3 <- mut_freqs2[mut_freqs2$Freq >= freq_inf & mut_freqs2$Freq <= freq_sup ,]
  
  m <- l[l$mutation %in% mut_selec,]
  sites <- unique(mut_freqs3$Var1)
  m2 <- m[m$mutation %in% sites ,]
  fig <- plot_ly(m2,x = as.Date(m2$week_date,format="%Y-%m-%d"),
                 y = ~freq, type = 'scatter', name=~mutation,mode = 'lines') %>% layout(xaxis = list(range = as.POSIXct(c(xmin,xmax))),
                                                                                        yaxis = list(range=c(ymin,ymax))) %>% layout(title = title)
  print(fig)
  print(mut_freqs3)
  write.csv(mut_freqs3,paste0(gene,"_","mutations.csv"),row.names=FALSE)
  write.csv(m2,paste0(gene,"_","mutations_frequencies.csv"),row.names=FALSE)
  print(dim(m))
  print(names(m2))
  print(head(m2))
  print(sites)
}

#example#
#mutations_time_nuc(data=e,xmin="2022-06-15",xmax="2022-09-15",ymin=0,ymax=8,freq_sup=11,freq_inf=2,title="Monkey_Peru",gene="nuc")#

######################################
##### CODE 2 : mutations_time_aa #####
######################################

mutations_time_aa <- function(data,xmin,xmax,ymin,ymax,freq_sup,freq_inf,title,gene){
    data = data 
    xmin = xmin 
    xmax = xmax
    ymin = ymin
    ymax = ymax
    freq_sup = freq_sup
    freq_inf = freq_inf
    title = title
    gene = gene
    region = paste0(gene,"_")
  
  #add new column dates#
  names(data) <- c("ID","epicodes","province","country","lineage","Collection.date","substitutions","aaSubstitutions")
  data <- add_epi_week(data,"Collection.date", system="cdc")
  data$week_date <- epi_week_date(data$epi_week,data$epi_year,system="cdc")
  data <- as.data.frame(data)
  data$dec_date <- decimal_date(ymd(data$week_date))
  data.frame(names(data))
  head(data)
  
  #extract unique mutation lists#
  #separate by commas#
  elements <- unlist(strsplit(data$aaSubstitutions, ","))
  a <- sort(unique(elements))
  b <- unique(gsub(":.*","",a))
  d <- as.data.frame(table(a))
  #write.csv(sort(a),"aa_mundo_mutations.csv",row.names = TRUE)#
  
  #unique(gsub(":.*","",a))#
  
  #generate table#
  data.frame(names(data))
  
  c <- 0
  d <- 0 
  e <- 0 
  f <- 0 
  g <- 0 
  h <- 0 
  j <- 0
  for (i in 1:length(a)){
    c <- append(c,data[grep(a[i],data$aaSubstitutions),1])
    d <- append(d,data[grep(a[i],data$aaSubstitutions),2])
    e <- append(e,data[grep(a[i],data$aaSubstitutions),3])
    f <- append(f,data[grep(a[i],data$aaSubstitutions),4])
    g <- append(g,data[grep(a[i],data$aaSubstitutions),5])
    h <- append(h,data[grep(a[i],data$aaSubstitutions),12])
    j <- append(j,rep(a[i],length(grep(a[i],data$aaSubstitutions))))
  }
  
  k <- data.frame(c,d,e,f,g,h,j)
  
  l <- k%>%group_by(h,j)%>%summarise(n=n())
  #plot#
  l <- as.data.frame(l)
  l$date <- date_decimal(l$h)
  
  names(l) <- c("decimal_date","mutation","freq","week_date")
  
  mut_freqs <- data.frame(table(k$j))
  mut_freqs$perc <- (mut_freqs$Freq/nrow(data))*100
  mut_selec <- mut_freqs[mut_freqs$perc >= freq_inf | mut_freqs$perc <= freq_sup ,1]
  mut_freqs2 <- mut_freqs[order(-mut_freqs$perc),]
  mut_freqs3 <- mut_freqs2[mut_freqs2$Freq >= freq_inf & mut_freqs2$Freq <= freq_sup ,]
  
  m <- l[l$mutation %in% mut_selec,]
  sites <- unique(mut_freqs3$Var1)
  m2 <- m[m$mutation %in% sites ,]
  fig <- plot_ly(m2,x = as.Date(m2$week_date,format="%Y-%m-%d"),
                 y = ~freq, type = 'scatter', name=~mutation,mode = 'lines') %>% layout(xaxis = list(range = as.POSIXct(c(xmin,xmax))),
                                                                                        yaxis = list(range=c(ymin,ymax))) %>% layout(title = title)
                                                                                        
nextstrain <- c("NBT03_gp001","NBT03_gp002","NBT03_gp003","NBT03_gp004","NBT03_gp005","NBT03_gp006","NBT03_gp007","NBT03_gp008","NBT03_gp009",
"NBT03_gp010","NBT03_gp011","NBT03_gp012","NBT03_gp013","NBT03_gp014","NBT03_gp015","NBT03_gp016","NBT03_gp017","NBT03_gp018","NBT03_gp019","NBT03_gp020",
"NBT03_gp021","NBT03_gp022","NBT03_gp023","NBT03_gp024","NBT03_gp025","NBT03_gp026","NBT03_gp027","NBT03_gp028","NBT03_gp029","NBT03_gp030","NBT03_gp031",
"NBT03_gp032","NBT03_gp033","NBT03_gp034","NBT03_gp035","NBT03_gp036","NBT03_gp037","NBT03_gp038","NBT03_gp039","NBT03_gp040","NBT03_gp041","NBT03_gp042",
"NBT03_gp043","NBT03_gp044","NBT03_gp045","NBT03_gp046","NBT03_gp047","NBT03_gp048","NBT03_gp049","NBT03_gp050","NBT03_gp051","NBT03_gp052","NBT03_gp053",
"NBT03_gp054","NBT03_gp055","NBT03_gp056","NBT03_gp057","NBT03_gp058","NBT03_gp059","NBT03_gp060","NBT03_gp061","NBT03_gp062","NBT03_gp063","NBT03_gp064",
"NBT03_gp065","NBT03_gp066","NBT03_gp067","NBT03_gp068","NBT03_gp069","NBT03_gp070","NBT03_gp071","NBT03_gp072","NBT03_gp073","NBT03_gp074","NBT03_gp075",
"NBT03_gp076","NBT03_gp077","NBT03_gp078","NBT03_gp079","NBT03_gp080","NBT03_gp081","NBT03_gp082","NBT03_gp083","NBT03_gp084","NBT03_gp085","NBT03_gp086",
"NBT03_gp087","NBT03_gp088","NBT03_gp089","NBT03_gp090","NBT03_gp091","NBT03_gp092","NBT03_gp093","NBT03_gp094","NBT03_gp095","NBT03_gp096","NBT03_gp097",
"NBT03_gp098","NBT03_gp099","NBT03_gp100","NBT03_gp101","NBT03_gp102","NBT03_gp103","NBT03_gp104","NBT03_gp105","NBT03_gp106","NBT03_gp107","NBT03_gp108",
"NBT03_gp109","NBT03_gp110","NBT03_gp111","NBT03_gp112","NBT03_gp113","NBT03_gp114","NBT03_gp115","NBT03_gp116","NBT03_gp117","NBT03_gp118","NBT03_gp119",
"NBT03_gp120","NBT03_gp121","NBT03_gp122","NBT03_gp123","NBT03_gp124","NBT03_gp125","NBT03_gp126","NBT03_gp127","NBT03_gp128","NBT03_gp129","NBT03_gp130",
"NBT03_gp131","NBT03_gp132","NBT03_gp133","NBT03_gp134","NBT03_gp135","NBT03_gp136","NBT03_gp137","NBT03_gp138","NBT03_gp139","NBT03_gp140","NBT03_gp141",
"NBT03_gp142","NBT03_gp143","NBT03_gp144","NBT03_gp145","NBT03_gp146","NBT03_gp147","NBT03_gp148","NBT03_gp149","NBT03_gp150","NBT03_gp151","NBT03_gp152",
"NBT03_gp153","NBT03_gp154","NBT03_gp155","NBT03_gp156","NBT03_gp157","NBT03_gp158","NBT03_gp159","NBT03_gp160","NBT03_gp161","NBT03_gp162","NBT03_gp163",
"NBT03_gp164","NBT03_gp165","NBT03_gp166","NBT03_gp167","NBT03_gp168","NBT03_gp169","NBT03_gp170","NBT03_gp171","NBT03_gp172","NBT03_gp173","NBT03_gp174",
"NBT03_gp175","NBT03_gp176","NBT03_gp177","NBT03_gp178")

nextclade <- c("OPG001","OPG002","OPG003","OPG015","OPG019","OPG021","OPG022","OPG023","OPG024","OPG025","OPG027","OPG029","OPG030","OPG031","OPG034","OPG035",
"OPG036","OPG037","OPG038","OPG039","OPG040","OPG042","OPG043","OPG044","OPG045","OPG046","OPG047","OPG048","OPG049","OPG050","OPG051","OPG052","OPG053",
"OPG054","OPG055","OPG056","OPG057","OPG058","OPG059","OPG060","OPG061","OPG062","OPG063","OPG064","OPG065","OPG066","OPG068","OPG069","OPG070","OPG071",
"OPG072","gene_biotype=protein_coding","OPG074","OPG075","OPG076","OPG077","OPG078","OPG079","OPG080","OPG081","OPG082","OPG083","OPG084","OPG085","OPG086",
"OPG087","OPG088","OPG089","OPG090","OPG091","OPG092","OPG093","OPG094","OPG095","OPG096","OPG097","OPG098","OPG099","OPG100","OPG101","OPG102","OPG103","OPG104",
"OPG105","OPG106","OPG107","OPG108","OPG109","OPG110","OPG111","OPG112","OPG113","OPG114","OPG115","OPG116","OPG117","OPG118","OPG119","OPG120","OPG121","OPG122",
"OPG123","OPG124","OPG125","OPG126","OPG127","OPG128","OPG129","OPG130","OPG131","OPG132","OPG133","OPG134","OPG135","OPG136","OPG137","OPG138","OPG139","OPG140",
"OPG141","OPG142","OPG143","OPG144","OPG145","OPG146","OPG147","OPG148","OPG149","OPG150","OPG151","OPG153","OPG154","OPG155","OPG156","OPG157","OPG158","OPG159",
"OPG160","OPG161","OPG162","OPG163","OPG164","OPG165","OPG167","OPG170","OPG171","OPG172","OPG173","OPG174","OPG175","OPG176","OPG178","OPG180","OPG181","OPG185",
"OPG187","OPG188","OPG189","OPG190","OPG191","OPG192","OPG193","OPG195","OPG197","OPG198","OPG199","OPG200","OPG204","OPG205","OPG208","OPG209","OPG210","OPG005",
"gene_biotype=protein_coding","gene_biotype=protein_coding","OPG003","OPG002","OPG001")

genes <- data.frame(nextstrain,nextclade)

names(mut_freqs3) <- c("mutation","frequency","percentage")
library(tidyr)
mut_freqs4 <- separate(mut_freqs3,"mutation",c("nextclade","mutation"),sep="[:]")
mut_freqs5 <- merge(mut_freqs4,genes,by="nextclade",all.x=TRUE)
mut_freqs6 <- mut_freqs5[order(mut_freqs5$frequency),]

  print(fig)
  print(mut_freqs3)
  print(mut_freqs6)
  write.csv(mut_freqs3,paste0(gene,"_","mutations.csv"),row.names=FALSE)
  write.csv(mut_freqs6,paste0(gene,"_","mutations_nextstrain.csv"),row.names=FALSE)
  print(dim(m))
  print(names(m2))
  print(sites)
  
}

#example#
#mutations_time_aa(data=e,xmin="2022-06-15",xmax="2022-09-15",ymin=0,ymax=8,freq_sup=11,freq_inf=2,title="Monkey : aasubstitutions",gene="prot")#


######################################
##### CODE 3 : plot_time #####
######################################

plot_time <- function(data=data,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,freq_sup=freq_sup,freq_inf=freq_inf,title=title,gene=gene){
  data = data 
  xmin = xmin 
  xmax = xmax
  ymin = ymin
  ymax = ymax
  freq_sup = freq_sup
  freq_inf = freq_inf
  title = title
  gene = gene
fig <- plot_ly(data,x = as.Date(data$week_date,format="%Y-%m-%d"),
               y = ~freq, type = 'scatter', name=~mutation,mode = 'lines') %>% layout(xaxis = list(range = as.POSIXct(c(xmin,xmax))),
                                                                                      yaxis = list(range=c(ymin,ymax))) %>% layout(title = title)
print(fig)
}