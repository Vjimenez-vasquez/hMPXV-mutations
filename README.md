# hMPXV-mutations
```r
A collection of codes for monkeypox virus temporal mutation tracking. 
```

# Example
```r
## 1 : load libraries ##
library(lubridate)
library(plotly)
library(dplyr)
library(epical)

## 1.1 (optional) install packages ##
#install.packages("remotes")#
#remotes::install_github("chrismerkord/epical")#
#install.packages("plotly")#
#install.packages("dplyr")#
#install.packages("lubridate")#

## 2 : load data ##
dir()
source("mutation_commands.R")
m <- read.csv("metadata.tsv", header=TRUE, sep="\t") 
dim(m)
as.data.frame(names(m))

## 3 : Generate a subset metadata with the following information ##
#column 1 : ID 
#column 2 : epicodes 
#column 3 : province 
#column 4 : country
#column 5 : lineage
#column 6 : dates (Collection.date)
#column 7 : $substitutions
#column 8 : $aaSubstitutions

# the subset is called "mut" # 
mut <- m[,c(1,2,6,5,14,3,16,17)]
lin <- sort(unique(m$lineage))
out <- mut[mut$lineage %in% lin[8:22],]
sort(unique(out$lineage))
dim(out)
head(out)
```

# Usage
```r
mutations_time_nuc(data=out,xmin="2022-04-01",xmax="2023-01-30",ymin=1,ymax=200,freq_sup=2800,freq_inf=300,title="Monkey_mundo",gene="nuc_mun")

mutations_time_aa(data=out,xmin="2022-06-15",xmax="2022-09-15",ymin=0,ymax=8,freq_sup=11,freq_inf=2,title="Monkey : aasubstitutions",gene="prot")

# arguments 
data = data frame containing the desired columns
xmin = minimum date to consider in the plot  
xmax = maximum date to consider in the plot
ymin = minimum mutation frequency to consider in the plot 
ymax = maximum mutation frequency to consider in the plot
freq_sup = minimum mutation frequency to consider in the search
freq_inf = maximum mutation frequency to consider in the search
title = title
gene = prefix to consider for the output files
```

# Output
```r


![plot](https://user-images.githubusercontent.com/89874227/221016162-80396ece-d8ca-4335-98d7-05495bb6c8d3.png)
