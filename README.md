# hMPXV-mutations
A code designed by Victor Jimenez Vasquez - vr.jimenez.vs@gmail.com

![graphic1](https://user-images.githubusercontent.com/89874227/221016548-351495d7-5d88-41a4-9920-594d2a1a0948.jpg)

```r

 A collection of codes for monkeypox virus temporal mutation tracking
 
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
# 1 mutations_time_nuc : Estimates the temporal variation of nucleotide substitutions
mutations_time_nuc(data=out,xmin="2022-04-01",xmax="2023-01-30",ymin=1,ymax=200,freq_sup=2800,freq_inf=300,title="Monkey_mundo",gene="nuc_mun")

#2 mutations_time_aa : Estimates the temporal variation of aminoacid mutations
mutations_time_aa(data=out,xmin="2022-06-15",xmax="2022-09-15",ymin=0,ymax=8,freq_sup=11,freq_inf=2,title="Monkey : aasubstitutions",gene="prot")

#3 plot_time : Generates individual plots taking as input "*_mutations_frequencies.csv" file obtained with mutations_time_nuc or mutations_time_aa commands 
plot_time(data=n4,xmin="2022-05-01",xmax="2022-09-10",ymin=1,ymax=3.5,freq_sup=24,freq_inf=10,title="Monkey_B.1.14",gene="nuc_B.1.14")

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
1) *_mutations.csv : data frame with the list of identified mutations and their associated frequencies and percentages. 
2) *_mutations_frequencies.csv : data frame with the identified mutation by epidemiological.
3) a plot showing the temporal variation of the identified mutations. 
```
