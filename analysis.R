###CP216 Preferentially Essential Gene Analysis - DepMap Gene Effect - 9-15-2021
​
library(data.table)
library(dplyr)
​
setwd("C:/Users/penneyc/Downloads")
​
#Achilles_Gene_Effect
AGE <- fread("achilles_gene_effect.csv")
#Achilles_Common_Essentials
ACE <- read.csv("achilles_common_essentials.csv")
#Subsetting Achilles_Gene_Effect table to remove common-essential gene columns.
are_common_essentials <- which(colnames(AGE) %in% ACE$gene)
AGE_mod <- AGE %>% select(-all_of(are_common_essentials))
​
#Main code:
​
#Initialize empty DF
df<-AGE_mod[FALSE,2:length(colnames(AGE_mod))]
​
#Set Counter
i <- 1
​
#While i is less that the number of cell lines (rows)
#Take the colMeans of all gene-effects except the ith row
#Then get the difference between those means and the ith rows GE values
#Add that difference row to the empty df we created
#Add 1 to counter
​
while(i <= length(rownames(AGE_mod))){
  geneMeans <- t(colMeans(AGE_mod[-i,2:length(colnames(AGE_mod))], na.rm = TRUE))
  diffs <- geneMeans-AGE_mod[i,2:length(colnames(AGE_mod))]
  df<-rbind(df,diffs)
  print(i)
  i <- i + 1}
​
#View(df)
​
#Write file
#write.csv(final, "analysis_table_out.csv")
​
####EXTRA Checking step.
​
      #a <- AGE_mod[-i,2:length(colnames(AGE_mod))]
      #View(a)
      #b <- colMeans(a, na.rm = TRUE)
      #View(b)
      # #c <- t(b)
      # View(c)
      # d<- AGE_mod[i,2:length(colnames(AGE_mod))]
      # View(d)
      # e <- c-d
      # View(e)
​
####____________
​
#Now, clear out environment and reload that table we just created and stored. Transpose and wrangle so genes are rows and columns are cell lines.
​
file <- fread("analysis_table_out.csv", header = TRUE)
file <- t(file)
colnames(file) <- file[1,]
file <- file[-1,]
​
#Initialize empty list
empty <- c()
​
#Initialize counter
i<-1
​
#Go column by column and sort the absolute values for each column by descending, then grab the top 10 gene names (row names) from that column. Store each to the empty list we had so its a list of lists. 
while(i <= 625){
  sub <- as.data.frame(file[,i])
  View(sub)
  sub$gene <- rownames(sub)
  sub$abs <- abs(as.numeric(as.character(sub$`file[, i]`)))
  sub <- sub[order(-sub$abs),]
  top10genes <- sub$gene[1:10]
  #View(top10genes)
  empty[[i]] <- top10genes
  #View(empty)
  #View(sub)
  print(i)
  i <- i + 1
}
​
#rbind list of lists to df.
df <- do.call(rbind, empty)
​
#Add DepMap_IDs back in.
df <- cbind(colnames(file),df)
​
#Add in DepMap cell line info:
sample_info <- read.csv("sample_info.csv")
sample_info <- sample_info[,c(1,3)]
merge <- merge(sample_info,df,by.x = "DepMap_ID",by.y = "V1" )
​
#Correct one empty gene name cell (see https://depmap.org/portal/cell_line/ACH-001173 for info) 
merge$stripped_cell_line_name[534] <- "U251MGDM"
​
#Add informative column names and write
colnames(merge) <- c("DepMap_ID", "Stripped Cell Line Name", 
                  "Gene_1",
                  "Gene_2",
                  "Gene_3",
                  "Gene_4",
                  "Gene_5",
                  "Gene_6",
                  "Gene_7",
                  "Gene_8",
                  "Gene_9",
                  "Gene_10")
#write.csv(merge, "CP216_Preferentially_Essential_Genes_Analysis_Table.csv", col.names = FALSE, row.names = FALSE)
​
#Final checking:
check <- fread("CP216_Preferentially_Essential_Genes_Analysis_Table.csv")
View(check)
