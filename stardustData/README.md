# StardustData
In this repository are collected all the datasets discussed in "Stardust: improving spatial tranScripTomics data analysis through space awARe modularity optimization baseD clUSTering."  
The datasets are: Lymph Node, Human Brest Cancer 1, Human Brest Cancer 2, Mouse Kidney and Human Heart.  
For each dataset is included:
* Image of the tissue used in the 10X Visium protocol,
* 2D coordinates of each spot,
* a compressed version of the expression matrix (uncompression is required in order to let Stardust work properly)

Moreover, under directory "Results" is collected all the data used in the method validation (Wilcoxon rank sum statistical test over 100 positional shufllings per dataset). Results subdirectories have been cleaned from unnecessary files.

In order to reproduce violins plot figures you can use the following bash and R code:
```bash
# cd into the dataset directory you are interested on
cd Results/MouseKidney
R
```
```R
library(ggplot2)
clusters <- list.dirs(path = paste(getwd(),"/0/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
zero <- read.delim(paste(getwd(),"/0/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)

clusters <- list.dirs(path = paste(getwd(),"/0.25/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
zerotwenyfive <- read.delim(paste(getwd(),"/0.25/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)

clusters <- list.dirs(path = paste(getwd(),"/0.5/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
zerofifty <- read.delim(paste(getwd(),"/0.5/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)

clusters <- list.dirs(path = paste(getwd(),"/0.75/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
zeroseventyfive <- read.delim(paste(getwd(),"/0.75/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)

clusters <- list.dirs(path = paste(getwd(),"/1/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
one <- read.delim(paste(getwd(),"/1/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)

spots <- dim(zero)[1]
x = factor(c(rep("0:1",spots),rep("1:4",spots),rep("1:2",spots),rep("3:4",spots),rep("1:1",spots)),
          levels=c("0:1","1:4","1:2","3:4","1:1"))

df = data.frame(x=x,y=c(zero$V2,zerotwenyfive$V2,zerofifty$V2,zeroseventyfive$V2,one$V2),
                physicalDistance = c(rep("no",spots),rep("yes",spots*4)))

ggplot(df,aes(x=x, y=y, fill=physicalDistance)) + geom_violin(trim=TRUE) + 
  geom_boxplot(width=0.1, fill="white") + stat_summary(fun=mean, geom="line", aes(group=1)) + 
  stat_summary(fun=mean, geom="point")  +  labs(x="Stardust configuration - space:transcripts", y = "Cell Stability Score")

```

In order to reproduce statistical evaluation figures you can use the following bash and R code:
```bash
# cd into the dataset directory you are interested on
cd Results/MouseKidney
R
```

```R
wilcox0 = ""
wilcox0.25 = ""
wilcox0.5 = ""
wilcox0.75 = ""
wilcox1 = ""

clusters <- list.dirs(path = paste(getwd(),"/0/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
compare0 <- read.delim(paste(getwd(),"/0/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)

clusters <- list.dirs(path = paste(getwd(),"/0.25/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
compare0.25 <- read.delim(paste(getwd(),"/0.25/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)

clusters <- list.dirs(path = paste(getwd(),"/0.5/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
compare0.5 <- read.delim(paste(getwd(),"/0.5/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)

clusters <- list.dirs(path = paste(getwd(),"/0.75/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
compare0.75 <- read.delim(paste(getwd(),"/0.75/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)

clusters <- list.dirs(path = paste(getwd(),"/1/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
clusters <- strsplit(clusters,"/")
clusters <- tail(clusters[[1]], n=1)
compare1 <- read.delim(paste(getwd(),"/1/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)


for(j in c(0,0.25,0.5,0.75,1)){
  w_temp = 1:100
  compare = get(paste("compare",j,sep=""))
  for(i in 101:200){
    clusters <- list.dirs(path = paste(getwd(),"/",j,"null/",i,"/Results/filtered_expression_matrix/",sep=""), full.names = TRUE, recursive = TRUE)[2]
    clusters <- strsplit(clusters,"/")
    clusters <- tail(clusters[[1]], n=1)
    a <- read.delim(paste(getwd(),"/",j,"null/",i,"/Results/filtered_expression_matrix/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep=""), header=FALSE, row.names=1)
    k = wilcox.test(compare$V2, y = a$V2 ,alternative = "greater")
    w_temp[i-100] = k$p.value
  }
  assign(paste("wilcox",j,sep=""), w_temp)
}
frame = data.frame(configuration=factor(c(rep("0:1",100),rep("1:4",100),rep("1:2",100),rep("3:4",100),rep("1:1",100)),levels=c("0:1","1:4","1:2","3:4","1:1")))
frame$test = rep(c(rep("wilcoxon",100)),5)
pvalue_correction = 100
frame$significant = c(wilcox0 < (0.05/pvalue_correction),wilcox0.25 < (0.05/pvalue_correction),wilcox0.5 < (0.05/pvalue_correction),wilcox0.75 < (0.05/pvalue_correction),wilcox1 < (0.05/pvalue_correction))
count = table(frame)
count = t(count[,,2])
rownames(count) = c("Wilcoxon")
barplot(count, main="Number of significant tests over 100 per configuration",
        xlab="Stardust configuration - space:transcripts", col=c("#2fbfc4"), ylim=c(0,100),
        legend = rownames(count), beside=TRUE,args.legend = list( x = "top",bty = "n", ncol = 1,inset = -0.1))
```

If you want to reproduce the entire statistical evaluation for a single dataset you can use the following bash and R code (it could take more than 24h depending on your machine specifications)
```bash
mkdir stat_evaluation && cd stat_evaluation
# download the dataset you are interested
wget https://github.com/GiovanniCS/StardustData/raw/main/Datasets/MouseKidney/filtered_expression_matrix.txt.zip
unzip filtered_expression_matrix.txt.zip
rm -rf __MACOSX && rm filtered_expression_matrix.txt.zip
wget https://raw.githubusercontent.com/GiovanniCS/StardustData/main/Datasets/MouseKidney/spot_coordinates.txt

R
```

```R
install.packages("devtools")
devtools::install_github("GiovanniCS/rCASC")
library(rCASC)

root_dir <- getwd()
for(j in c(0,0.25,0.5,0.75,1)){
    for(i in 101:200){        
        system(paste("mkdir -p ",root_dir,"/",j,"null/",i,"/scratch",sep=""))
        setwd(paste(root_dir,"/",j,"null/",i,sep=""))
        system(paste("cp ",root_dir,"/filtered_expression_matrix.txt .",sep=""))
        system(paste("cp ",root_dir,"/spot_coordinates.txt .",sep=""))
        scratch.folder <- paste(getwd(),"/scratch",sep="")
        file <- paste(getwd(),"/filtered_expression_matrix.txt",sep="")
        tissuePosition <- paste(getwd(),"/spot_coordinates.txt",sep="")
        shuf <- read.delim(tissuePosition, row.names=1)
        rownames(shuf) = sample(rownames(shuf))
        write.table(shuf, tissuePosition, append = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE)
        StardustPermutation(group="docker",scratch.folder = scratch.folder,nPerm=80,
                            file=file, tissuePosition=tissuePosition, spaceWeight=j, 
                            permAtTime=20, percent=10, separator="\t",
                            logTen=0, pcaDimensions=5, seed=111)
        cluster.path <- paste(data.folder=dirname(file), "Results", strsplit(basename(file),"\\.")[[1]][1], sep="/")
        cluster <- as.numeric(list.dirs(cluster.path, full.names = FALSE, recursive = FALSE))
        permAnalysisSeurat(group="docker",scratch.folder = scratch.folder,file=file, nCluster=cluster,separator="\t",sp=0.8)
        system("rm filtered_expression_matrix.txt")
        system("rm Results/filtered_expression_matrix.txt")
        system("rm Results/spot_coordinates.txt")
    }
}
```