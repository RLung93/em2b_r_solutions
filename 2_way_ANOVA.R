library("scatterplot3d") # load
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
library(RColorBrewer)
n <- 15
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

FeedStock <- rep(c("I","II","III","IV","V"), each = 2, times = 3)
Process <- rep(c("A","B","C"), each = 10)
Process.FeedStock <- paste(Process, FeedStock, sep ="-")
Yield <- c(106,110,95,100,94,107,103,104,100,102,110,112,98,99,100,
           101,108,112,105,107,94,97,86,87,98,99,99,101,94,98)
full.chem.data <- data.frame(ID = FeedStock.Process,
                             FeedStock = factor(FeedStock),
                             Process = factor(Process),
                             Yield = Yield)
plot(data=full.chem.data, Yield~FeedStock)
chemmod <- lm(Yield~FeedStock+Process+FeedStock:Process)

# single SumOfSquares
# decomposition
CrossMeans <- rep(tapply(Yield, Process.FeedStock, mean), each = 2)
ProcessMeans <- rep(tapply(Yield, Process, mean), each = 10)
FeedStockMeans <- rep(tapply(Yield, FeedStock, mean), each = 2, times = 3)
null.mean <- mean(Yield)

SS.FeedStock <- sum((FeedStockMeans - null.mean)^2)
df.FeedStock <- length(unique(FeedStock))-1
MS.FeedStock <- SS.FeedStock/df.FeedStock

SS.Process <- sum((ProcessMeans - null.mean)^2)
df.Process <- length(unique(Process))-1
MS.Process <- SS.Process/df.Process

SS.Cross <- sum((CrossMeans - ProcessMeans - FeedStockMeans + null.mean)^2)
df.Cross <- length(unique(Process.FeedStock)) - (df.FeedStock + df.Process) - 1
MS.Cross <- SS.Cross/df.Cross

SS.error <- sum((Yield - CrossMeans)^2)
df.error <- length(Yield) - length(unique(Process.FeedStock))
MS.error <- SS.error/df.error

SS.total <- sum((Yield - mean(Yield))^2)
df.total <- length(Yield) - 1
  
F.FeedStock <- MS.FeedStock/MS.error
P.FeedStock <- 1-pf(F.FeedStock, df1 = df.FeedStock, df2 = df.error)
F.Process <- MS.Process/MS.error
P.Process <- 1-pf(F.Process, df1 = df.Process, df2 = df.error)
F.Cross <- MS.Cross/MS.error
P.Cross <- 1-pf(F.Cross, df1 = df.Cross, df2 = df.error)


data.frame(SumOfSquares = c(SS.FeedStock, SS.Process, SS.Cross, SS.error, SS.total),
           DegreesOfFreedom = c(df.FeedStock, df.Process, df.Cross, df.error, df.total),
           MeanSquare = c(MS.FeedStock, MS.Process, MS.Cross, MS.error, NA),
           FStatistic = c(F.FeedStock, F.Process, F.Cross, NA, NA),
           PValue = c(P.FeedStock, P.Process, P.Cross, NA, NA))

anova(chemmod)
