## Script to plot Nanodrop Data from ndv export
## J. Kabisch 2014
## kabisch@uni-greifswald.de
library(reshape)
library(ggplot)
header <- read.table("TALENJK0713033.ndv",nrows = 4, sep="\t",header=F, stringsAsFactors = FALSE)
values <- read.table("TALENJK0713033.ndv", skip=4, sep="\t",header=T, dec = ",", comment.char="?")
#transTest <- read.table("TALENJK0713033_test.csv", skip=4, sep="\t", header=T, dec=",")
curve.values <- subset(values, select = c(Sample.ID, X220:X303)) # extract values for plots
amount <- subset(values, select = c(Sample.ID, ng.ul)) # extract quantity
protein <- subset(values, select = c(X260.280)) # extract protein contamination
md <- melt(curve.values, id=(c("Sample.ID"))) # transform to ggplot format
md$variable<-substring(md$variable,2) # remove X in front of variables
ggplot(md) + geom_line(aes(x=variable, y=value, group=Sample.ID)) + facet_wrap(~Sample.ID) + ylim(-1, 6)

           