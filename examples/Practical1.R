#############################################################
### Target trial emulation and causal inference practical ###
#############################################################

#################
### Section 1 ###
#################

## load data
lalonde <- read.table(file = "~/Downloads/lalonde_nsw.csv") # change filepath accordingly

## extract variables
y <- lalonde$re78
t <- lalonde$treat