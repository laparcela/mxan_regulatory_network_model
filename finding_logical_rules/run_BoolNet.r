library(BoolNet) 

#perl truth_tables.pl 3 | grep -v '*' | grep -v 'x' > Rules-3-Regulators.txt 
#perl truth_tables.pl 4 | grep -v '*' | grep -v 'x' > Rules-4-Regulators.txt 
#perl truth_tables.pl 2 | grep -v '*' | grep -v 'x' > Rules-2-Regulators.txt 
#perl truth_tables.pl 1 | grep -v '*' | grep -v 'x' > Rules-1-Regulators.txt 

one.r   = read.csv("Rules-1-Regulators.txt", header = F, sep =" ")
two.r   = read.csv("Rules-2-Regulators.txt", header = F, sep =" ")
three.r = read.csv("Rules-3-Regulators.txt", header = F, sep =" ")
four.r  = read.csv("Rules-4-Regulators.txt", header = F, sep =" ")

pep.three = subset(three.r, three.r[,9] == '-')

sink("testNet.bn")
cat("targets, factors\n")
cat("MKAPA, PKTA2 & PKTD9\n") ###
cat("MKAPB, PKTA4 & PKTD1\n") ###
cat("MKAPC, PKTC2\n") 
cat("PKTA2, MRPC2\n") 
cat("PKTA4, MKAPA\n") 
cat("PKTC2, MKAPA & MKAPB\n") ###
cat("PKTD1, PEP & MKAPA & MKAPC\n") ###
cat("PKTD9, !STARV\n")
cat("PSKA5, PKTC2\n")
cat("STARV, STARV\n")
cat("PEP, !STARV\n")
cat("MRPC2, PSKA5\n")
sink()
# read file
net = loadNetwork("testNet.bn")
print(net)

ptm = proc.time()
for(mkapA in 1:nrow(two.r)){
    net$interactions$MKAPA$func = two.r[mkapA, 1:4]
    for(mkapB in 1:nrow(two.r)){
        net$interactions$MKAPB$func = two.r[mkapB, 1:4]
        for(pktc2 in 1:nrow(two.r)){
            net$interactions$PKTC2$func = two.r[pktc2, 1:4]
            for(pktd1 in 1:nrow(pep.three)){
                net$interactions$PKTD1$func = pep.three[pktd1, 1:8]
                attr = getAttractors(net)
                print(attr)
                print(c(mkapA, mkapB, pktc2, pktd1))
            }
        }
    }
}

proc.time() - ptm
