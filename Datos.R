setwd("F:/ETAPA 2/R/DEG/RSEM")

T0 <- as.data.frame(read.csv("res05UD-T0.csv"))
T1 <- as.data.frame(read.csv("res05UD-T1.csv"))
T24 <- as.data.frame(read.csv("res05UD-T24.csv"))
T72 <- as.data.frame(read.csv("res05UD-T72.csv"))

print("Positivo")
sum(T0$FDR<0.05 & T0$logFC>0)
print("Negativo")
sum(T0$FDR<0.05 & T0$logFC<0)
print("Total")
sum(T0$FDR<0.05)

print("Positivo")
sum(T1$FDR<0.05 & T1$logFC>0)
print("Negativo")
sum(T1$FDR<0.05 & T1$logFC<0)
print("Total")
sum(T1$FDR<0.05)

print("Positivo")
sum(T24$FDR<0.05 & T24$logFC>0)
print("Negativo")
sum(T24$FDR<0.05 & T24$logFC<0)
print("Total")
sum(T24$FDR<0.05)

print("Positivo")
sum(T72$FDR<0.05 & T72$logFC>0)
print("Negativo")
sum(T72$FDR<0.05 & T72$logFC<0)
print("Total")
sum(T72$FDR<0.05)


