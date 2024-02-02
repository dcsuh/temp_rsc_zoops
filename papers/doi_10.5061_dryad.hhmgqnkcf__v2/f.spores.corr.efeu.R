##################################################################
#####################    READ IN THE DATA    #####################
##################################################################

library(gplots)

setwd("C:/Users/straussa/Documents/Research/Hall Lab/f spores/data")
data <- read.csv("fits.summary.csv")
data <- data[data$index!="species" , ] #omit species data
data$round <- as.factor(t(matrix(unlist(strsplit(as.character(data$index), split=" ")), nrow=2))[,1])
data$clone <- as.factor(t(matrix(unlist(strsplit(as.character(data$index), split=" ")), nrow=2))[,2])

repeats <- unique(efeu[duplicated(efeu$clone)==T,]$clone)

length(unique(data$func))
efeu <- data[data$func=="exp.f.exp.u" , ]
repeats <- unique(efeu[duplicated(efeu$clone)==T,]$clone)
singles <- efeu[which(efeu$clone %in% repeats ==F),]$index
efeu$dup <- ifelse(efeu$clone %in% repeats ==T, 1,0)
efeu1 <- efeu[efeu$round==1, ]
efeu2 <- efeu[efeu$round==2, ]
efeu3 <- efeu[efeu$round==3, ]

vol <- 15
z.seq <- seq(from=0,to=393, length.out=50) # sequence to plot for x axis
f.min <- 3.6 # this is the bottom 2.5% 


##################################################################
####################   CORRELATION PLOTS    ######################
##################################################################

par(mfrow=c(2,2),mar=c(0,0,1,1), oma=c(3,3,0,0)) 

plot(efeu$f.0.hat, efeu$alpha, col="white", xlab = "", ylab = "", 
     xlim = c(5,20), ylim = c(-0.012,0.0035), 
     cex=2, pch=21, bg="white", xaxt="n", yaxt="n")
axis(side=2, at=c(0,-0.005,-0.01), cex.axis=1.2)
points(efeu1$f.0.hat, efeu1$alpha, cex=1.6, pch=21, bg="black")
points(efeu2$f.0.hat, efeu2$alpha, cex=1.6, pch=21, bg="black")
points(efeu3$f.0.hat, efeu3$alpha, cex=1.6, pch=21, bg="black")
points(efeu3[efeu3$clone=="Mid.252", ]$f.0.hat, efeu3[efeu3$clone=="Mid.252", ]$alpha, cex=1.6, pch=21, bg="green3") #pull to front
points(efeu3[efeu3$clone=="Mid.244", ]$f.0.hat, efeu3[efeu3$clone=="Mid.244", ]$alpha, cex=1.6, pch=21, bg="dodgerblue") #pull to front
points(efeu2[efeu2$clone=="Bris.10", ]$f.0.hat, efeu2[efeu2$clone=="Bris.10", ]$alpha, cex=1.6, pch=21, bg="dark orange") #pull to front

plot(efeu$w, efeu$alpha, col="white", xlab = "", ylab = "", 
     xlim = c(-0.001,0.004),ylim = c(-0.012,0.0035), 
     cex=0, pch=21, bg="white", xaxt="n", yaxt="n",
     panel.first = rect(0, -10,10,10, col='snow2', border=NA))
points(efeu1$w, efeu1$alpha, cex=1.6, pch=21, bg="black")
points(efeu2$w, efeu2$alpha, cex=1.6, pch=21, bg="black")
points(efeu3$w, efeu3$alpha, cex=1.6, pch=21, bg="black")
points(efeu3[efeu3$clone=="Mid.252", ]$w, efeu3[efeu3$clone=="Mid.252", ]$alpha, cex=1.6, pch=21, bg="green3") #pull to front
points(efeu3[efeu3$clone=="Mid.244", ]$w, efeu3[efeu3$clone=="Mid.244", ]$alpha, cex=1.6, pch=21, bg="dodgerblue") #pull to front
points(efeu2[efeu2$clone=="Bris.10", ]$w, efeu2[efeu2$clone=="Bris.10", ]$alpha, cex=1.6, pch=21, bg="dark orange") #pull to front

plot(efeu$f.0.hat, efeu$u, col="white", xlab = "", ylab = "", 
     xlim = c(5,20), ylim = c(-0.0005,0.01), 
     cex=2, pch=21, bg="white", xaxt="n", yaxt="n")
axis(side=1, at=c(6,12,18), cex.axis=1.2)
axis(side=2, at=c(0,0.005, 0.01), cex.axis=1.2)
points(efeu1$f.0.hat, efeu1$u, cex=1.6, pch=21, bg="black")
points(efeu2$f.0.hat, efeu2$u, cex=1.6, pch=21, bg="black")
points(efeu3$f.0.hat, efeu3$u, cex=1.6, pch=21, bg="black")
points(efeu3[efeu3$clone=="Mid.252", ]$f.0.hat, efeu3[efeu3$clone=="Mid.252", ]$u, cex=1.6, pch=21, bg="green3") #pull to front
points(efeu3[efeu3$clone=="Mid.244", ]$f.0.hat, efeu3[efeu3$clone=="Mid.244", ]$u, cex=1.6, pch=21, bg="dodgerblue") #pull to front
points(efeu2[efeu2$clone=="Bris.10", ]$f.0.hat, efeu2[efeu2$clone=="Bris.10", ]$u, cex=1.6, pch=21, bg="dark orange") #pull to front

plot(efeu$w, efeu$u, col="white", xlab = "", ylab = "", 
     xlim = c(-0.001,0.004), ylim = c(-0.0005,0.01), 
     cex=0, pch=21, bg="white", xaxt="n", yaxt="n",
     panel.first = rect(0,-10,10,10, col='snow2', border=NA))
axis(side=1, at=c(0.0, 0.002,0.004), cex.axis=1.2)
points(efeu1$w, efeu1$u, cex=1.6, pch=21, bg="black")
points(efeu2$w, efeu2$u, cex=1.6, pch=21, bg="black")
points(efeu3$w, efeu3$u, cex=1.6, pch=21, bg="black")
points(efeu3[efeu3$clone=="Mid.252", ]$w, efeu3[efeu3$clone=="Mid.252", ]$u, cex=1.6, pch=21, bg="green3") #pull to front
points(efeu3[efeu3$clone=="Mid.244", ]$w, efeu3[efeu3$clone=="Mid.244", ]$u, cex=1.6, pch=21, bg="dodgerblue") #pull to front
points(efeu2[efeu2$clone=="Bris.10", ]$w, efeu2[efeu2$clone=="Bris.10", ]$u, cex=1.6, pch=21, bg="dark orange") #pull to front


#############################################################################################
# relationships among other parameters?

par(mfrow=c(1,2),mar=c(0,8,1,1), oma=c(3,0,0,0)) 

plot(efeu$alpha, efeu$u, col="white", xlab = "", ylab = "", 
     xlim = c(-0.012,0.0035), ylim = c(-0.0005,0.01), 
     cex=2, pch=21, bg="white", xaxt="n", yaxt="n")
axis(side=1, at=c(0, -0.005, -0.01), cex.axis=1)
axis(side=2, at=c(0,0.005, 0.01), cex.axis=1)
points(efeu1$alpha, efeu1$u, cex=1.5, pch=21, bg="black")
points(efeu2$alpha, efeu2$u, cex=1.5, pch=21, bg="black")
points(efeu3$alpha, efeu3$u, cex=1.5, pch=21, bg="black")
points(efeu3[efeu3$clone=="Mid.252", ]$alpha, efeu3[efeu3$clone=="Mid.252", ]$u, cex=1.6, pch=21, bg="green3") #pull to front
points(efeu3[efeu3$clone=="Mid.244", ]$alpha, efeu3[efeu3$clone=="Mid.244", ]$u, cex=1.6, pch=21, bg="dodgerblue") #pull to front
points(efeu2[efeu2$clone=="Bris.10", ]$alpha, efeu2[efeu2$clone=="Bris.10", ]$u, cex=1.6, pch=21, bg="dark orange") #pull to front

plot(efeu$w, efeu$f.0.hat, col="white", xlab = "", ylab = "", 
     xlim = c(-0.001,0.004), ylim = c(5,20), 
     cex=0, pch=21, bg="white", xaxt="n", yaxt="n",
     panel.first = rect(0, -10, 30, 30, col='snow2', border=NA))
axis(side=2, at=c(6,12,18), cex.axis=1)
axis(side=1, at=c(0,0.002,0.004), cex.axis=1)
points(efeu1$w, efeu1$f.0.hat, cex=1.5, pch=21, bg="black")
points(efeu2$w, efeu2$f.0.hat, cex=1.5, pch=21, bg="black")
points(efeu3$w, efeu3$f.0.hat, cex=1.5, pch=21, bg="black")
points(efeu3[efeu3$clone=="Mid.252", ]$w, efeu3[efeu3$clone=="Mid.252", ]$f.0.hat, cex=1.6, pch=21, bg="green3") #pull to front
points(efeu3[efeu3$clone=="Mid.244", ]$w, efeu3[efeu3$clone=="Mid.244", ]$f.0.hat, cex=1.6, pch=21, bg="dodgerblue") #pull to front
points(efeu2[efeu2$clone=="Bris.10", ]$w, efeu2[efeu2$clone=="Bris.10", ]$f.0.hat, cex=1.6, pch=21, bg="dark orange") #pull to front


####################################################################
########## how consistent are clones among rounds? #################
####################################################################

efeu.dup <- efeu[efeu$dup==1,]
unique(efeu.dup$clone)
efeu.dup$col <- "black"
efeu.dup$col <- ifelse(efeu.dup$clone=="A45", "purple", efeu.dup$col)
efeu.dup$col <- ifelse(efeu.dup$clone=="Bris.10", "dark orange", efeu.dup$col)
efeu.dup$col <- ifelse(efeu.dup$clone=="Bris.6", "green3", efeu.dup$col)
efeu.dup$col <- ifelse(efeu.dup$clone=="Bris.112", "dodgerblue", efeu.dup$col)
efeu.dup$col <- ifelse(efeu.dup$clone=="War.5", "white", efeu.dup$col)
efeu.dup$col <- ifelse(efeu.dup$clone=="Dog.4", "black", efeu.dup$col)

efeu.dup1 <- efeu.dup[efeu.dup$round==1, ]
efeu.dup2 <- efeu.dup[efeu.dup$round==2, ]
efeu.dup3 <- efeu.dup[efeu.dup$round==3, ]


par(mfrow=c(1,2),mar=c(0,8,1,1), oma=c(3,0,0,0)) 

plot(efeu.dup$alpha, efeu.dup$u, col="white", xlab = "", ylab = "", 
     xlim = c(-0.008,0.001), ylim = c(-0.0005,0.005), 
     cex=2, pch=21, bg="white", xaxt="n", yaxt="n")
axis(side=1, at=c(0, -0.002, -0.004), cex.axis=1)
axis(side=2, at=c(0,0.002, 0.004), cex.axis=1)
points(efeu.dup1$alpha, efeu.dup1$u, cex=1.5, pch=22, bg=efeu.dup1$col)
points(efeu.dup2$alpha, efeu.dup2$u, cex=1.5, pch=24, bg=efeu.dup2$col)
points(efeu.dup3$alpha, efeu.dup3$u, cex=1.5, pch=25, bg=efeu.dup3$col)

plot(efeu.dup$w, efeu.dup$f.0.hat, col="white", xlab = "", ylab = "", 
     xlim = c(-0.001,0.003), ylim = c(5,17), 
     cex=0, pch=21, bg="white", xaxt="n", yaxt="n",
     panel.first = rect(0, -10, 30, 30, col='snow2', border=NA))
axis(side=2, at=c(6,11,16), cex.axis=1)
axis(side=1, at=c(0,0.001,0.002), cex.axis=1)
points(efeu.dup1$w, efeu.dup1$f.0.hat, cex=1.5, pch=22, bg=efeu.dup1$col)
points(efeu.dup2$w, efeu.dup2$f.0.hat, cex=1.5, pch=24, bg=efeu.dup2$col)
points(efeu.dup3$w, efeu.dup3$f.0.hat, cex=1.5, pch=25, bg=efeu.dup3$col)

##################################################################
########################  REACTION NORMS  ########################
##################################################################

par(mfrow=c(4,1),mar=c(4,3,1,1), oma=c(0,0,0,0)) 

#############
###   F   ###
#############

Fs=matrix(ncol=length(z.seq), nrow=nrow(efeu))
rownames(Fs) <- efeu$index
for (i in 1:nrow(efeu)){
  Fs[i,]=pmax(efeu$f.0.hat[i]*exp(efeu$alpha[i]*z.seq*(efeu$length[i]^2))*(efeu$length[i]^2), f.min)} 

plot(z.seq, Fs[1,], type="l",xlab = "", ylab = "", xaxt="n", yaxt="n", ylim=c(0,25)) # A45 R1
axis(side=2, at=c(0,10,20), cex.axis=1.2)
axis(side=1, at=c(0,200,400), cex.axis=1.2)
lines(z.seq, Fs[2,]) # Bris.10 R1
lines(z.seq, Fs[3,]) # Bris.112 R1
lines(z.seq, Fs[4,]) # Bris.6 R1
lines(z.seq, Fs[5,]) # Cback.256 R1
lines(z.seq, Fs[6,]) # Dog.4 R1
lines(z.seq, Fs[7,]) # Mid.276 R1
lines(z.seq, Fs[8,]) # Mid.277 R1
lines(z.seq, Fs[9,]) # War.5 R1
lines(z.seq, Fs[10,]) # A45 R2
lines(z.seq, Fs[11,]) # A48 R2
lines(z.seq, Fs[13,]) # Dog.4 R2
lines(z.seq, Fs[14,]) # Is.278 R2
lines(z.seq, Fs[15,]) # Mid.263 R2
lines(z.seq, Fs[16,]) # Mid.273 R2
lines(z.seq, Fs[17,]) # A43 R3
lines(z.seq, Fs[18,]) # A45 R3
lines(z.seq, Fs[19,]) # Bris.10 R3
lines(z.seq, Fs[20,]) # Bris.111 R3
lines(z.seq, Fs[21,]) # Bris.112 R3
lines(z.seq, Fs[22,]) # Bris.6 R3
lines(z.seq, Fs[23,]) # Cback.276 R3
lines(z.seq, Fs[26,]) # Mid.281 R3
lines(z.seq, Fs[27,]) # War.5 R3
lines(z.seq, Fs[12,], lwd=3, col="dark orange") # Bris.10 R2
lines(z.seq, Fs[24,], lwd=3, col="dodgerblue") # Mid.244 R3
lines(z.seq, Fs[25,], lwd=3, col="green3") # Mid.252 R3

#############
###   U   ###
#############

# first, spore consumption rate for each genotype (the x axes):
Zs=matrix(ncol=length(z.seq), nrow=nrow(efeu))
rownames(Zs) <- efeu$index
for (i in 1:nrow(efeu)){
  Zs[i,]=pmax(efeu$f.0.hat[i]*exp(efeu$alpha[i]*z.seq*(efeu$length[i]^2))*(efeu$length[i]^2), f.min)*z.seq} 

# then, U per spore consumed:
Us=matrix(ncol=length(z.seq), nrow=nrow(efeu))
rownames(Us) <- efeu$index
for (i in 1:nrow(efeu)){
  Us[i,]=efeu$u[i]*exp(efeu$w[i]*z.seq*pmax(efeu$f.0.hat[i]*exp(efeu$alpha[i]*z.seq*(efeu$length[i]^2))*(efeu$length[i]^2), f.min))} 

plot(Zs[1,], Us[1,], type="l",xlab = "", ylab = "", xaxt="n", yaxt="n", 
     ylim=c(0,3.5e-3), xlim=c(0,8000))
axis(side=2, at=c(0,2e-3), cex.axis=1.2)
axis(side=1, at=c(0,3000,6000), cex.axis=1.2)
lines(Zs[2,], Us[2,]) # Bris.10 R1
lines(Zs[3,], Us[3,]) # Bris.112 R1
lines(Zs[4,], Us[4,]) # Bris.6 R1
lines(Zs[5,], Us[5,]) # Cback.256 R1
lines(Zs[6,], Us[6,]) # Dog.4 R1
lines(Zs[7,], Us[7,]) # Mid.276 R1
lines(Zs[8,], Us[8,]) # Mid.277 R1
lines(Zs[9,], Us[9,]) # War.5 R1
lines(Zs[10,], Us[10,]) # A45 R2
lines(Zs[11,], Us[11,]) # A48 R2
lines(Zs[13,], Us[13,]) # Dog.4 R2
lines(Zs[14,], Us[14,]) # Is.278 R2
lines(Zs[15,], Us[15,]) # Mid.263 R2
lines(Zs[16,], Us[16,]) # Mid.273 R2
lines(Zs[17,], Us[17,]) # A43 R3
lines(Zs[18,], Us[18,]) # A45 R3
lines(Zs[19,], Us[19,]) # Bris.10 R3
lines(Zs[20,], Us[20,]) # Bris.111 R3
lines(Zs[21,], Us[21,]) # Bris.112 R3
lines(Zs[22,], Us[22,]) # Bris.6 R3
lines(Zs[23,], Us[23,]) # Cback.276 R3
lines(Zs[26,], Us[26,]) # Mid.281 R3
lines(Zs[27,], Us[27,]) # War.5 R3
lines(Zs[12,], Us[12,], lwd=3, col="dark orange") # Bris.10 R2
lines(Zs[24,], Us[24,], lwd=3, col="dodgerblue") # Mid.244 R3
lines(Zs[25,], Us[25,], lwd=3, col="green3") # Mid.252 R3

#############
### BETA  ###
#############

betas=matrix(ncol=length(z.seq), nrow=nrow(efeu))
rownames(betas) <- efeu$index
for (i in 1:nrow(efeu)){
  betas[i,]=pmax(efeu$f.0.hat[i]*exp(efeu$alpha[i]*z.seq*(efeu$length[i]^2))*(efeu$length[i]^2), f.min)*
    efeu$u[i]*exp(efeu$w[i]*z.seq*pmax(efeu$f.0.hat[i]*exp(efeu$alpha[i]*z.seq*(efeu$length[i]^2))*(efeu$length[i]^2), f.min))} 

plot(z.seq, betas[1,], type="l",xlab = "", ylab = "", xaxt="n", yaxt="n", ylim=c(0,0.06)) # A45 R1
axis(side=2, at=c(0,0.03), cex.axis=1.2)
axis(side=1, at=c(0,200,400), cex.axis=1.2)
lines(z.seq, betas[2,]) # Bris.10 R1
lines(z.seq, betas[3,]) # Bris.112 R1
lines(z.seq, betas[4,]) # Bris.6 R1
lines(z.seq, betas[5,]) # Cback.256 R1
lines(z.seq, betas[6,]) # Dog.4 R1
lines(z.seq, betas[7,]) # Mid.276 R1
lines(z.seq, betas[8,]) # Mid.277 R1
lines(z.seq, betas[9,]) # War.5 R1
lines(z.seq, betas[10,]) # A45 R2
lines(z.seq, betas[11,]) # A48 R2
lines(z.seq, betas[13,]) # Dog.4 R2
lines(z.seq, betas[14,]) # Is.278 R2
lines(z.seq, betas[15,]) # Mid.263 R2
lines(z.seq, betas[16,]) # Mid.273 R2
lines(z.seq, betas[17,]) # A43 R3
lines(z.seq, betas[18,]) # A45 R3
lines(z.seq, betas[19,]) # Bris.10 R3
lines(z.seq, betas[20,]) # Bris.111 R3
lines(z.seq, betas[21,]) # Bris.112 R3
lines(z.seq, betas[22,]) # Bris.6 R3
lines(z.seq, betas[23,]) # Cback.276 R3
lines(z.seq, betas[26,]) # Mid.281 R3
lines(z.seq, betas[27,]) # War.5 R3
lines(z.seq, betas[12,], lwd=3, col="dark orange") # Bris.10 R2
lines(z.seq, betas[24,], lwd=3, col="dodgerblue") # Mid.244 R3
lines(z.seq, betas[25,], lwd=3, col="green3") # Mid.252 R3

#############
# Infection #
#############

trs=matrix(ncol=length(z.seq), nrow=nrow(efeu))
rownames(trs) <- efeu$index
for (i in 1:nrow(efeu)){
  trs[i,]=pmax(efeu$f.0.hat[i]*exp(efeu$alpha[i]*z.seq*(efeu$length[i]^2))*(efeu$length[i]^2), f.min)*
    efeu$u[i]*exp(efeu$w[i]*z.seq*pmax(efeu$f.0.hat[i]*exp(efeu$alpha[i]*z.seq*(efeu$length[i]^2))*(efeu$length[i]^2), f.min))*
    z.seq} 

plot(z.seq, trs[1,], type="l",xlab = "", ylab = "", xaxt="n", yaxt="n", ylim=c(0,4)) # A45 R1
axis(side=1, at=c(0,200,400), cex.axis=1.2)
axis(side=2, at=c(0,2), cex.axis=1.2)
lines(z.seq, trs[2,]) # Bris.10 R1
lines(z.seq, trs[3,]) # Bris.112 R1
lines(z.seq, trs[4,]) # Bris.6 R1
lines(z.seq, trs[5,]) # Cback.256 R1
lines(z.seq, trs[6,]) # Dog.4 R1
lines(z.seq, trs[7,]) # Mid.276 R1
lines(z.seq, trs[8,]) # Mid.277 R1
lines(z.seq, trs[9,]) # War.5 R1
lines(z.seq, trs[10,]) # A45 R2
lines(z.seq, trs[11,]) # A48 R2
lines(z.seq, trs[13,]) # Dog.4 R2
lines(z.seq, trs[14,]) # Is.278 R2
lines(z.seq, trs[15,]) # Mid.263 R2
lines(z.seq, trs[16,]) # Mid.273 R2
lines(z.seq, trs[17,]) # A43 R3
lines(z.seq, trs[18,]) # A45 R3
lines(z.seq, trs[19,]) # Bris.10 R3
lines(z.seq, trs[20,]) # Bris.111 R3
lines(z.seq, trs[21,]) # Bris.112 R3
lines(z.seq, trs[22,]) # Bris.6 R3
lines(z.seq, trs[23,]) # Cback.276 R3
lines(z.seq, trs[26,]) # Mid.281 R3
lines(z.seq, trs[27,]) # War.5 R3
lines(z.seq, trs[12,], lwd=3, col="dark orange") # Bris.10 R2
lines(z.seq, trs[24,], lwd=3, col="dodgerblue") # Mid.244 R3
lines(z.seq, trs[25,], lwd=3, col="green3") # Mid.252 R3


