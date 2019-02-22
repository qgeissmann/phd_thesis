rm(list=ls())

h_min <- 0.9
h_max <- 1.2
r1 <- 0.01
r2 <- 0.92

t <- seq(from = 0, to = 72, length.out=400)
c <- (1 - sin(t * 2*pi / 24 ))  /2 

# sleeps for 12h, every 12h, at night
p <- rep(F,length(t))
s <- rep(0,length(t))

for(i in 2:length(t)){
	if(!p[i -1]) {
		s[i] <- s[i-1] + (1- s[i-1]) * r1
	}
	else{
		s[i] <- s[i-1] * r2
	}
	if(p[i-1])
		p[i] <- ifelse(s[i] + c[i] > h_min, T, F)
	else
		p[i] <- ifelse(s[i] + c[i] > h_max, T, F)
	}
	
dt <- data.frame(t=t,c=c,s=s, p=p)#,s3=s3)


lim <-  scale_y_continuous(lim=c(0,1))
pl1 <- ggplot(dt, aes(x=t))  + geom_line(aes(y= c)) + lim
pl2 <- ggplot(dt, aes(x=t))  + geom_line(aes(y= s)) + lim
pl3 <- ggplot(dt, aes(x=t))  + geom_line(aes(y= (s + c)/2, col=as.numeric(p)))  + geom_hline(yintercept=c(h_min,h_max)/2, linetype=3)  + lim +  theme(legend.position="none")

pp_baseline <- plot_grid(pl1,pl2,pl3, ncol=1)
pp_baseline



p <- rep(F,length(t))
s <- rep(0,length(t))

for(i in 2:length(t)){
	if(!p[i -1]) {
		s[i] <- s[i-1] + (1- s[i-1]) * r1
	}
	else{
		s[i] <- s[i-1] * r2
	}
	if(p[i-1])
		p[i] <- ifelse(s[i] + c[i] > h_min, T, F)
	else
		p[i] <- ifelse(s[i] + c[i] > h_max & t[i] > 24, T, F)
	}
	
dt <- data.frame(t=t,c=c,s=s, p=p)#,s3=s3)



pl1 <- ggplot(dt, aes(x=t))  + geom_line(aes(y= c)) + lim
pl2 <- ggplot(dt, aes(x=t))  + geom_line(aes(y= s)) + lim
pl3 <- ggplot(dt, aes(x=t/24))  + geom_line(aes(y= (s + c)/2, col=as.numeric(p)))  + geom_hline(yintercept=c(h_min,h_max)/2, linetype=3)  + lim +  theme(legend.position="none")


pp_sd <-  plot_grid(pl1,pl2,pl3, ncol=1)

pdf("two-processes.pdf", w=12, h=8)
plot_grid(pp_baseline, pp_sd, ncol=2)
dev.off()


