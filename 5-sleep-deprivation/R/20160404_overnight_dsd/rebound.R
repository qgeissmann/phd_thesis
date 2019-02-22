rm(list=ls());gc()

library(behavr)
METADATA <- "metadata.csv"
CACHE <- "./cache/"
RESULT_DIR <- "/data/ethoscope_results"
REMOTE_DIR <- "ftp://nas.lab.gilest.ro/auto_generated_data/ethoscope_results/"
SD_START_END <- days(1:2/2)
REF_SLEEP_WINDOW <- hours(c(0, 3))
TARGET_SLEEP_WINDOW <- REF_SLEEP_WINDOW + days(1)
PDF_NAME <- "overnight_dsd.pdf"

source("../trunk.R")

dt <- dt[xmv(interval) %in% c("20", "Control")]
make_plot()





dt_ord_2 <- dt[,
			{
			dd <- copy(.SD[, .(t, x_rel, moving, interactions)])
			dd[, t := (floor(t/60))*60]
			dd[,
				.(
					x_rel=mean(x_rel),
					walked_dist = sum(abs(diff(x_rel))),
					behaviour = ifelse(all(!moving),1,2),
					interactions = sum(interactions)
					),
					by="t"]
			},
			by="id"]


dt_ord_2[, behaviour := ifelse(behaviour != 1, (walked_dist > .25) + 2  , 1)]
dt_ord_2[, behaviour:= ordered(c("q","m","w")[behaviour], levels = c("q","m","w"))]
dt_ord_2[, micro.mov. := (behaviour == "m")]
dt_ord_2[, walking := (behaviour == "w")]
dt_ord_2[, quiescent := (behaviour == "q")]

stat_rebound_dt <- rejoin(
						dt_ord_2[,
						.(
						
						quiet_baseline_night = mean(quiescent[t %between% (SD_START_END - days(1))]),
						quiet_baseline_day3h = mean(quiescent[t %between% REF_SLEEP_WINDOW]),
						quiet_sd_night = mean(quiescent[t %between% SD_START_END]),
						quiet_rebound_day3h = mean(quiescent[t %between% TARGET_SLEEP_WINDOW]),
						interactions = sum(interactions)
						)
						,by = id]
						)


mod <- lm(quiet_rebound_day3h ~ quiet_baseline_day3h * sex, stat_rebound_dt[sdi == 0])
stat_rebound_dt[, quiet_rebound_day3h_pred := predict(mod, stat_rebound_dt)]
stat_rebound_dt[, quiet_rebound_day3h_diff := quiet_rebound_day3h - quiet_rebound_day3h_pred]



mod <- lm(quiet_sd_night ~ quiet_baseline_night * sex, stat_rebound_dt[sdi == 0])
stat_rebound_dt[, quiet_sd_night_pred := predict(mod, stat_rebound_dt)]
stat_rebound_dt[, quiet_sd_night_diff := quiet_sd_night - quiet_sd_night_pred]

stat_rebound_dt <- stat_rebound_dt[interval != "Control"]

stat_rebound_dt[, interval2 := as.numeric(as.character(interval))]

layers <- list( geom_jitter(alpha=.3),
				facet_grid(. ~ sex),
				stat_summary(fun.y = "mean", geom="point", shape = 4, size=2.5, colour="black"),
				stat_summary(fun.data = "mean_cl_boot", geom="errorbar",colour="black")
				)

		
pl1 <- ggplot(stat_rebound_dt, aes(interval, interactions, colour=interval2)) + 
		layers +
		scale_y_sqrt(name=expression(N[stimuli]), breaks=c(10, 50, 100, 500, 1000)) 
		


pl2 <- ggplot(stat_rebound_dt, aes(interval, quiet_rebound_day3h_diff * 180, colour= interval2)) + 
		layers +
		scale_y_continuous(name=expression(Delta~immobility[rebound]~(min)), limits=c(-180, 180)) +
		geom_hline(yintercept=0, colour="red", linetype=2) 


sel_dt <- stat_rebound_dt[ interval %in% c("20","60", "220", "420", "680", "1000")]
 + scale_colour_gradient(limits=c(0, 1000), high="red", low="blue")

	pl3 <- ggplot(sel_dt, aes( interactions, quiet_sd_night, colour= interval2)) + 
				scale_y_continuous(name=expression(P(quiescent)[SD])) +
				scale_x_sqrt(name=expression(N[stimuli]), breaks=c(10, 50, 100, 500, 1000)) +
				geom_point(alpha=.3) + geom_smooth(aes(group=interval),method = "lm") + facet_grid(. ~ sex)
			
pl4 <- ggplot(sel_dt, 
		aes( quiet_sd_night_diff * 12 * 60, quiet_rebound_day3h_diff * 3 * 60, colour= interval2)) + 
		geom_point(alpha=.3) + 
		geom_smooth(aes(group=interval), method="lm")+ 
		facet_grid(. ~ sex) +
		geom_hline(yintercept=0, colour="red", linetype=2) 
		
library(cowplot)


pdf("overnight_sd_intervals", w=12, h=8)
	print(
		plot_grid(pl1,pl2,pl3,pl4, ncol=2, labels=LETTERS[1:4])
		)
dev.off()


ggplot(stat_rebound_dt[interval %in% c("20", "120", "300", "540", "840")], aes( quiet_baseline_night, as.numeric(quiet_rebound_day3h_diff> 0), colour= interval)) + 
		geom_point() + geom_smooth(method = "glm", method.args = list(family= "binomial") )+ facet_grid(. ~ sex) 

ggplot(stat_rebound_dt[interval %in% c("20", "120", "300", "540", "840")], aes( quiet_sd_night, as.numeric(quiet_rebound_day3h_diff> 0), colour= interval)) + 
		geom_point() + geom_smooth(method = "glm", method.args = list(family= "binomial") )+ facet_grid(. ~ sex) 


 

ggplot(stat_rebound_dt[interval %in% c("20", "120", "300", "540", "840")], aes( quiet_baseline_night, quiet_rebound_day3h_diff, colour= interval)) + 
		geom_point() + 
		geom_smooth()+ 
		facet_grid(interval ~ sex) 


ggplot(stat_rebound_dt[interval %in% c("20", "120", "300", "540", "840")], aes( interactions, quiet_rebound_day3h_diff, colour= interval)) + 
		geom_point() + 
		scale_x_sqrt(name=expression(N[stimuli]), breaks=c(10, 50, 100, 500, 1000)) +
		geom_smooth(method="lm")+ 
		facet_grid(interval ~ sex) 
		
ggplot(stat_rebound_dt[interval %in% c("20", "120", "300", "540", "840")], aes( interactions, quiet_baseline_night, colour= interval)) + 
		scale_y_continuous(name=expression(P(quiescent)[baseline])) +
		scale_x_sqrt(name=expression(N[stimuli]), breaks=c(10, 50, 100, 500, 1000)) +
		geom_point() + geom_smooth(method = "lm") + facet_grid(. ~ sex)

#~ ggplot(stat_rebound_dt[interval %in% c("20", "120", "300", "540", "840")], aes( quiet_baseline_night, quiet_sd_night, colour= interval)) + 
#~ #		scale_y_continuous(name=expression(P(quiescent)[baseline])) +
#~ #		scale_x_sqrt(name=expression(N[stimuli]), breaks=c(10, 50, 100, 500, 1000)) +
#~ 		geom_point() + geom_smooth(method = "lm") + facet_grid(. ~ sex)


#~ ggplot(stat_rebound_dt[interval %in% c("20", "120", "300", "540", "840")], aes( interactions, quiet_rebound_day3h_diff * 180, colour= interval)) + 
#~ 		scale_y_continuous(name=expression(Delta~immobility[rebound]~(min))) +
#~ 		scale_x_sqrt(name=expression(N[stimuli]), breaks=c(10, 50, 100, 500, 1000)) +
#~ 		geom_hline(yintercept=0, colour="red", linetype=2) + facet_grid(. ~ sex) +
#~ 		geom_point() + geom_smooth(method = "lm")


#~ lm(quiet_rebound_day3h_diff ~quiet_sd_night_diff, stat_rebound_dt[sex=="M"])$coefficients
#~ foo <- function(d, i){
#~ 	lm(quiet_rebound_day3h_diff ~quiet_sd_night_diff, d, subset=i)$coefficients
#~ 	}
	
#~ boot::boot(data=stat_rebound_dt[sex=="M"], statistic=foo, R=1000)


