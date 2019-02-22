rm(list=ls());gc()
options(nwarnings = 1000)
library(scopr)
library(ggetho)
library(sleepr)
library(gtools)
library(cowplot)
#source("../ggplot_themes.R")

no_legend = list(theme(legend.position="none"))

layer_barpl <- function(){
	list(
		facet_grid( . ~ sex),
		geom_jitter(alpha=.3, height=0),
		stat_summary(fun.y = "mean", geom="point", shape = 4, size=2.5, colour="black"),
		stat_summary(fun.data = "mean_cl_boot", geom="errorbar",colour="black"),
		scale_fill_manual(values=CONTROL_SD_PALETTE),
		scale_colour_manual(values=CONTROL_SD_PALETTE)
		)
	}
	

layers <- function(palette = CONTROL_SD_PALETTE, annotate=TRUE, sd_start_end = SD_START_END){
  out <- list(
    stat_pop_etho(method= mean_cl_boot),
    facet_grid( sex ~ .),
    stat_ld_annotations(),
    #coord_cartesian(xlim = c(days(-1),days(2))),
    scale_y_continuous(limits = c(NA,1)),
    scale_fill_manual(values=palette),
    scale_colour_manual(values=palette)
    )
  if(annotate)
    out <- c(out, list(
        annotate("segment",y = .9, yend = .9,   x = sd_start_end[1], xend = sd_start_end[2],   colour = "black",alpha=0.5,size=3),
        annotate("text",y=0.95,x=mean(sd_start_end), label="treatment")
    ))
    out
}

FEMALE_MALE_PALETTE <- c("#be2828ff", "#282896ff")
CONTROL_SD_PALETTE <- c( "#969696ff", "#3caa3cff")
METADATA <- "metadata.csv"
CACHE <- "./cache/"
#~ RESULT_DIR <- "./raw_results/"
RESULT_DIR <- "/data/ethoscope_results/"
REMOTE_DIR <- "ftp://nas.lab.gilest.ro/auto_generated_data/ethoscope_results/"

SD_START_END <- days(c(.5,10))
REF_SLEEP_WINDOW <- hours(c(0, 3))
TARGET_SLEEP_WINDOW <- REF_SLEEP_WINDOW + days(10)

met <- fread(METADATA)
met <- met[status == "OK"]

#~ met <- link_ethoscope_metadata_remote(met,
#~                                       remote_dir =  REMOTE_DIR,
#~                                       result_dir = RESULT_DIR,
#~                                       verbose = TRUE)

met <- link_ethoscope_metadata(met,
                                      result_dir = RESULT_DIR)                                                                            
dt <- load_ethoscope(met,
					   max_time=days(20),
					   reference_hour=9.0,
					   cache = CACHE,
					   FUN = sleep_annotation,
					   ncores=1)

summary(dt)
                                      
                             
dt <- stitch_on(dt, on="uid")   

curate_data <- function(data){
  
	data <- data[is_interpolated == F]
	data[, treatment := as.factor(ifelse(sdi == 0, "Control", "SD")), meta=T]
	data[, interval := round(((11-sdi) ^ 1.7)) * 20, meta=T]	
    data[, interval := plyr::revalue(as.factor(interval),c("1180"="Control")), meta=T]
	data[, death_date := ifelse(!is.na(as.Date(death_date)), as.Date(death_date),Inf), meta=T]
	data[, lifespan := as.numeric(death_date) - as.numeric(as.Date(datetime)), meta=T]	
	data <- data[t < days(xmv(lifespan)- xmv(baseline_days))]
	data[, t := t - days(xmv(baseline_days))]
	data <- data[t < days(14)]
	data[,x_rel:=ifelse(xmv(region_id) > 10, 1-x, x)]

	norm_x <- function(x){
		min_x <- quantile(na.omit(x),probs=0.01)
		x <- x - min_x
		max_x <- quantile(na.omit(x),probs=0.99)
		x / max_x
		}
	data[, x_rel := norm_x(x_rel), by=id]
	data

}

# then we apply this function to our data
dt <- curate_data(dt)





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



				
overview_pl <- ggetho(dt_ord_2, aes(t, quiescent, colour=treatment)) + 
					facet_grid(sex ~ .) +
					scale_fill_manual(values=CONTROL_SD_PALETTE) +
					scale_colour_manual(values=CONTROL_SD_PALETTE) +
					stat_pop_etho(method = mean_cl_boot) + 
					scale_y_continuous(limits=c(0,1)) +
					scale_x_days() + theme(legend.position = "none")

###############


tdt <- dt_ord_2[t >  days(1) & t <days(10)]


etho_stimuli <- ggetho(tdt,
					aes(y = interactions * 30, fill=treatment),
					summary_FUN = mean, time_wrap=hours(24)) +
					layers(annotate=FALSE)+
					scale_y_continuous(limits = c(NA,NA), name=expression(N[stimuli]~(30~min^{-1})))  +
					no_legend + facet_grid( . ~ sex)
etho_stimuli				
l <- list(quiescent ="P(quiescent)",
          micro.mov. ="P(micro-moving)",
	      walking ="P(walking)",
		  x_rel ="Position (rel.)")

all_data<- lapply(1:length(l), function(i){
	print(i)
	out <- ggetho(tdt, aes_string(y=names(l)[i], colour="treatment"), time_wrap=hours(24))$data
	setnames(out, names(l)[i], "value")
	behav <- tdt[tdt[,names(l)[i], with=F][[1]], behaviour][1]
	if(names(l)[i] == "x_rel")
		behav <- "x_rel"
		
	out[, behaviour := behav]
	})



all_data <- rbindlist(all_data)
all_data_in_sd <- copy(all_data)

main_plot <- ggplot(all_data, aes(t , value, fill =treatment, colour=treatment)) + 
				facet_grid(behaviour ~ sex) +
				scale_fill_manual(values=CONTROL_SD_PALETTE) +
				scale_colour_manual(values=CONTROL_SD_PALETTE) +
				stat_pop_etho(method = mean_cl_boot) + 
				scale_y_continuous(limits=c(0,1)) +
				scale_x_hours() + theme(legend.position = "none")
#############

tdt <- dt_ord_2[t >  days(10) & t <  days(13)]

all_data<- lapply(1:length(l), function(i){
	print(i)
	out <- ggetho(tdt, aes_string(y=names(l)[i], colour="treatment"))$data
	setnames(out, names(l)[i], "value")
	behav <- tdt[tdt[,names(l)[i], with=F][[1]], behaviour][1]
	if(names(l)[i] == "x_rel")
		behav <- "x_rel"
		
	out[, behaviour := behav]
	})



 all_data <- rbindlist(all_data)


main_plot_reb <- ggplot(all_data, aes(t , value, fill =treatment, colour=treatment)) + 
				facet_grid(behaviour ~ sex) +
				scale_fill_manual(values=CONTROL_SD_PALETTE) +
				scale_colour_manual(values=CONTROL_SD_PALETTE) +
				stat_pop_etho(method = mean_cl_boot) + 
				scale_y_continuous(limits=c(0,1)) +
				scale_x_days() + theme(legend.position = "none")



stat_rebound_dt <- rejoin(
						dt_ord_2[,
						.(
						quiet_baseline_day3h = mean(quiescent[t %between% REF_SLEEP_WINDOW]),
						quiet_rebound_day3h = mean(quiescent[t %between% TARGET_SLEEP_WINDOW])
						#interactions = sum(interactions)
						)
						,by = id]
						)
stat_rebound_dt[, treatment := as.character(ifelse(interval =="Control", "Control", "SD"))]

									
summ_stat <- stat_rebound_dt[,{
	w = wilcox.test(quiet_rebound_day3h_diff, mu=0, alternative="greater")
	list(pval=w$p.value, n=.N)},
	keyby="sex,interval,treatment"]
print(summ_stat)			

summ_stat[, text:=sprintf("%s\nN=%03d", gtools::stars.pval(pval), n)]

bar_quiet_reb_day3h_diff <- ggplot(stat_rebound_dt, aes(x=treatment, quiet_rebound_day3h_diff * 3 * 60, colour=treatment)) + layer_barpl() +
			geom_hline(yintercept=0, colour="red", linetype=2) + scale_y_continuous(name=expression(Delta~immobility[rebound]~(min)), limits=c(-180, 180)) + 
			geom_label(data= summ_stat, aes(x=treatment, label=text), y=-Inf, inherit.aes=F, vjust = -0.1) +
			scale_x_discrete(name="Treatment")

pl_row2 <- plot_grid( main_plot, main_plot_reb,
					  etho_stimuli , bar_quiet_reb_day3h_diff + no_legend,
					  labels=LETTERS[c(2,4,3,5)],
						  nrow=2, rel_heights=c(3,1))
		
pl <- plot_grid(overview_pl, pl_row2, nrow=2, rel_heights=c(1,4), labels=c("A", ""))		  
pdf("long-sd.pdf", w=12, h=10)
pl
dev.off()
						  




library(survminer)
library(survival)
library(randomForest)

#surv_data <- fread("/tmp/prolonged_sd_stat_dt.csv")
surv_data <- dt[meta=T]
surv_data[, dead := 1+ !is.infinite(lifespan)] # dead = 2, censored = 1
surv_data[, lifespan_baseline := lifespan - baseline_days]
surv_data[, lifespan_baseline := ifelse(is.infinite(lifespan), 10,lifespan - baseline_days)]
s <- survfit(Surv(lifespan_baseline, dead) ~ sex + treatment, data = surv_data)
ggsurv <- ggsurvplot(s, data=surv_data, conf.int = TRUE, palette=rep(CONTROL_SD_PALETTE, 2))
pl <- ggsurv$plot  + facet_grid( sex ~ .) + theme_grey()
pl_surv <- pl + annotate("segment",y = .70, yend = .70,   x = .5, xend = 10,   colour = "black",alpha=0.5,size=3) +
        annotate("text",y=0.75,x= 9.5/2, label="treatment") +
        scale_y_continuous(labels = scales::percent)	
        

summ_dt_global <- all_data_in_sd[, .(global_avg = mean(value)), by="id,behaviour"]
corr_data <- summ_dt_global[surv_data[, list(id, lifespan_baseline)], on="id"]
corr_data <- corr_data[dt[, .(interactions = sum(interactions)), by=id], on = "id"]
corr_data <- na.omit(corr_data[dt[meta=T], on="id"][is.finite(lifespan)])


behav_lifespan <- ggplot(corr_data, aes(global_avg, lifespan, colour=treatment)) + geom_point() +
						facet_grid(behaviour ~ sex) + geom_smooth(method="lm") +
						scale_fill_manual(values=CONTROL_SD_PALETTE) +
						scale_colour_manual(values=CONTROL_SD_PALETTE)
			
sink("lifespan_stats.txt")
mod <- lm(lifespan ~ (global_avg : behaviour)  * sex, corr_data[treatment == "Control"])
summary(mod)

mod <- lm(lifespan ~ (global_avg : behaviour)  * sex, corr_data[treatment == "SD"])
summary(mod)

corr_data[, sex := as.factor(sex)]
rf <- randomForest(lifespan ~ global_avg * behaviour  * sex * treatment, corr_data, ntree=5000)
rf
col = CONTROL_SD_PALETTE[2]
corr_data[, k_interactions := interactions/1000]
interaction_lifespan_pl <- ggplot(corr_data[treatment == "SD"], aes(k_interactions, lifespan)) + geom_point(colour=col) +
								facet_grid(. ~ sex) + geom_smooth(method="lm",colour=col, fill=col) 
			
mod <- lm(lifespan ~  sex * k_interactions, corr_data[treatment == "SD"])
summary(mod)
sink()

pl <- plot_grid(
				plot_grid(
					pl_surv + no_legend,
					interaction_lifespan_pl, labels=LETTERS[1:2], nrow=2),
				behav_lifespan + no_legend, labels=c("", "C"), ncol=2)		  

pdf("long-sd-lifespan.pdf", w=12, h=8)
pl
dev.off()
				

#~ ###################################
#~ all_pl_objs <- list()


#~ # a set of layers or our next big plots
#~ layers <- function(palette = CONTROL_SD_PALETTE, annotate=TRUE){
#~   out <- list(
#~     stat_pop_etho(method= mean_cl_boot),
#~     facet_grid( sex ~ .),
#~     stat_ld_annotations(),
#~     coord_cartesian(xlim = c(days(-3),days(12))),
#~     scale_y_continuous(limits = c(NA,1)),
#~     scale_fill_manual(values=palette),
#~     scale_colour_manual(values=palette), 
#~     ethogram_theme
#~     )
#~   if(annotate)
#~     out <- c(out, list(
#~         annotate("segment",y = .9, yend = .9,   x = days(.5), xend = days(10),   colour = "black",alpha=0.5,size=3),
#~         annotate("text",y=0.95,x=days(9.5/2), label="treatment")
#~     ))
#~     out
#~ }





#~ zoomed_cartesian <-  coord_cartesian(xlim = c(days(9.5),days(11)))

#~ all_pl_objs$etho_sleep <- ggetho(dt[xmv(sdi) %in% c(0,10)],
#~                                     aes(y = asleep, fill=treatment)) +
#~ 									layers()

#~ all_pl_objs$etho_sleep_zoom <- all_pl_objs$etho_sleep +
#~                                   zoomed_cartesian

#~ all_pl_objs$etho_quiet <- ggetho(dt[xmv(sdi) %in% c(0,10)],
#~                                     aes(y = !moving, fill=treatment)) +
#~ 									layers()
#~ all_pl_objs$etho_quiet_zoom <- all_pl_objs$etho_quiet +
#~                                     zoomed_cartesian


#~ all_pl_objs$etho_stimuli <- ggetho(dt[xmv(sdi) %in% c(0,10)],
#~                                     aes(y = interactions, fill=treatment),
#~                                     summary_FUN = sum) +
#~                                   layers() + scale_y_continuous(limits = c(NA,NA))

#~ all_pl_objs$etho_stimuli_zoom <- all_pl_objs$etho_stimuli +
#~                                         zoomed_cartesian
                                        

#~ dt[, interval := as.numeric(as.character(xmv(interval)))]



#~ all_pl_objs$etho_stimuli_rel_ovrw <- ggetho(dt[t > days(1) &  t < days(10) & xmv(treatment) == "SD"],
#~                                     aes(z = (interactions / 10) / (1/as.numeric(interval))),
#~                                     summary_FUN = mean) +
#~                                     stat_tile_etho()
                                    

#~ all_pl_objs$etho_stimuli_rel <- ggetho(dt[t > days(1) &  t < days(10) & xmv(treatment) == "SD"],
#~                                     aes(y = (interactions / 10) / (1/as.numeric(interval))),
#~                                     summary_FUN = mean, time_wrap= hours(24)) +
#~                                   layers(annotate=F) + scale_y_continuous(name="N_interactions / N_maximum" ,limits = c(NA,NA)) + coord_cartesian(xlim = c(days(0),days(1)))


#~ all_pl_objs$etho_stimuli_rel_day_7_to_10 <- ggetho(dt[t > days(7) &  t < days(10) & xmv(treatment) == "SD"],
#~                                     aes(y = (interactions / 10) / (1/as.numeric(interval))),
#~                                     summary_FUN = mean, time_wrap= hours(24)) +
#~                                   layers(annotate=F) + scale_y_continuous(name="N_interactions / N_maximum" ,limits = c(NA,NA)) + coord_cartesian(xlim = c(days(0),days(1)))


#~ dt[, distance := abs(c(0, diff(x))), by=id]

#~ all_pl_objs$etho_distance <- ggetho(dt[xmv(sdi) %in% c(0,10)],
#~                                     aes(y = distance, fill=treatment),
#~                                     summary_FUN = sum) +
#~                                   layers() + scale_y_continuous(limits = c(NA,NA))

#~ all_pl_objs$etho_distance_wrapped <- ggetho(dt[t > days(1) &  t < days(10) & xmv(sdi) %in% c(0,10)],
#~                                     aes(y = distance / 9, fill=treatment),
#~                                     summary_FUN = sum,  time_wrap= hours(24)) +  
#~ 									layers() + scale_y_continuous(limits = c(NA,NA)) +
#~ 									coord_cartesian(xlim = c(days(0),days(1)))

#~ all_pl_objs$etho_distance_wrapped_day_7_to_10 <- ggetho(dt[t > days(7) &  t < days(10) & xmv(sdi) %in% c(0,10)],
#~                                     aes(y = distance / 3, fill=treatment),
#~                                     summary_FUN = sum,  time_wrap= hours(24)) +  
#~ 									layers() + scale_y_continuous(limits = c(NA,NA)) +
#~ 									coord_cartesian(xlim = c(days(0),days(1)))







all_pl_objs$etho_stimuli_zoom <- all_pl_objs$etho_stimuli +
                                        zoomed_cartesian
                                        
                                        
                                        

#~ ################### scalar stats here and barplots (1 value/animal)

#~ # shared layers
#~ layer_barpl <- function(){
#~ 	list(
#~ 		facet_grid(sex ~ .),
#~ 		geom_jitter(alpha=.3, height=0),
#~ 		stat_summary(fun.y = "mean", geom="point", shape = 4, size=2.5, colour="black"),
#~ 		stat_summary(fun.data = "mean_cl_boot", geom="errorbar",colour="black"),
#~ 		scale_fill_manual(values=CONTROL_SD_PALETTE),
#~ 		scale_colour_manual(values=CONTROL_SD_PALETTE), generic_theme
#~ 		)
#~ 	}
	
	
#~ stat_rebound_dt <- na.omit(rejoin(
#~ 						dt[,
#~ 						.(
#~ 						overall_sleep = mean(asleep),
#~ 						sleep_baseline_day = mean(asleep[t %between%  hours(c(0, 12))]),
#~ 						sleep_baseline_night = mean(asleep[t %between% hours(c(-12, 0))]),
#~ 						sleep_baseline_day3h = mean(asleep[t %between% hours(c(0, 3))]),
#~ 						quiet_baseline_day3h = mean(!moving[t %between% hours(c(0, 3))]),
#~ 						sleep_baseline_day6h = mean(asleep[t %between% hours(c(0, 6))]),
#~ 						sleep_baseline_all = mean(asleep[t %between% hours(c(-12, 12))]),
#~ 						sleep_rebound_day3h = mean(asleep[t %between% (days(9) + hours(c(24, 24 +3)))]),
#~ 						quiet_rebound_day3h = mean(!moving[t %between% (days(9) + hours(c(24, 24 +3)))]),
#~ 						sleep_rebound_day6h = mean(asleep[t %between% (days(9) + hours(c(24, 24 +6)))]),
#~ 						quiet_rebound_day6h = mean(!moving[t %between% (days(9) + hours(c(24, 24 +6)))]),
#~ 						interactions = sum(interactions)
#~ 						)
#~ 						,by = id]
#~ 						))


	
#~ all_pl_objs$bar_sleep_reb_day3h <- ggplot(stat_rebound_dt, aes(interval, sleep_rebound_day3h, colour=treatment)) + layer_barpl()
#~ all_pl_objs$bar_sleep_reb_day6h <- ggplot(stat_rebound_dt, aes(interval, sleep_rebound_day6h, colour=treatment)) + layer_barpl()
#~ all_pl_objs$bar_quiet_reb_day3h <- ggplot(stat_rebound_dt, aes(interval, quiet_rebound_day3h, colour=treatment)) + layer_barpl()
#~ all_pl_objs$bar_quiet_reb_day6h <- ggplot(stat_rebound_dt, aes(interval, quiet_rebound_day6h, colour=treatment)) + layer_barpl()
#~ all_pl_objs$bar_interactions <- ggplot(stat_rebound_dt, aes(interval, interactions, colour=treatment)) + layer_barpl() + scale_y_sqrt()

#~ # quiet is more linear!!
#~ all_pl_objs$bar_quiet_reb_day3h <- ggplot(stat_rebound_dt, aes(interval, quiet_rebound_day3h, colour=treatment)) + layer_barpl()
#~ all_pl_objs$bar_quiet_reb_day3h <- ggplot(stat_rebound_dt, aes(interval, quiet_rebound_day3h, colour=treatment)) + layer_barpl()
#~ all_pl_objs$bar_quiet_reb_day6h <- ggplot(stat_rebound_dt, aes(interval, quiet_rebound_day6h, colour=treatment)) + layer_barpl()






#~ ggplot(stat_rebound_dt, aes(x=quiet_baseline_day3h, y=quiet_rebound_day3h)) + geom_point(alpha=.3) + facet_grid(sex ~ .) + scale_x_sqrt()+ scale_y_sqrt()

#~ ggplot(stat_rebound_dt[sdi %in% c(0, 10)], aes(x=quiet_baseline_day3h, y=quiet_rebound_day3h, colour=treatment)) + geom_point(alpha=1) + facet_grid(sex ~ .) + 
#~ 			#geom_smooth(data=stat_rebound_dt[treatment=="Control"])
#~ 			geom_smooth(method = "lm")




mod <- lm(quiet_rebound_day3h ~ quiet_baseline_day3h * sex, stat_rebound_dt[sdi == 0])

stat_rebound_dt[, quiet_rebound_day3h_pred := predict(mod, stat_rebound_dt)]
stat_rebound_dt[, quiet_rebound_day3h_diff := quiet_rebound_day3h - quiet_rebound_day3h_pred]
stat_rebound_dt[, quiet_rebound_day3h_pred_sign := quiet_rebound_day3h_diff > 0]



#~ all_pl_objs$bar_quiet_reb_day3h_min_baseline <- ggplot(stat_rebound_dt, aes(interval, quiet_rebound_day3h - quiet_baseline_day3h, colour=treatment)) + layer_barpl()
#~ all_pl_objs$bar_quiet_reb_day3h_diff <- ggplot(stat_rebound_dt, aes(interval, quiet_rebound_day3h_diff * 3 * 60, colour=treatment)) + layer_barpl() +
#~ 						geom_hline(yintercept=0, colour="red", linetype=2) + scale_y_continuous(name="Extra quiescence in 3h (min)")


#~ summ_stat <- stat_rebound_dt[,{
#~ 				w = wilcox.test(quiet_rebound_day3h_diff, mu=0, alternative="greater")
#~ 				list(pval=w$p.value, n=.N)},
#~ 				keyby="sex,interval,treatment"]
#~ summ_stat[, text:=sprintf("%s\nN=%03d", stars.pval(pval), n)]
#~ all_pl_objs$bar_quiet_reb_day3h_diff <- all_pl_objs$bar_quiet_reb_day3h_diff + geom_label(data= summ_stat, aes(x=interval, label=text),  y=-0.25 * 180, colour="black")


#~ all_pl_objs$bar_quiet_has_reb <- ggplot(stat_rebound_dt, aes(interval, as.numeric(quiet_rebound_day3h_pred_sign > 0), colour=treatment)) + layer_barpl() + 
#~ 						scale_y_continuous(labels = scales::percent, name="Animals with positive rebound")	 +
#~ 						geom_hline(yintercept=0.5, colour="red", linetype=2)

#~ summ_stat <- stat_rebound_dt[,{
#~ 				w= binom.test(sum(quiet_rebound_day3h_pred_sign),.N, alternative="greater")		
#~ 				list(
#~ 				pval=w$p.value, 
#~ 				n=.N, s=sum(quiet_rebound_day3h_pred_sign))},
#~ 				keyby="sex,interval,treatment"]
				
#~ summ_stat[, text:=sprintf("%s\nN=%03d", stars.pval(pval), n)]

#~ all_pl_objs$bar_quiet_has_reb <- all_pl_objs$bar_quiet_has_reb  + geom_label(data= summ_stat, aes(x=interval, label=text),  y=0.25, colour="black")


#~ pdf("rebound_quantif.pdf", w=8,h=8)
#~ all_pl_objs$bar_quiet_reb_day3h 
#~ all_pl_objs$bar_quiet_reb_day3h_min_baseline
#~ all_pl_objs$bar_quiet_reb_day3h_diff 
#~ all_pl_objs$bar_quiet_has_reb
#~ dev.off()



library(survminer)
library(survival)
library(randomForest)

#surv_data <- fread("/tmp/prolonged_sd_stat_dt.csv")
surv_data <- dt[meta=T]
surv_data[, dead := 1+ !is.infinite(lifespan)] # dead = 2, censored = 1
surv_data[, lifespan_baseline := lifespan - baseline_days]
surv_data[,   := ifelse(is.infinite(lifespan), 10,lifespan - baseline_days)]
s <- survfit(Surv(lifespan_baseline, dead) ~ sex + treatment, data = surv_data)
ggsurv <- ggsurvplot(s, data=surv_data, conf.int = TRUE, palette=rep(CONTROL_SD_PALETTE, 2))
pl <- ggsurv$plot  + facet_grid( sex ~ .) + theme_grey()
pl_surv <- pl + annotate("segment",y = .70, yend = .70,   x = .5, xend = 10,   colour = "black",alpha=0.5,size=3) +
        annotate("text",y=0.75,x= 9.5/2, label="treatment") +
        scale_y_continuous(labels = scales::percent)	

#~ boot_km_dt <- surv_data[,{
#~ 					s <-Surv(lifespan_baseline, dead)
#~ 					o <- Hmisc::bootkm(s, B=1000)
#~ 					out <- quantile(o, p=c(0.05, 0.95))
#~ 					#print(median(s))
#~ 					as.list(out)
#~ 					}, by= "sex,treatment"]
					


summ_dt_global <- all_data_in_sd[, .(global_avg = mean(value)), by="id,behaviour"]
corr_data <- summ_dt_global[surv_data[, list(id, lifespan_baseline)], on="id"]
corr_data <- corr_data[dt[, .(interactions = sum(interactions)), by=id], on = "id"]
corr_data <- na.omit(corr_data[dt[meta=T], on="id"][is.finite(lifespan)])


behav_lifespan <- ggplot(corr_data, aes(global_avg, lifespan, colour=treatment)) + geom_point() +
						facet_grid(behaviour ~ sex) + geom_smooth(method="lm") +
						scale_fill_manual(values=CONTROL_SD_PALETTE) +
						scale_colour_manual(values=CONTROL_SD_PALETTE)
			

mod <- lm(lifespan ~ (global_avg : behaviour)  * sex, corr_data[treatment == "Control"])
summary(mod)

mod <- lm(lifespan ~ (global_avg : behaviour)  * sex, corr_data[treatment == "SD"])
summary(mod)

corr_data[, sex := as.factor(sex)]
rf <- randomForest(lifespan ~ global_avg * behaviour  * sex * treatment, corr_data, ntree=5000)
rf
col = CONTROL_SD_PALETTE[2]
corr_data[, k_interactions := interactions/1000]
interaction_lifespan_pl <- ggplot(corr_data[treatment == "SD"], aes(k_interactions, lifespan)) + geom_point(colour=col) +
								facet_grid(. ~ sex) + geom_smooth(method="lm",colour=col, fill=col) 
			
mod <- lm(lifespan ~  k_interactions, corr_data[treatment == "SD"])
summary(mod)




#~ lm(quiet_rebound_day3h_diff ~quiet_sd_night_diff, stat_rebound_dt[sex=="M"])$coefficients
foo <- function(d, i){
	#lm(quiet_rebound_day3h_diff ~quiet_sd_night_diff, d, subset=i)$coefficients
	lm(lifespan ~ k_interactions, d, subset=i)$coefficients
	}
	
bci <- boot::boot(data=corr_data[treatment == "SD"], statistic=foo, R=1000)




col = CONTROL_SD_PALETTE[2]
ggplot(corr_data[treatment == "SD"], aes(global_avg, lifespan, colour=treatment)) + geom_point(colour=col) +facet_grid(behaviour ~ sex) + geom_smooth(method="lm", colour=col)

mod <- lm(lifespan ~ (global_avg : behaviour)  * sex, corr_data[treatment == "Control"])
summary(mod)



corr_data[, sex := as.factor(sex)]
rf <- randomForest(lifespan ~ global_avg * behaviour  * sex, corr_data[treatment == "Control"])
summary(rf)



pdf("prolonged_sd_correlations.pdf")
ggplot(corr_data[treatment == "Control"], aes(x=overall_sleep, y=lifespan_baseline, colour=sex, shape=sex)) + 
				geom_point() + 
				geom_smooth(method="lm") 
				

corr_dt <- corr_data[treatment == "Control" & sex=="M"]				
cor.test(corr_dt$lifespan_baseline, corr_dt$overall_sleep, method="spearman")
corr_dt <- corr_data[treatment == "Control" & sex=="F"]				
cor.test(corr_dt$lifespan_baseline, corr_dt$overall_sleep, method="spearman")


ggplot(corr_data[treatment != "Control"], aes(x=interactions, y=lifespan_baseline, colour=sex, shape=sex)) + 
				geom_point() + 
				geom_smooth(method="lm") 

corr_dt <- corr_data[treatment != "Control" & sex=="M"]				
cor.test(corr_dt$lifespan_baseline, corr_dt$interactions, method="spearman")
corr_dt <- corr_data[treatment != "Control" & sex=="F"]				
cor.test(corr_dt$lifespan_baseline, corr_dt$interactions, method="spearman")

dev.off()


stim_dt <- dt[t > days(1) &  t < days(10) & xmv(treatment) == "SD"]

stim_dt <-behavr::bin_apply_all(stim_dt, y=interactions, x_bin_length=mins(30), FUN=sum)

stim_simple_dt <- rejoin(stim_dt)[, .(N_stimuli = mean(interactions)), by="t,sex"]

ggplot(stim_simple_dt, aes(x=t, y=N_stimuli, colour=sex)) + geom_line()


tmp_ts <- ts(stim_simple_dt[sex=="M", N_stimuli],  frequency=48) # 48 obervations a day
s <- stl(tmp_ts,s.window="per")
apply(s$time.series,2, var) / var(tmp_ts)
plot(s)
ts <- s$time.series
dt_f <- data.table(interactions = ts[, "trend"], t = time(ts) * days(1), sex="M")

tmp_ts <- ts(stim_simple_dt[sex=="F", N_stimuli],  frequency=48) # 48 obervations a day
s <- stl(tmp_ts,s.window="per")
apply(s$time.series,2, var) / var(tmp_ts)
plot(s)
ts <- s$time.series
dt_m <- data.table(interactions = ts[, "trend"], t = time(ts) * days(1), sex="F")


all_pl_objs$etho_stimuli_2 <- ggetho(dt[xmv(sdi) == 10 & t > days(1) & t < days(10)],
                                    aes(y = interactions, fill=sex),
                                    summary_FUN = sum) +
                                    stat_pop_etho(method= mean_cl_boot, linetype=2) +
									stat_ld_annotations() + scale_y_continuous(limits = c(NA,NA)) +
									scale_fill_manual(values=FEMALE_MALE_PALETTE) +scale_colour_manual(values=FEMALE_MALE_PALETTE) +
									geom_line(data=rbind(dt_f, dt_m),  size=2) +
									ethogram_theme

all_pl_objs$etho_stimuli_2




pdf("prolonged_sd_rebound.pdf", w=12,h=6)
for(p_name in names(all_pl_objs)){
	print(p_name)
    pl <- all_pl_objs[[which(names(all_pl_objs) == p_name)]]
    print(pl + ggtitle(p_name))
}
dev.off()


