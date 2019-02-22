options(nwarnings = 1000)
library(cowplot)
library(scopr)
library(ggetho)
library(sleepr)

#FEMALE_MALE_PALETTE <- c("#be2828ff", "#282896ff")
#FEMALE_MALE_PALETTE <- c("#f57d75ff", "#74a2ceff")
CONTROL_SD_PALETTE <- c( "#969696ff", "#3caa3cff")

                      
   
curate_data <- function(data){
  data[, t := t - days(xmv(baseline_days))]
	# first we remove animals that do not have enought data points
	valid_animals <- data[,.(t_range = max(t) - min(t)), by=id][t_range >= days(5)]$id
	data <- data[t > days(-3) & 
				 t < days(+2) &
				 id %in% valid_animals]
	
	data[, treatment := as.factor(ifelse(sdi == 0, "Control", "SD")), meta=T]
	data[, interval := round(((11-sdi) ^ 1.7)) * 20, meta=T]	
    data[, interval := plyr::revalue(as.factor(interval),c("1180"="Control")), meta=T]
	
	data[,x_rel:=ifelse(xmv(region_id) > 10, 1-x, x)]

	norm_x <- function(x){
		min_x <- quantile(na.omit(x),probs=0.01)
		x <- x - min_x
		max_x <- quantile(na.omit(x),probs=0.99)
		x / max_x
		}
	data[, x_rel := norm_x(x_rel), by=id]
	
	# We currate furter the data by removing individuals that do not have a matching pair 
	# with a different treatment in the same experiment.
	# this means we reduce polution of the data by a large number of controls not necessarily done 
	# at the same time/ same flies..
	males_to_keep <- meta(data)[sex == "M",
						.(
						  n_conditions = length(unique(sdi)),
						  sex = "M"
						  ), 
						by = .(datetime, machine_id)]
	           
	females_to_keep <- meta(data)[sex == "F",
						.(
						  n_conditions = length(unique(sdi)),
						  sex = "F"
						  ), 
						by = .(datetime, machine_id)]
						           
	experiments_to_keep <- rbind(males_to_keep, females_to_keep)
	experiments_to_keep <- experiments_to_keep[n_conditions > 1, -"n_conditions", with=FALSE]
	id_to_keep <- meta(data)[experiments_to_keep, on=c(names(experiments_to_keep))]$id
	out <- data[ id %in% id_to_keep, verbose=T]
	out
}

####### population ethogrames here

# a set of layers or our next big plots
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
	
	


make_plot <- function(){	
	summary(dt)

	dt_etho <- dt[xmv(sdi) %in% c(0,10) & t %between% days(c(-1,2))]

	all_pl_objs <- list()


	all_pl_objs$etho_quiet <- ggetho(dt_etho,
										aes(y = !moving, fill=treatment)) +
										layers() 				+
										scale_y_continuous(limits = c(0,1), name=expression(P(immobile))) 
															
	all_pl_objs$etho_quiet

	all_pl_objs$etho_quiet_zoom <- all_pl_objs$etho_quiet + 
											coord_cartesian(xlim = c(SD_START_END[2]-hours(2),SD_START_END[2] + hours(6)))
	all_pl_objs$etho_quiet_zoom                            
	all_pl_objs$etho_stimuli <- ggetho(dt_etho[t %between% SD_START_END],
										aes(y = interactions, fill=treatment),
										summary_FUN = sum) +
										layers(annotate=FALSE)+
										scale_y_continuous(limits = c(NA,50), name=expression(N[stimuli]~(30~min^{-1}))) 
										 
	all_pl_objs$etho_stimuli
	
	
#~ 	all_pl_objs$etho_stimuli <- ggetho(dt_etho,
#~ 										aes(y = x_rel, fill=treatment)) +
#~ 										layers(annotate=FALSE)+
#~ 										#scale_y_continuous(limits = c(NA,NA), name=expression(N[stimuli]~(30~min^{-1}))) 
	

	################### scalar stats here and barplots (1 value/animal)




dt_ord_2 <- dt[,
			{
			dd <- copy(.SD[, .(t, x_rel, moving)])
			dd[, t := (floor(t/60))*60]
			dd[,
				.(
					x_rel=mean(x_rel),
					walked_dist = sum(abs(diff(x_rel))),
					behaviour = ifelse(all(!moving),1,2)
					),
					by="t"]
			},
			by="id"]


	dt_ord_2[, behaviour := ifelse(behaviour != 1, (walked_dist > .25) + 2  , 1)]
	dt_ord_2[, behaviour:= ordered(c("q","m","w")[behaviour], levels = c("q","m","w"))]
	dt_ord_2[, micro.mov. := (behaviour == "m")]
	dt_ord_2[, walking := (behaviour == "w")]
	dt_ord_2[, quiescent := (behaviour == "q")]


	# shared layers
	stat_rebound_dt <- rejoin(
							dt_ord_2[,
							.(
							quiet_baseline_day3h = mean(quiescent[t %between% REF_SLEEP_WINDOW]),
							quiet_rebound_day3h = mean(quiescent[t %between% TARGET_SLEEP_WINDOW])
							)
							,by = id]
							)
	stat_rebound_dt[, treatment := as.character(ifelse(interval =="Control", "Control", "SD"))]

	
	mod <- lm(quiet_rebound_day3h ~ quiet_baseline_day3h * sex, stat_rebound_dt[sdi == 0])

	stat_rebound_dt[, quiet_rebound_day3h_pred := predict(mod, stat_rebound_dt)]
	stat_rebound_dt[, quiet_rebound_day3h_diff := quiet_rebound_day3h - quiet_rebound_day3h_pred]


	summ_stat <- stat_rebound_dt[,{
					w = wilcox.test(quiet_rebound_day3h_diff, mu=0, alternative="greater")
					list(pval=w$p.value, n=.N)},
					keyby="sex,interval,treatment"]
	print(summ_stat)			
	
	summ_stat[, text:=sprintf("%s\nN=%03d", gtools::stars.pval(pval), n)]

	all_pl_objs$bar_quiet_reb_day3h_diff <- ggplot(stat_rebound_dt, aes(x=treatment, quiet_rebound_day3h_diff * 3 * 60, colour=treatment)) + layer_barpl() +
							geom_hline(yintercept=0, colour="red", linetype=2) + scale_y_continuous(name=expression(Delta~immobility[rebound]~(min)), limits=c(-180, 180)) + 
							geom_label(data= summ_stat, aes(x=treatment, label=text), y=-Inf, inherit.aes=F, vjust = -0.1) +
							scale_x_discrete(name="Treatment")
							


	no_legend = list(theme(legend.position="none"))



tdt <- dt_ord_2[t > SD_START_END[1] - days(1) & t < SD_START_END[2] +days(1)]


l <- list(quiescent ="P(quiescent)",
          micro.mov. ="P(micro-moving)",
	      walking ="P(walking)",
		  x_rel ="Position (rel.)")

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

main_plot <- ggplot(all_data, aes(t, value, fill =treatment, colour=treatment)) + 
				facet_grid(behaviour ~ sex) +
				scale_fill_manual(values=CONTROL_SD_PALETTE) +
				scale_colour_manual(values=CONTROL_SD_PALETTE) +
				stat_pop_etho(method = mean_cl_boot) + 
				scale_y_continuous(limits=c(0,1)) +
				scale_x_days() + theme(legend.position = "none")

pl_row2 <- plot_grid( all_pl_objs$etho_stimuli + no_legend + facet_grid( . ~ sex),
						  all_pl_objs$bar_quiet_reb_day3h_diff + no_legend,
						  labels=LETTERS[2:3],
						  nrow=1)
						  
#~ 	pl_row3 <- plot_grid( all_pl_objs$bar_interactions + no_legend,  
#~ 						  all_pl_objs$bar_quiet_reb_day3h_diff + no_legend,
#~ 						  labels=LETTERS[4:5],
#~ 						  nrow=1)


	pdf(PDF_NAME, w=12, h=8)
 	print(plot_grid(main_plot,
 				pl_row2,
 				labels=c("A", ""), 
 				ncol=1, rel_heights=c(2,1)))
	dev.off()
	return(stat_rebound_dt)
}

met <- fread(METADATA)
met <- met[status == "OK"]
met <- link_ethoscope_metadata(met, result_dir = RESULT_DIR)



dt <- load_ethoscope(met,
					   max_time=days(7),
					   reference_hour=9.0,
					   cache = CACHE,
					   FUN = sleep_annotation,
					   ncores=1)
									
# then we apply this function to our data
dt <- curate_data(dt)
# we have a look at our resulting data
