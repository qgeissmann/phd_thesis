rm(list=ls());gc()
options(nwarnings = 1000)
library(scopr)
library(ggetho)
library(sleepr)
library(gtools)
library(magrittr)
library(corrr)
library(dplyr)
source("../ggplot_themes.R")


#~ FEMALE_MALE_PALETTE <- c("#be2828ff", "#282896ff")
FEMALE_MALE_PALETTE <- c("#f57d75ff", "#74a2ceff")
CONTROL_SD_PALETTE <- c( "#969696ff", "#3caa3cff")
METADATA <- "metadata_complete.csv"
CACHE <- "./cache/"
RESULT_DIR <- "./raw_results/"
REMOTE_DIR <- "ftp://nas.lab.gilest.ro/auto_generated_data/ethoscope_results/"
BEHAVIOUR_STATES_PALETTE <- c("#999999ff", "#4e9a06ff", "#0070b0ff")

met <- fread(METADATA)
met <- met[status == "OK"]
met <- link_ethoscope_metadata_remote(met,
                                      remote_dir =  REMOTE_DIR,
                                      result_dir = RESULT_DIR,
                                      verbose = TRUE)
                                      
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
	data <- data[id %in% dt[,.(max(t), max(t)/.N),by=id][V1 > days(13) & V2 < 11][,id]]
	data <- data[t %between% c(days(1), days(13))]
	data
}

# then we apply this function to our data
dt <- curate_data(dt)

pl <- ggetho(dt, aes(y=paste(as.Date(datetime),  id), z=asleep)) + 
			stat_tile_etho() 


dt[, LD := as.factor(ifelse(t %% hours(24) > hours(12), "D", "L"))]

stat_dt <- dt[, 
				.(sleep_days_02_07 = mean(asleep[t>days(2) & t < days(7)]),
				  sleep_days_08_13 = mean(asleep[t>days(8) & t < days(13)]),
				  sleep_days_11_13 = mean(asleep[t>days(11) & t < days(13)]),
				  sleep_days_08_10 = mean(asleep[t>days(8) & t < days(10)]),
				  sleep_days_05_07 = mean(asleep[t>days(5) & t < days(7)]),
				  sleep_days_02_04 = mean(asleep[t>days(2) & t < days(4)]),
				  quiet_days_02_07 = mean(!moving[t>days(2) & t < days(7)]),
				  quiet_days_08_13 = mean(!moving[t>days(8) & t < days(13)]),
				  quiet_days_11_13 = mean(!moving[t>days(11) & t < days(13)]),
				  quiet_days_08_10 = mean(!moving[t>days(8) & t < days(10)]),
				  quiet_days_05_07 = mean(!moving[t>days(5) & t < days(7)]),
				  quiet_days_02_04 = mean(!moving[t>days(2) & t < days(4)])
				  ),
				  by = c("id")]
				  
stat_dt <- rejoin(stat_dt)


pdf("tube_change_correlations.pdf", w=12,h=6)
# female correlation mat

as.data.table(select(stat_dt[sex=="F"], starts_with("sleep_")) %>% correlate(method = "spearman"))
as.data.table(select(stat_dt[sex=="M"], starts_with("sleep_")) %>% correlate(method = "spearman"))
as.data.table(select(stat_dt[sex=="F"], starts_with("quiet_")) %>% correlate(method = "spearman"))
as.data.table(select(stat_dt[sex=="M"], starts_with("quiet_")) %>% correlate(method = "spearman"))



		 
ggplot(stat_dt, aes(sleep_days_02_07, sleep_days_08_13, colour = sex, shape=sex, fill=sex)) +
			geom_smooth(method="lm")+
			geom_point(size=2, alpha=.5) + labs(x="Sleep before tube change", y="Sleep after tube change") +
			generic_theme
			
		 
ggplot(stat_dt, aes(quiet_days_02_07, quiet_days_08_13, colour = sex, shape=sex, fill=sex)) +
			geom_smooth(method="lm")+
			geom_point(size=2, alpha=.5) + labs(x="Quiescence before tube change", y="Quiescence after tube change") +
			generic_theme

ggplot(stat_dt, aes(quiet_days_02_04, quiet_days_08_13, colour = sex, shape=sex, fill=sex)) +
			geom_smooth(method="lm")+
			geom_point(size=2, alpha=.5) + labs(x="Quiescence before tube change", y="Quiescence after tube change") +
			generic_theme			
			
#~ 			scale_fill_manual(values=FEMALE_MALE_PALETTE) +
#~ 			scale_colour_manual(values=FEMALE_MALE_PALETTE)			
			
dev.off()


################# behavioural paths


dt[,x_rel:=ifelse(xmv(region_id) > 10, 1-x, x)]

norm_x <- function(x){
	min_x <- quantile(na.omit(x),probs=0.01)
	x <- x - min_x
	max_x <- quantile(na.omit(x),probs=0.99)
	x / max_x
	}
dt[, x_rel := norm_x(x_rel), by=id]

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




#~ dt_ord_2[, befr :=  ifelse(t < days(7) | t > days(8), before, NA)]

dt_ord_2[, day := floor(t / days(1))]


tern_dt <- na.omit(dt_ord_2)[,
			{
			dd <- copy(.SD[, .(t, x_rel,behaviour)])
			dd[, hour :=  (floor(t/hours(.25)) * .25) %% 24]
			dd[,
				.(
					value = c(as.vector(table(.SD[,behaviour])/.N), mean(x_rel)),
					behaviour = c(levels(behaviour),"x_rel")
					),
					by="hour"]
			},
			by="id,day"]
			
			

tern_dt[, behaviour := factor(behaviour, levels = c("q","m","w","x_rel"))]
setkeyv(tern_dt,"id")
setbehavr(tern_dt, dt_ord_2[meta=T])




tern_dt_wide <- behavr(data.table(reshape(tern_dt[!(behaviour  %in% c("x_rel", "entropy"))], 
							timevar = "behaviour", 
							idvar = c("id", "hour", "day"),
							direction = "wide"), key="id"),
					   	tern_dt[meta=T])

setnames(tern_dt_wide, c("value.q", "value.m", "value.w"),c("q","m","w"))					   	
tern_dt_wide <- rejoin(tern_dt_wide)
					   	

day_set_map <- c("before_1", "before_2","after_1", "after_2")[c(1,2,1,2,1,2,NA,4,3,4,3,4,3)]

tern_dt_wide[, set := day_set_map[day]]



#~ tern_dt_before <- tern_dt[befr==T]
#~ setkeyv(tern_dt_before,"id")
#~ setbehavr(tern_dt_before, dt_ord_2[meta=T])
#~ tern_dt_before[, befr := NULL]

#~ tern_dt_after <- tern_dt[befr==F]
#~ setkeyv(tern_dt_after,"id")
#~ setbehavr(tern_dt_after, dt_ord_2[meta=T])
#~ tern_dt_after[, befr := NULL]


#~ tern_dt_wide_before <- behavr(data.table(reshape(tern_dt_before[!(behaviour  %in% c("x_rel", "entropy"))], 
#~ 							timevar = "behaviour", 
#~ 							idvar = c("id", "hour"),
#~ 							direction = "wide"), key="id"),
#~ 					   	tern_dt_before[meta=T])
					   	
#~ tern_dt_wide_after <- behavr(data.table(reshape(tern_dt_after[behaviour != "x_rel"], 
#~ 							timevar = "behaviour", 
#~ 							idvar = c("id", "hour"),
#~ 							direction = "wide"), key="id"),
#~ 					   	tern_dt_after[meta=T])
					   	

#~ setnames(tern_dt_wide_before, c("value.q", "value.m", "value.w"),c("q","m","w"))					   	
#~ setnames(tern_dt_wide_after, c("value.q", "value.m", "value.w"),c("q","m","w"))					   	


#~ tern_dt_wide_after <- rejoin(tern_dt_wide_after)
#~ tern_dt_wide_before <- rejoin(tern_dt_wide_before)

### clustering:


path_divergeance <- function(d1, d2){
	valid_t <- intersect(d1$hour, d2$hour)
	d1 <- d1[hour %in% valid_t][,.(q,m,w),keyby=hour][, -"hour"]
	d2 <- d2[hour %in% valid_t][,.(q,m,w),keyby=hour][, -"hour"]
	a <- (rowSums(abs(d1 - d2)))
	bc <- rowSums(sqrt(d1 * d2))
	bc <- ifelse(bc == 0, 1, bc)
	bd <- -log(bc)
	return(mean(bd))
	}


tern_dt_wide <- tern_dt_wide[,.(
								q = mean(q),
								m = mean(m),
								w = mean(w)
								),by="id,set,hour"] 
								

tern_dt_wide <- dt[,.(id,sex),meta=T][tern_dt_wide]


set.seed(3212)

make_dist_data <- function(s, n_sim =5000){
	tdt <- tern_dt_wide[sex==s]
	uids <- unique(tdt$id)

	print(s)

	self_diff_set_dists <- function(set_1, set_2){
								print(c(set_1, set_2))
								sapply(uids, function(i){
									path_divergeance(tdt[set == set_1 & id==i], tdt[set == set_2 & id==i])
								})
						}

	out <- list(
			b1_b2  = self_diff_set_dists("before_1", "before_2"),
			a1_a2  = self_diff_set_dists("after_1", "after_2"),
			a1_b1  = self_diff_set_dists("after_1", "before_1"),
			a2_b2  = self_diff_set_dists("after_2", "before_2")
		)

	samples <-  matrix(sample(uids, replace=T, size=n_sim * 2), nrow=2)

	out$h_a1_b1 <- mapply(FUN =function(i, j){
						path_divergeance(tdt[set == "before_1" & id==i], 
										 tdt[set == "after_1" & id==j])
					}, i = samples[1,], j= samples[2,]
				)		
#~ 	out$h_a2_b2 <- mapply(FUN =function(i, j){
#~ 						path_divergeance(tdt[set == "before_2" & id==i], 
#~ 										 tdt[set == "after_2" & id==j])
#~ 					}, i = samples[1,], j= samples[2,]
#~ 				)		


#~ 	out$h_b1_b2 <- mapply(FUN =function(i, j){
#~ 						path_divergeance(tdt[set == "before_1" & id==i], 
#~ 										 tdt[set == "before_2" & id==j])
#~ 					}, i = samples[1,], j= samples[2,]
#~ 				)		
	out$h_a1_a2 <- mapply(FUN =function(i, j){
						path_divergeance(tdt[set == "after_1" & id==i], 
										 tdt[set == "after_2" & id==j])
					}, i = samples[1,], j= samples[2,]
				)		




	out2 <- rbindlist(
		lapply(names(out), function(n){
			data.table(distance = out[[n]],
						set = n,
						sex=s)
			})
	)
	out2
}

m <- make_dist_data(s="M")
f <- make_dist_data(s="F")
d <- rbind(m, f)

d[ ,group := ifelse(set %like% "h_", "expected", "observed")]
d[ ,group := factor(group,level=c("observed", "expected"))]
d[, set := gsub("h_", "", set)]


pl_dst <- ggplot(d[set %like% "a1_"], aes(y=distance, x= interaction(group, set), fill=sex, alpha=group, linetype=group)) + 
					geom_violin() +
					#geom_jitter(data=d[group=="observed"]) +
					stat_summary(fun.y = "mean", geom="point", shape = 4, size=1, colour="black")+
					stat_summary(fun.data = "mean_cl_boot", geom="errorbar",colour="black", width=0.2)+
					facet_grid(sex ~.) + 
					scale_y_log10(breaks=c(.02, .05, .1,.2, .5), name = "Distance")  +
					scale_x_discrete(name="",
									labels=c("observed.a1_a2" = bquote(A[1]*','*A[2]),
											 "expected.a1_a2" = bquote(H[0](A[1]*','*A[2])),
											 "observed.a1_b1" = bquote(A[1]*','*B[1]),
											 "expected.a1_b1" = bquote(H[0](A[1]*','*B[1]))
											 )
										 )+
					 scale_alpha_manual(values=c(1,0.5))

					 scale_fill_manual()
					 values=c("grey40", "grey70"))


library(cowplot)
library(boot)

boot_dt <- d[d[, .(expected = mean(distance)), by="group,sex,set"][group=="expected", -"group"], on=c("sex","set")][group=="observed"]

boot_dt[, diff_dist := distance - expected]

mean_fun <- function(data, i){
  -mean(boot_dt[i, diff_dist])
}

bo <- boot_dt[, {
	bo <- boot(.SD, statistic=mean_fun, R=5000)
	out <- as.data.table(boot.ci(bo, conf=0.95, type="basic")$basic)
	out
}
	, by="sex,set"]

bo[boot_dt[, .(dist=-mean(diff_dist)), by="sex,set"], on=c("sex","set")]


bo_res <- bo[boot_dt[, .(dist=-mean(diff_dist)), by="sex,set"], on=c("sex","set")]


pl_a <- ggetho(dt, aes(y=!moving, fill=sex))+ 
			stat_pop_etho(method= mean_cl_boot) + 
			scale_y_continuous(name = "P(immobile)") +
			geom_rect(data= data.table(xmin=days(7), xmax=days(8), ymin=0, ymax=.9), 
					  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), size=3, alpha=.5, inherit.aes=F) + 
		  stat_ld_annotations()  + facet_grid(status ~ .) +
		  scale_x_days(limits=c(days(1),days(13)),  breaks=days(1:6) * 2) 
			

tern_dt[, set := day_set_map[day]]


tern_dt

x_given_q_dt <- dt_ord_2[quiescent==T, .(behaviour="x_given_q",value = mean(x_rel)), by="id,day"]

day_dt <- tern_dt[,.(value = mean(value)), by="id,day,behaviour"]
day_lst_dt <- split(day_dt, by="behaviour")
day_lst_dt$ x_given_q <- x_given_q_dt
day_dt <- rbindlist(day_lst_dt)

day_dt <- dt[meta=T][day_dt, on="id"]
ref <- day_dt[day==1, .(id=id, ref=value, behaviour=behaviour)]
day_dt <- na.omit(ref[day_dt, on=c("id","behaviour")])


to_boot <- function(data, i){
	rho <- cor.test(data[i]$value,data[i]$ref, method="spearman", exact=FALSE)$estimate
	}

foo <- function(data){
	
	rho <- to_boot(data, 1:nrow(data))
	bo <- boot(data, statistic=to_boot, R=1000)
	out <- as.data.table(boot.ci(bo, conf=0.95, type="basic")$basic)
	out[, rho := rho]
	}
	
data.table(day=1, rho=1)
all_corr_dt <- day_dt[day > 1 & day < 13,foo(.SD), by="day,sex,behaviour"]
day_1_cor <- data.table(day=1, conf = NA, V2=NA, V3=NA, V4=1, V5=1, rho=1)
all_corr_dt <- all_corr_dt[,rbind(.SD, day_1_cor), by="sex,behaviour"]

pl_corr <- ggplot(all_corr_dt[behaviour != "x_rel"], aes(day * days(1) + hours(12), rho, colour = sex, fill=sex)) + 
				#geom_ribbon(aes(ymin=V4, ymax=V5), alpha=.3, colour=NA ) +
				#all_corr_dt
				geom_bar(stat="identity", position="dodge") +
				geom_errorbar(aes(ymin=V4, ymax=V5), alpha=.3, position="dodge", colour="black") + 
				geom_hline(yintercept=0, linetype = 2) +
				facet_grid(behaviour ~ .) +
				scale_y_continuous(name=expression(rho)) +
				scale_x_days(limits=c(days(1),days(13)),  breaks=days(1:6) * 2) +
				scale_fill_manual(values=FEMALE_MALE_PALETTE) +
				scale_colour_manual(values=FEMALE_MALE_PALETTE)
				
#pl_corr
			 





pl_q <- ggetho(dt_ord_2, aes(y=quiescent))
pl_m <- ggetho(dt_ord_2, aes(y=micro.mov.))
pl_w <- ggetho(dt_ord_2, aes(y=walking))

setnames(pl_q$data, "quiescent", "value")
setnames(pl_m$data, "micro.mov.", "value")
setnames(pl_w$data, "walking", "value")

pl_q$data[, behaviour := "q"]
pl_m$data[, behaviour := "m"]
pl_w$data[, behaviour := "w"]

dt_plot<- rbindlist(list(pl_q$data, pl_m$data, pl_w$data))

pl_behav <- ggplot(dt_plot, aes(y = value, x=t, fill=behaviour, colour=behaviour)) +
				facet_grid(sex ~ .) +
				stat_pop_etho(method= mean_cl_boot) + 
				scale_y_continuous(name = "P(behaviour)") +
				geom_rect(data= data.table(xmin=days(7), xmax=days(8), ymin=0, ymax=.9), 
						  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), size=3, alpha=.5, inherit.aes=F) + 
				stat_ld_annotations()  +  
				scale_x_days(limits=c(days(1),days(13)),  breaks=days(1:6) * 2) 



pl_out <- plot_grid(pl_behav + theme(legend.position="None"), 
					pl_corr  + theme(legend.position="None"), 
					nrow=2, labels = LETTERS[1:2], rel_heights = c(2,4))

pdf("tube_change_correlation.pdf", w=12, h=8)
print(pl_out)
dev.off()

#~ dt_scal <- na.omit(tern_dt)[, .(value= mean(value)), by="set,behaviour,id"]

#~ dt_scal_wide <- reshape(dt_scal,timevar="set",idvar=c("id","behaviour"), direction="wide")
#~ dt_scal_wide <- dt_scal_wide[dt[meta=T], on="id"]




#~ tdt <- copy(dt_scal_wide)
#~ tdt2 <- copy(dt_scal_wide)
#~ tdt[, y:=value.before_2]
#~ tdt[, group:="before_2"]

#~ tdt2[, y:=value.after_1]
#~ tdt2[, group:="after_1"]
#~ tdt <- rbind(tdt,tdt2)

#~ layers <- c(geom_point(),
#~ geom_smooth(method="lm"),
#~ facet_grid(group ~ behaviour),
#~ scale_x_continuous(limits=c(0,1)),
#~ scale_y_continuous(limits=c(0,1)))

#~ pl_corr <- ggplot(tdt, aes(value.before_1, y, colour=sex)) + layers

#~ pl_corr <-  plot_grid(pl_b1_b2, pl_b1_a1, ncol=2)

#~ plot_grid(pl_a, pl_dst, rel_heights = c(1,2), nrow=2)
#~ plot_grid(plot_grid(pl_a, pl_dst, rel_heights = c(1,2), nrow=2),
#~ 	pl_corr, nrow=2,rel_widths=c(4,2)
#~ )
tern_dt
