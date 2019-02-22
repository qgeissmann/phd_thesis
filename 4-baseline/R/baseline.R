rm(list=ls());gc()
options(nwarnings = 1000)
library(scopr)
library(ggetho)
library(sleepr)
library(cowplot)
#source("../ggplot_themes.R")
#~ library(ggtern)
#~ library(grid)
#~ library(Gmisc)

FEMALE_VIDEO_TABLE <- "female_videos.csv"
#FEMALE_MALE_PALETTE <- c("#be2828ff", "#282896ff")
FEMALE_MALE_PALETTE <- c("#f57d75ff", "#74a2ceff")
CONTROL_SD_PALETTE <- c( "#969696ff", "#3caa3cff")
METADATA <- "metadata.csv"
CACHE <- "./cache/"
BEHAVIOUR_STATES_PALETTE <- c("#999999ff", "#4e9a06ff", "#0070b0ff")
#~ RESULT_DIR <- "./raw_results/"
RESULT_DIR <- "/run/media/quentin/Maxtor/auto_generated_data"
#REMOTE_DIR <- "ftp://nas.lab.gilest.ro/auto_generated_data/ethoscope_results/"

# remover 0453|2016−04−04_17−39−52_035aee|17 ?!

THR <- 1
THR_W <- 2.5


met <- fread(METADATA)
met <- met[status == "OK"]

#~ met <- link_ethoscope_metadata_remote(met,
#~                                       remote_dir =  REMOTE_DIR,
#~                                       result_dir = RESULT_DIR,
#~                                       verbose = TRUE)

met_init <- link_ethoscope_metadata(met, result_dir = RESULT_DIR)
                                      
                                                                            
dt <- load_ethoscope(met_init,
					   max_time=days(7),
					   reference_hour=9.0,
					   cache = CACHE,
					   FUN = sleep_annotation,
					   ncores=1)

summary(dt)
                 
                 

   
curate_data <- function(data){
  data[, t := t - days(xmv(baseline_days))]
  data <- data[is_interpolated == F]
	# first we remove animals that do not have enought data points
	valid_animals <- data[,.(t_range = max(t) - min(t)), by=id][t_range >= days(5)]$id
	data <- data[t > days(-2) & 
				 t < days(+2) &
				 id %in% valid_animals]
	data[, t := t+ days(2)]
	data <- data[xmv(sdi)==0]
	data[,x_rel:=ifelse(xmv(region_id) > 10, 1-x, x)]

	norm_x <- function(x){
		min_x <- quantile(na.omit(x),probs=0.01)
		x <- x - min_x
		max_x <- quantile(na.omit(x),probs=0.99)
		x / max_x
		}
	data[, x_rel := norm_x(x_rel), by=id]
}

# then we apply this function to our data
dt <- curate_data(dt)


dt_dam <-  bin_apply_all(dt, beam_crosses, x_bin_length = 60, FUN=sum)
dt_dam[ ,moving:=beam_crosses >0]


dt_dam[,asleep := sleepr:::sleep_contiguous(moving,
                                        1/60,
                                        min_valid_time = 300), by=id]



dt_dam[, h := floor((t %% hours(24)) / hours(.5)) * hours(.5)]
dt[, h := floor((t %% hours(24)) / hours(.5)) * hours(.5)]

stat_dt_dam <- dt_dam[, .(
					mean_moving_d = mean(moving),
					mean_sleep_d = mean(asleep)
					) ,by="id"]
stat_dt <- dt[, .(
					mean_moving = mean(moving),
					mean_sleep = mean(asleep)
					) ,by="id"]

					


stat_dt_dam_h <- dt_dam[, .(
					mean_moving_d = mean(moving),
					mean_sleep_d = mean(asleep)
					) ,by="id,h"]

stat_dt_h <- dt[, .(
					mean_moving = mean(moving),
					mean_sleep = mean(asleep)
					) ,by="id,h"]

#stat_dt_dam_h <- rejoin(stat_dt_dam_h)


stat_dt_dam_h <- stat_dt_dam_h[stat_dt_h, on =c("id", "h")][dt[, meta=T], on="id"]

ggplot(stat_dt_dam_h,aes(mean_sleep_d, mean_sleep, colour=sex)) + geom_point(alpha=.2) + facet_wrap(~h, ncol=6) + 
			scale_x_continuous(limits=c(0,1)) +
			scale_y_continuous(limits=c(0,1)) + 
			geom_abline(slope=1, intercept=0) + geom_smooth(method="lm")
			


rho <- function(x,y){
	cor.test(x, y, method="spearman")$estimate
}
stat_dt_dam_h_summ <-stat_dt_dam_h[, .(
						rho_m = rho(.SD$mean_moving_d, .SD$mean_moving),
						rho_s = rho(.SD$mean_sleep_d, .SD$mean_sleep)
						),by="h,sex"]

stat_dt_dam_h_summ[, dummy := ""]

pl_rho_m_over_h <- ggplot(stat_dt_dam_h_summ, aes(x=h, y=rho_m, colour=sex)) + geom_point() +
						scale_fill_manual(values=FEMALE_MALE_PALETTE)+
						scale_colour_manual(values=FEMALE_MALE_PALETTE) +
				geom_line() + scale_x_hours(limits=c(0,hours(24))) + stat_ld_annotations()  + facet_grid(dummy~.) +
				scale_y_continuous(limits=c(0,1), name= expression(rho[DAM~","~Etho.]))

pl_rho_s_over_h <- ggplot(stat_dt_dam_h_summ, aes(x=h, y=rho_s, colour=sex)) + geom_point() +
						scale_fill_manual(values=FEMALE_MALE_PALETTE)+
						scale_colour_manual(values=FEMALE_MALE_PALETTE) +
				geom_line() + scale_x_hours(limits=c(0,hours(24))) + stat_ld_annotations()  + facet_grid(dummy~.) +
				scale_y_continuous(limits=c(0,1), name= expression(rho[DAM~","~Etho.]))


dmet <- copy(dt_dam[meta=T])
dt_dam_etho <- copy(dt_dam)
dmet[ ,id := paste("d_", id, sep="")]
dt_dam_etho[ ,id := paste("d_", id, sep="")]
setkeyv(dt_dam_etho, "id")
setkeyv(dmet, "id")
dt_dam_etho <- setbehavr(dt_dam_etho,dmet)

dt_dam_etho <- bind_behavr_list(l=list(dt_dam_etho, dt[, .(id,t, beam_crosses, moving, asleep,h)]))
dt_dam_etho[, method := ifelse(id %like% "d_","DAM", "Ethoscope" ),meta=T]

#ggetho(dt_dam, aes(y=!moving, colour=sex)) + stat_pop_etho() + facet_grid( method ~ .)

#~ ggetho(dt_dam, aes(y=!moving, colour=sex)) + stat_pop_etho() + facet_grid( method ~ .) +
#~ 			scale_fill_manual(values=FEMALE_MALE_PALETTE)+
#~ 			scale_colour_manual(values=FEMALE_MALE_PALETTE)
			

pl_activity_zt <- ggetho(dt_dam_etho, aes(y=!moving, colour=sex), time_wrap=hours(24)) + 
						stat_pop_etho(method=mean_cl_boot) + 
						facet_grid( method ~ .) +
						scale_fill_manual(values=FEMALE_MALE_PALETTE)+
						scale_colour_manual(values=FEMALE_MALE_PALETTE) +
						scale_y_continuous(limits=c(0,1), name="P(immobile)")
			
pl_sleep_zt <- ggetho(dt_dam_etho, aes(y=asleep, colour=sex), time_wrap=hours(24)) + stat_pop_etho(method=mean_cl_boot) + 
						facet_grid( method ~ .) +
						scale_fill_manual(values=FEMALE_MALE_PALETTE)+
						scale_colour_manual(values=FEMALE_MALE_PALETTE) +
						scale_y_continuous(limits=c(0,1), name="P(asleep)")
						


bdt_m <- bout_analysis(moving, dt_dam_etho)
bdt_s <- bout_analysis(asleep, dt_dam_etho)

bdt_m[, phase := ifelse(t %% hours(24) < hours(12), "L", "D")]
bdt_s[, phase := ifelse(t %% hours(24) < hours(12), "L", "D")]

bdt_m[, phase := factor(phase, levels =c("L", "D"))]
bdt_s[, phase := factor(phase, levels =c("L", "D"))]

summ_sleep_bout <-  bdt_s[meta=T][bdt_s[asleep==T][, .(mean_length=mean(duration),n=.N), by="id,phase"]]
summ_immob_bout <-  bdt_m[meta=T][bdt_m[moving==F ][, .(mean_length=mean(duration),n=.N), by="id,phase"]]


summ_sleep_bout2 <-  bdt_s[meta=T][bdt_s[asleep==T][, .(mean_length=mean(duration),n=.N), by="id"]]
summ_immob_bout2 <-  bdt_m[meta=T][bdt_m[moving==F ][, .(mean_length=mean(duration),n=.N), by="id"]]


#~ ggplot(summ_immob_bout, aes(mean_length, n, colour=sex)) + geom_point() + facet_grid(method ~ .) + scale_x_log10() + scale_y_sqrt()
#~ ggplot(summ_sleep_bout, aes(mean_length, n, colour=sex)) + geom_point() + facet_grid(method ~ .) + scale_x_log10() + scale_y_sqrt()


pl_sleep_fract <-  ggplot(summ_sleep_bout, aes(mean_length, n, colour=sex)) + 
															  geom_point(alpha=.3, size=.3) + 
															  facet_grid(method ~ phase) + 
															  scale_x_log10(name="Average bout duration (s)") + 
															  scale_y_sqrt(name=expression(N[bout]), breaks=c(5, 20, 50, 100, 200)) + 
															  stat_density2d(bins=4, size=.5)
															  
pl_immob_fract <- ggplot(summ_immob_bout, aes(mean_length, n, colour=sex)) + 
															geom_point(alpha=.3, size=.3) + 
															facet_grid(method ~ phase) + 
															scale_x_log10(name="Average bout duration (s)") + 
															scale_y_sqrt(name=expression(N[bout]), breaks=c(50, 200, 500, 1000, 2000)) + 
															stat_density2d(bins=4, size=.5)
 
left <- plot_grid( pl_activity_zt + theme(legend.position="none"),
			   pl_rho_m_over_h + theme(legend.position="none"),
			   pl_immob_fract  + theme(legend.position="none"),
			   nrow=3, labels=LETTERS[1:3], rel_heights=c(2,1,2))

right <- plot_grid( pl_sleep_zt  + theme(legend.position="none"),
			   pl_rho_s_over_h + theme(legend.position="none"),
			   pl_sleep_fract + theme(legend.position="none"),
			   nrow=3, labels=LETTERS[4:6], rel_heights=c(2,1,2))
		
pl <- plot_grid(left, right, ncol=2)

pdf("activity_sleep_dam_vs_etho.pdf", w=12, h=8)
print(pl)
dev.off()


#ggetho(bdt_m[moving==F], aes(x=t, y= duration, colour=sex), summary_FUN=length, time_wrap=hours(24)) + stat_pop_etho() + facet_grid(method~.) + scale_y_log10()
############################## second order behaviour state	


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

# test
#~ pdf("walking_threshold_validation.pdf", w=6, h=6)
#~ ggplot(dt[behaviour != "q" & xmv(sex) == "F"], aes(walked_dist * 60)) + 
#~ 			geom_histogram(aes(y=100 * ..count../sum(..count..)), bins=200) + 
#~ 			scale_x_sqrt(limits = 60 * c(0, 3), name=expression(Total~distance~moved~(mm.min^{-1})), breaks= 60 * c(0.1, 0.2, 0.5,1,2)) +
#~ 			scale_y_continuous(name= "Density (%)") +
#~ 			geom_vline(xintercept= 60 *.25, size=2, colour="red") 
#~ dev.off()




set.seed(20)
rep_ids <-	sapply( c("M", "F"), 
				function(s){
					sample(dt_ord_2[meta=T][ sex== s, id],
						   size=3, 
						   replace=F)
						   }
			   )
	   
dt_2_rep <- copy(dt_ord_2[id %in% rep_ids & t %between% days(c(2,4))])
dt_full_rep <- copy(dt[id %in% rep_ids & t %between% days(c(2,4))])
dt_2_rep[, label := as.factor(as.numeric(id))]
dt_2_rep[, sex := xmv(sex)]
pl_exple_behaviour  <- ggplot(dt_2_rep, aes(t - days(2) ,x_rel)) + 
							geom_rect(mapping=aes(xmin=(t - days(2)) , xmax=(60+t- days(2)),fill=behaviour), ymin=-1, ymax=1, alpha=.90)+
							#geom_line(data=dt_full_rep, size = .1) +
							geom_line(size = .5) +
							scale_fill_manual(values=BEHAVIOUR_STATES_PALETTE) +
							geom_hline(yintercept=.5, linetype=2, alpha=.5, size=2) +
							facet_grid(sex + label ~ .)  + 
							scale_x_hours() + scale_y_continuous(breaks=c(0,  0.5), name="Position (rel.)") + 
							theme(panel.spacing = unit(0, "cm"), strip.text.y = element_text(angle = -0)) +
							coord_cartesian(ylim = c(0,1))

			

dt_ord_2[, sex := xmv(sex)]
n_bins <- 50
dt_ord_2[, x_rel_cut := floor(x_rel * n_bins) / n_bins]
dt_ord_2[, x_rel_cut := ifelse(x_rel_cut > 1, 1, x_rel_cut)]
dt_ord_2[, x_rel_cut := ifelse(x_rel_cut <0, 0, x_rel_cut)]



hist_pos_dt  <- dt_ord_2[,.(n=.N),  by="sex,behaviour,x_rel_cut"]
hist_pos_dt[,N := sum(n),  by="sex"]
hist_pos_dt[,density := n_bins * n/N]
pl_pos_density <- ggplot(hist_pos_dt, aes(x_rel_cut, density , fill=behaviour)) + 
						geom_density(stat="identity", alpha=.33) + 
						facet_grid(sex ~ .) +
						scale_fill_manual(values=BEHAVIOUR_STATES_PALETTE, name="Behaviour")  +
						scale_x_continuous(name = "Position (rel.)")  +
						scale_y_continuous(name = "P(behaviour, position)")  
						
						
			

			

lev <- dt_ord_2[, levels(behaviour)]
l <- sapply(lev, function(x)dt_ord_2[,behaviour == x])
dt_ord_2 <- cbind(dt_ord_2,l)
setbehavr(dt_ord_2, dt[meta=T])
layers <- function(colour){
	
		list(
			stat_pop_etho(colour=colour),
			facet_grid(. ~ sex),
			scale_y_continuous(limits=c(0,.8)),
			scale_x_hours()
			)
		}
pl_tmpl <- ggetho(dt_ord_2, aes(y=q), time_wrap = hours(24)) + 
			facet_grid(sex ~ .)+
			scale_y_continuous(limits=c(0,1), name="P(behaviour | time)")+
			scale_x_hours() + stat_ld_annotations()

pl_data_q <- ggetho(dt_ord_2, aes(y=q), time_wrap = hours(24))$data
pl_data_m <- ggetho(dt_ord_2, aes(y=m), time_wrap = hours(24))$data
pl_data_w <- ggetho(dt_ord_2, aes(y=w), time_wrap = hours(24))$data

pl_behav_vs_t <- pl_tmpl + 	stat_pop_etho(mapping = aes( y = q), data=pl_data_q, colour = BEHAVIOUR_STATES_PALETTE[1], method=mean_cl_boot) +
								stat_pop_etho(mapping = aes( y = m), data=pl_data_m, colour = BEHAVIOUR_STATES_PALETTE[2], method=mean_cl_boot) +
								stat_pop_etho(mapping = aes( y = w), data=pl_data_w, colour = BEHAVIOUR_STATES_PALETTE[3], method=mean_cl_boot)


pl2 <- plot_grid(	pl_exple_behaviour +theme(legend.position="bottom"),
					plot_grid(pl_behav_vs_t, pl_pos_density + theme(legend.position="none"),  labels=LETTERS[2:3], nrow=1, rel_widths = c(2,1)),
			nrow=2, labels=LETTERS[1])
			
pdf("behavioural_state.pdf", w=12, h=8)
print(pl2)
dev.off()
##########
tern_dt <- dt_ord_2[,
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
			by="id"]
			
			
stat_dt_ord2 <- dt_ord_2[, .(q=mean(q),m=mean(m), w=mean(w)), by=id]

tern_dt[, behaviour := factor(behaviour, levels = c("q","m","w","x_rel"))]

tern_dt_wide <- behavr(data.table(reshape(tern_dt[!(behaviour  %in% c("x_rel", "entropy"))], 
							timevar = "behaviour", 
							idvar = c("id", "hour"),
							direction = "wide"), key="id"),
					   	tern_dt[meta=T])
setnames(tern_dt_wide, c("value.q", "value.m", "value.w"),c("q","m","w"))					   	
tern_dt_wide <- rejoin(tern_dt_wide)

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


	
set.seed(2)
clust_dt <- copy(tern_dt_wide[, .(id,hour, q,m,w)])
clust_dt <- na.omit(clust_dt)
setkeyv(clust_dt, c("id", "hour"))

valid_ids <- c( as.character(sample(dt[sex=="M",id, meta=T], 400, replace=F)),
				as.character(sample(dt[sex=="F",id, meta=T], 400, replace=F)))
				
clust_dt <- clust_dt[id %in% valid_ids]
list_dt <- split(clust_dt, as.character(clust_dt$id))
system.time(distances <- proxy::dist(list_dt, method=path_divergeance))
hcl <- hclust(distances, method="average")


#fixme
#dt[,mean_asleep:= runif(.N), meta=T]





library(ggdendro)
library(dendextend)
library(ape)
#pdf("unrooted_dendro.pdf")
dend <- as.dendrogram(hcl)

names(FEMALE_MALE_PALETTE) <- c("F", "M")

colbranches <- function(n){
	tmp_labs <- labels(n)
#	print(tmp_labs)
	new_col <- unique(FEMALE_MALE_PALETTE[dt[meta=T][tmp_labs,sex]])
	print(new_col)
    if(length(new_col) >1)
        return(n)
  a <- attributes(n) # Find the attributes of current node
  attr(n, "edgePar") <- c(a$edgePar, list(col=new_col, lwd=2))
  n
  }

dend2 <- dendrapply(dend, colbranches)

N <- length(valid_ids)
l=numeric(N*2)
l[1] <- "black"
i=2
a <- dendrapply(dend2, function(x){
    o = attributes(x)$edgePar$col
    print(o)
    o <- ifelse(is.null(o), "black", o)
    l[i] <<- o
    i <<- i+1
    })


clust_res <- letters[1:4][stats::cutree(hcl, k=4)]
tree <- as.phylo(dend2)
labs <- tree$tip.label
tree$tip.label <- clust_res
cols <- FEMALE_MALE_PALETTE[dt[meta=T][labs,sex]]

pdf("unrooted_dendro.pdf")
plot(tree, show.tip.label=T, lab4ut="axial", tip.color = cols , edge.color=l[3:(N*2)], cex=1.5, type="unrooted")
add.scale.bar(cex = 1.5, font = 2)
dev.off()



stat_dt_ord2 <- dt_ord_2[, .(q=mean(q),m=mean(m), w=mean(w),
								m_12_21 = mean(m[t %% hours(24) %between% hours(c(12,21))]),
								m_12_21_by_food = mean(x_rel[t %% hours(24) %between% hours(c(12,21))] < .06 &  m[t %% hours(24) %between% hours(c(12,21))]),
								q_21_24 = mean(q[t %% hours(24) %between% hours(c(21,24))])
								), by=id]


stat_dt <- dt[, .(
					mean_moving = mean(moving),
					mean_sleep = mean(asleep),
					mean_sleep_21_24 = mean(asleep[t %% hours(24) > hours(21)])
					) ,by="id"]


stat_dt_ord2 <- rejoin(stat_dt)[stat_dt_ord2]

stat_dt_ord2[, round(mean_cl_boot(mean_sleep) * 60*24,2),by="sex"]


stat_dt_ord2[, round(mean_cl_boot(m_12_21 * (21 - 12)/ (m *24)) * 100,2),by="sex"]
stat_dt_ord2[, round(mean_cl_boot(100*q_21_24 * 3)/(q*24),2),by="sex"]
stat_dt_ord2[, round(mean_cl_boot(m_12_21_by_food/m_12_21) * 100,2),by="sex"]
100 * sum(stat_dt_ord2[sex=="F",mean_sleep] < 0.05) / nrow(stat_dt_ord2[sex=="F"])






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
dt_ord_2[,sex :=xmv(sex)]

pdf("walked_dist_density.pdf", w=12, h=8)
pl <- ggplot(dt_ord_2[quiescent == F], aes(walked_dist * 60,colour=sex, fill=sex)) + 
		geom_density(alpha=.3) +
		scale_x_sqrt(limits = c(0, 300), breaks = c(10, 20, 50, 100, 200)) + 
		geom_vline(xintercept=15) 
pl 

dt_ord_2[, ZT := floor((t / hours(1)) %% 24)]
pl + facet_wrap( ~ ZT, nrow=6)		
		
dev.off()
