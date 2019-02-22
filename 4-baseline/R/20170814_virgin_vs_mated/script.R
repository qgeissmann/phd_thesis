rm(list=ls());gc()
source("../trunk.R")
library(scopr)
library(ggetho)
library(sleepr)


METADATA <- "metadata.csv"
CACHE <- "./cache/"
RESULT_DIR <- "/data/ethoscope_results"
SD_START_END <- days(1:2/2)
REF_SLEEP_WINDOW <- hours(c(0, 3))
TARGET_SLEEP_WINDOW <- REF_SLEEP_WINDOW + days(1)
PDF_NAME <- "virgin_mated_states.pdf"
CONTROL_SD_PALETTE <- c( "#969696ff", "#3caa3cff")
                        
   
curate_data <- function(data){
  data[, t := t - days(xmv(baseline_days))]
	# first we remove animals that do not have enought data points
	valid_animals <- data[,.(t_range = max(t) - min(t)), by=id][t_range >= days(5)]$id
	data <- data[t > days(-2.5) & 
				 t < days(+3) &
				 id %in% valid_animals]
	
	data[,x_rel:=ifelse(xmv(region_id) > 10, 1-x, x)]

	norm_x <- function(x){
		min_x <- quantile(na.omit(x),probs=0.01)
		x <- x - min_x
		max_x <- quantile(na.omit(x),probs=0.99)
		x / max_x
		}
	data[, x_rel := norm_x(x_rel), by=id]
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
dt[, treatment := ifelse(mated, "Mated", "Control"), meta=T]

#~ dt_etho <- dt[xmv(sdi) %in% c(0,10) & t %between% days(c(-1,2))]



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



tdt <- dt_ord_2[t > days(-1.5) & t < days(2)]


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
				facet_grid(behaviour ~ .) +
				scale_fill_manual(values=CONTROL_SD_PALETTE) +
				scale_colour_manual(values=CONTROL_SD_PALETTE) +
				stat_pop_etho(method = mean_cl_boot) + 
				scale_y_continuous(limits=c(0,1)) +
				scale_x_hours() + theme(legend.position = "none")

pdf("mated_females_plot.pdf", w=12,h=6)

plot_grid(plot_grid(main_plot, NULL, ncol=1, labels=LETTERS[1:2]),
	     NULL,
	     nrow=1,
	     rel_widths=c(2,1), 
	     labels=c("","C"))

dev.off()
