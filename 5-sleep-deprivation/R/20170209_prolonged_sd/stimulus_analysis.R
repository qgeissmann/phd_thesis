rm(list=ls());gc()
options(nwarnings = 1000)
library(scopr)
library(ggetho)
library(sleepr)
library(gtools)

source("../ggplot_themes.R")

FEMALE_MALE_PALETTE <- c("#be2828ff", "#282896ff")
CONTROL_SD_PALETTE <- c( "#969696ff", "#3caa3cff")
METADATA <- "metadata.csv"
CACHE <- "./cache/"
#~ RESULT_DIR <- "./raw_results/"
RESULT_DIR <- "/data/ethoscope_results"

REMOTE_DIR <- "ftp://nas.lab.gilest.ro/auto_generated_data/ethoscope_results/"

met <- fread(METADATA)
met <- met[status == "OK" & sdi != 0]
#~ met <- link_ethoscope_metadata_remote(met,
#~                                       remote_dir =  REMOTE_DIR,
#~                                       result_dir = RESULT_DIR,
#~                                       verbose = TRUE)

met <- link_ethoscope_metadata(met, result_dir = RESULT_DIR)
   
#met <- met[machine_name == "ETHOSCOPE_009"]
                                   
                                                                            

max_velocity <- function(sdt,  velocity_correction_coef =3e-3,
                                   masking_duration=0){
  d <- copy(sdt)
  d[,dt := c(NA,diff(t))]
  #d[,surface_change := xor_dist * 1e-3]
  d[,dist := 10^(xy_dist_log10x1000/1000) ]
  d[,velocity := dist/dt]

  a = velocity_correction_coef

  # masking here
  if(!"has_interacted" %in% colnames(d)){
    if(masking_duration >0)
      warning("Data does not contain an `has_interacted` column.
              Cannot apply masking!.
              Set `masking_duration = 0` to ignore masking")
    d[, has_interacted := 0]
  }

  d[,interaction_id := cumsum(has_interacted)]
  d[,
    masked := t < (t[1] + masking_duration),
    by=interaction_id
    ]
  d[ ,velocity := ifelse(masked & interaction_id != 0, 0, velocity)]
  d[,interaction_id := NULL]
  d[,masked := NULL]
  # end of masking

  d[, velocity_corrected :=  velocity  * dt  /a]
  d
}

process_interactions <- function(sdt){
	sdt <- max_velocity(sdt)
	sdt[, interaction_id := cumsum(has_interacted)]
	int_sdt <- sdt[has_interacted==T]
	setnames(int_sdt, "t", "interaction_t")
	sdt <- sdt[int_sdt[, .(interaction_id, interaction_t)], on="interaction_id"]
	sdt <- sdt[, t_rel := t- interaction_t]
	sdt <- sdt[ t_rel < 20]
	sdt
}

dt <- load_ethoscope(met,
					   max_time=days(20),
					   reference_hour=9.0,
					   cache = CACHE,
					   FUN=process_interactions,
					   ncores=1)
					   
dt <- stitch_on(dt, on="uid")   


# first we remove animals that do not have enought data points
dt[, t := t - days(xmv(baseline_days))]
dt[, interaction_t := t - t_rel]
dt[,  interaction_id := cumsum(t_rel == 0), by=id]

#dt[, interaction_t := interaction_t - days(xmv(baseline_days))]
dt[, treatment := as.factor(ifelse(sdi == 0, "Control", "SD")), meta=T]
dt[, interval := round(((11-sdi) ^ 1.7)) * 20, meta=T]	
dt[, interval := plyr::revalue(as.factor(interval),c("1180"="Control")), meta=T]





int_dt <- dt[, .(
				t=unique(interaction_t),
				max_velo = max(velocity_corrected[t_rel > 10 &t_rel <=20]),
				med_velo = median(velocity_corrected[t_rel > 10 &t_rel <=20])) , 
				
			  by="id,interaction_id"]			  
setkeyv(int_dt, "id")
int_dt <- behavr(int_dt, dt[meta=T])

int_dt
ggetho(int_dt[xmv(sdi) ==10], aes(z=max_velo>1, y=paste0(sex,"|",id))) + stat_tile_etho() 

pdf("prolonged_sd_stimulus_efficiency.pdf", w=12,h=6)
ggetho(int_dt, aes(y=max_velo>1,color=sex)) + 
			stat_ld_annotations(height=1, colour=NA, alpha=.3,ld_colours = c("white", "grey")) +
			stat_pop_etho() + 
			facet_grid( sex ~ .) +
			scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Moving response") 
			
ggetho(int_dt[t>days(1)], aes(y=max_velo>1,color=sex), time_wrap=hours(24)) + 
			stat_pop_etho() + 
			facet_grid( sex ~ .) +
			scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Moving response") + stat_ld_annotations()
			
ggetho(int_dt, aes(y=max_velo>2.5,color=sex)) + 
			stat_ld_annotations(height=1, colour=NA, alpha=.3,ld_colours = c("white", "grey")) +
			stat_pop_etho() + 
			facet_grid( sex ~ .) +
			scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Walking response") 


ggetho(int_dt[t>days(1)], aes(y=max_velo>2.5,color=sex), time_wrap=hours(24)) + 
			stat_pop_etho() + 
			facet_grid( sex ~ .) +
			scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Walking response") + 
			stat_ld_annotations()

ggetho(int_dt, aes(y=med_velo,color=sex), summary_FUN=median) + 
			stat_ld_annotations(height=1, colour=NA, alpha=.3,ld_colours = c("white", "grey")) +
			stat_pop_etho() + 
			facet_grid( sex ~ .) +
			scale_y_continuous(name="Median velocity post stim") + scale_y_sqrt() + geom_hline(yintercept=c(1,2.5), alpha=.5, col="red")
			
ggetho(int_dt[t>days(1)], aes(y=med_velo,color=sex), summary_FUN=median, time_wrap=hours(24)) + 
			stat_pop_etho() + 
			facet_grid( sex ~ .) +
			scale_y_continuous(name="Median velocity post stim") + scale_y_sqrt() + geom_hline(yintercept=c(1,2.5), alpha=.5, col="red")+ 
			stat_ld_annotations(colour="black")


dev.off()
