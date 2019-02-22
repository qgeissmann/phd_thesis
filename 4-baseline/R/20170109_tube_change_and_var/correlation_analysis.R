rm(list=ls());gc()
library(behavr)
library(scopr)
library(sleepr)
library(ggetho)


QUERY <- "query_complete.csv"
CACHE <- "~/Desktop/ethoscope_cache/"

q <- parse_query(QUERY, "/data/ethoscope_results/")
q <- q[status == "OK"]

dt <- query_ethoscopes(q,
					   reference_hour = 9.0,
#					   max_time=days(6),
					   columns=c("x", "has_interacted", "xy_dist_log10x1000"),
					   cache = CACHE,
					   FUN = sleep_annotation
					   )


valid_animals <- dt[,.(max_t = max(t)), by=id][max_t >= days(7)]$id


# this updates metadata
dt <- dt[id %in% valid_animals, verbose=T]

dt_stchd <- stitch_on(dt[xmv(uid) != 86], on="uid")

valid_animals <- dt_stchd[,.(max_t = max(t)), by=id][max_t >= days(14)]$id
dt_stchd <- dt_stchd[id %in% valid_animals, verbose=T]

pl <- ggetho(dt, aes(z=asleep)) + stat_tile_etho()
pl2 <- ggetho(dt_stchd, aes(z=asleep, y = id)) + stat_tile_etho()




summary_dt <- dt_stchd[, 
				.(sleep_days_02_07 =  mean(asleep[t>days(2) & t < days(7)])
				  ),
				  by = id]
setmeta(dt_stchd, rejoin(summary_dt))
sleep_days_02_06_str := sprintf("%04d",round(1000 * sleep_days_02_06)
dt_stchd[, 
    tag_sex_sleep_uid := sprintf("%s|%04d|uid",sex,round(1000 * sleep_days_02_07), uid),
    meta=T]

pl2 <- ggetho(dt_stchd, aes(z=asleep, y = tag_sex_sleep_uid)) + 
            stat_tile_etho()




dt_stchd[, LD := as.factor(ifelse(t %% hours(24) > hours(12), "D", "L"))]


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
#~ 				  sleep_days_03_04 = mean(!moving[t>days(3) & t < days(4)]),
#~ 				  sleep_days_04_05 = mean(!moving[t>days(3) & t < days(4)]),
#~ 				  sleep_days_05_06 = mean(!moving[t>days(5) & t < days(6)]),
#~ 				  sleep_days_07_08 = mean(!moving[t>days(7) & t < days(8)]),
#~ 				  sleep_days_09_10 = mean(!moving[t>days(9) & t < days(10)])
				  ),
				  by = c("id")]
				  
stat_dt <- rejoin(stat_dt)


pdf("tube_change_correlations.pdf", w=12,h=6)


cor.test(stat_dt[sex=="F"]$sleep_days_02_04,
		 stat_dt[sex=="F"]$sleep_days_05_07, method="spearman")

cor.test(stat_dt[sex=="F"]$sleep_days_05_07,
		 stat_dt[sex=="F"]$sleep_days_08_10, method="spearman")

cor.test(stat_dt[sex=="F"]$sleep_days_08_10,
		 stat_dt[sex=="F"]$sleep_days_11_13, method="spearman")
		 


cor.test(stat_dt[sex=="F"]$sleep_days_02_07,
		 stat_dt[sex=="F"]$sleep_days_08_13, method="spearman")

cor.test(stat_dt[sex=="M"]$sleep_days_02_07,
		 stat_dt[sex=="M"]$sleep_days_08_13, method="spearman")

		 
		 
ggplot(stat_dt, aes(sleep_days_02_07, sleep_days_08_13, colour = sex, shape=sex, fill=sex)) +
			geom_smooth(method="lm")+
			geom_point(size=2, alpha=.5) + labs(x="Sleep before tube change", y="Sleep after tube change") 
			#+
#~ 			generic_theme
#~ 			scale_fill_manual(values=FEMALE_MALE_PALETTE) +
#~ 			scale_colour_manual(values=FEMALE_MALE_PALETTE)			
			
dev.off()

