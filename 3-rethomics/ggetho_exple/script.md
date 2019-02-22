---
title: "Untitled"
output: pdf_document
---









```r
dt[ ,treatment := ifelse(sdi == 0, "Control", "Treated"), meta=T]
pl <- 
ggetho(dt, aes(x = t,
               y = moving,
               colour = sex)) +
    stat_pop_etho(method = ggplot2::mean_cl_boot) +
    stat_ld_annotations() +
    scale_x_days() +
    facet_grid( treatment ~ .) 
```

```
## Scale for 'x' is already present. Adding another scale for 'x', which
## will replace the existing scale.
```


```r
pdf("plot.pdf", w=8, h=4)
pl
dev.off()
```

```
## png 
##   2
```
