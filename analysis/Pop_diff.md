---
title: "Population differantiation"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---







# Assesing population structure using Hendrick's Gst


```r
#Gst
Gst_Hedrick(clonecorrect(ultimhier, ~County/Season, keep = 1:2))
```

```
## $per.locus
##         Py28         Py62         Py69         Py30         Py55 
##  0.514023333  0.064589738 -0.006064058  0.126695184  0.629411956 
##         Py57 
##  0.481884660 
## 
## $global
## [1] 0.3201126
```

# Genetic distance


## County and season using Nei's distance

```r
set.seed(999)
ultimhier %>%
genind2genpop(pop = ~County/Season) %>%
aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = nei.dist)
```

```
## 
##  Converting data from a genind to a genpop object... 
## 
## ...done.
```

![](Pop_diff_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```
## 
## Phylogenetic tree with 6 tips and 5 internal nodes.
## 
## Tip labels:
## [1] "kalamazoo_fall-11"   "kent_fall-11"        "wayne_spring-13"    
## [4] "kalamazoo_spring-13" "wayne_fall-12"       "kalamazoo_fall-12"  
## Node labels:
## [1] 100.0    NA    NA  78.8  59.5
## 
## Rooted; includes branch lengths.
```

## Tree using Bruvo's distance


```r
#Repeat length
ssr.reps <- c(3, 3, 6, 6, 2, 3)

#Bruvo
setPop(ult.cc) <- ~Season
ult.tree <- bruvo.boot(ult.cc, 
                       replen = ssr.reps, 
                       sample = 1000,
                       tree ="nj", 
                       cutoff = 50, 
                       quiet = TRUE)
```

```
## Warning in bruvo.boot(ult.cc, replen = ssr.reps, sample = 1000, tree =
## "nj", : Some branch lengths of the tree are negative. Normalizing branches
## according to Kuhner and Felsenstein (1994)
```

![](Pop_diff_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
ult.tree$tip.label <- stringr::str_trim(ult.tree$tip.label) %>% 
  stringr::str_replace_all(" ","")


Py_metadata <- read.csv("../data/Py_metadata.csv") %>%
  select(Isolate, Origin, Season, Year, Species) %>%
  column_to_rownames("Isolate")

ult.tr2 <- ggtree(ult.tree, branch.length = "none") + 
  geom_tiplab() + 
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50), hjust =-.3) + 
  xlim_tree(10) +
  theme_tree()
```



```r
colors <- c("#cf4661","#bf7658","#ce5d2a",
            "#c4a643","#6e7a32","#6ab544",
            "#53aa78","#4cb8d1","#7385cb",
            "#7e63cb","#c851b8","#bf679b")

#Adding heatmap
ult.tr3 <- gheatmap(ult.tr2, Py_metadata, colnames = FALSE, width = 0.6, offset = 1)
```

```
## Warning: attributes are not identical across measure variables;
## they will be dropped
```

```r
#Formatting colors and legend for heatmap
lbl <- get_heatmap_column_position(ult.tr3, by="top")
ult.tr3 + scale_fill_manual(breaks=c("Kalamazoo","Kent","Wayne","fall-11","fall-12","spring-13","2011","2012","2013",
                                       "P. ultimum ", "P. ultimum var. ultimum"), values = colors) +
  geom_text(data=lbl, aes(x, y, label=label), 
            nudge_y = 1.5, 
            nudge_x = 0.5, 
            angle=45,
            size=4)
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill',
## which will replace the existing scale.
```

```
## Warning in fun(x, ...): NAs introduced by coercion
```

```
## Warning in fun(x, ...): NAs introduced by coercion
```

```
## Warning in FUN(X[[i]], ...): NAs introduced by coercion

## Warning in FUN(X[[i]], ...): NAs introduced by coercion
```

<img src="Pop_diff_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />



# AMOVA


```r
#AMOVA
data(ultimhier)
ultimhier
splithierarchy(ultimhier) <- ~County/Season
ultimhier
setpop(ultimhier)<- ~Season 
ultimhier
#Replacing mising data with 0
nanzero <- missingno(ultimhier, "zero")
nanind <- missingno(ultimhier, "geno")
table(gethierarchy(ultimhier, ~County))
table(gethierarchy(nanind, ~County/Season, combine = FALSE))
ultimhieramova <- poppr.amova(ultimhier, ~County/Season)
ultimhieramova0 <- poppr.amova(nanind, ~County/Season)
ultimzhienamova <-poppr.amova(ultimhier, ~County/Season, within = FALSE)
ultimamovacc <- poppr.amova(ultimhier, ~County/Season, clonecorrect = TRUE)
ultimamovacc0 <- poppr.amova(nanind, ~County/Season, clonecorrect = TRUE)
ultimhieramova
ultimhieramova0
ultimzhienamova
ultimamovacc
ultimamovacc0
write.table(ultimhieramova$componentsofcovariance, sep = ",", file = "~/Documents/ultimhierAMOVA.csv")
write.table(ultimhier$statphi, sep = ",", file = "~/Documents/ultimhierphiAMOVA.csv")
#significance
ultimhier
set.seed(1999)
ultmsignif   <- randtest(ultimhieramova0, nrepet = 999)
Aeutccsignif <- randtest(ultimamovacc0, nrepet = 1000)
plot(ultmsignif)
ultmsignif
pdf(file = "~/Documents/significance AMOVA.pdf", width = 5, height = 5)
plot (Aeutccsignif)
dev.off()
Aeutccsignif
write.table((ultmsignif), sep = ",", file = "~/Documents/significancewithoutAMOVA.csv")
plot(Aeutccsignif)
```


