# Clone correction
Alejandro Rojas  






##Clone correction

```r
#Clone correction
mcc_TY <- clonecorrect(ultimhier, strata = ~County/Season, keep = c(1,2))
mcc_TY
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##    65 original multilocus genotypes 
##    79 diploid individuals
##     6 codominant loci
## 
## Population information:
## 
##     2 strata - County, Season
##     6 populations defined - 
## kalamazoo_fall-11, kent_fall-11, wayne_spring-13, kalamazoo_spring-13, wayne_fall-12, kalamazoo_fall-12
```

```r
#Locus table before and after clone correction
setPop(ultimhier) <- ~County/Season
cc <- locus_table(mcc_TY, info = FALSE)
mp <- locus_table(ultimhier, info = FALSE)

#Locus Differences after clone correction
locus_diff <- mp - cc
barplot(locus_diff[, "1-D"], ylab = "Change in Simpson's Index", xlab = "Locus",
        main = "Comparison of clone-corrected vs. uncorrected data")
```

![](Clone_correction_files/figure-html/clones-1.png)<!-- -->


Clone correction by population
------------------------------


```r
plot_simp_diff <- function(pop_name, clone_corrected, un_corrected){
  # Step 1: calculate diversity for clone-corrected data
  cc <- locus_table(clone_corrected, pop = pop_name, info = FALSE)
  # Step 2: calculate diversity for uncorrected data
  uc <- locus_table(un_corrected, pop = pop_name, info = FALSE)
  # Step 3: Take the difference
  res <- uc - cc
  # Step 4: Plot Simpson's index.
  barplot(res[, "1-D"], main = pop_name, ylab = "Change in Simpson's Index", xlab = "Locus")
}

par(mfrow = c(2, 3)) # Set up the graphics to have two rows and three columns

for (i in popNames(ultimhier)){
  plot_simp_diff(i, mcc_TY, ultimhier)
}
```

![](Clone_correction_files/figure-html/clone_correction_pop-1.png)<!-- -->

Clone correction by season
--------------------------


```r
#clone correction by season
mcc_TS <- clonecorrect(ultimhier, strata = ~Season, keep = 1:1)
mcc_TS
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##    65 original multilocus genotypes 
##    76 diploid individuals
##     6 codominant loci
## 
## Population information:
## 
##     2 strata - County, Season
##     3 populations defined - fall-11, spring-13, fall-12
```

```r
#Compare before and after correction
#First Setting uncorrected pop to season
setPop(ultimhier) <- ~Season
plot_simp_diff("all",mcc_TS,ultimhier)
```

![](Clone_correction_files/figure-html/clone_correc_season-1.png)<!-- -->

```r
#harmonic mean to test the effective population size 
locus_table(ultimhier, pop = "fall-11")
```

```
## 
## allele = Number of observed alleles
## 1-D = Simpson index
## Hexp = Nei's 1978 gene diversity
## ------------------------------------------
```

```
##       summary
## locus  allele  1-D Hexp Evenness
##   Py28   6.00 0.74 0.75     0.80
##   Py62   2.00 0.38 0.38     0.80
##   Py69   2.00 0.49 0.50     0.99
##   Py30   3.00 0.59 0.59     0.85
##   Py55   4.00 0.74 0.74     0.97
##   Py57   3.00 0.21 0.21     0.56
##   mean   3.33 0.53 0.53     0.83
```

```r
locus_table(ultimhier, pop = "fall-12")
```

```
## 
## allele = Number of observed alleles
## 1-D = Simpson index
## Hexp = Nei's 1978 gene diversity
## ------------------------------------------
```

```
##       summary
## locus  allele  1-D Hexp Evenness
##   Py28   6.00 0.63 0.64     0.70
##   Py62   3.00 0.33 0.34     0.60
##   Py69   2.00 0.49 0.49     0.97
##   Py30   3.00 0.48 0.49     0.71
##   Py55   5.00 0.65 0.66     0.71
##   Py57   2.00 0.45 0.46     0.92
##   mean   3.50 0.51 0.51     0.77
```

```r
locus_table(ultimhier, pop = "spring-13")
```

```
## 
## allele = Number of observed alleles
## 1-D = Simpson index
## Hexp = Nei's 1978 gene diversity
## ------------------------------------------
```

```
##       summary
## locus  allele  1-D Hexp Evenness
##   Py28   4.00 0.60 0.62     0.73
##   Py62   2.00 0.31 0.32     0.72
##   Py69   2.00 0.39 0.40     0.81
##   Py30   3.00 0.40 0.42     0.63
##   Py55   4.00 0.69 0.72     0.87
##   Py57   2.00 0.48 0.50     0.97
##   mean   2.83 0.48 0.50     0.79
```

