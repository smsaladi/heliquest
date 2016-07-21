

# Helical wheel plots in R

## Installation


```r
# Until this package is on CRAN:
devtools::install_github("smsaladi/heliquest")
```

# Examples


```r
library(heliquest)

draw_helical_wheel("KPKPKPKPKPKPKPKPKPKPKPKPKPKPK")
```

```
## Loading required package: httr
```

```
## Warning: package 'httr' was built under R version 3.2.5
```

<img src="README_files/figure-html/examples-1.png" width="576" />

```r
draw_helical_wheel("ACDEFGHIKLMNPQRSTVWY",
                   FactC = 0.05, FONT1 = 3, FONT2 = 5,
                   CEXFT = 1, CEXTEXT = 0.8,
                   FlFH = 0, ANGT = 100,
                   NBMIN = 18, NBM2 = 36, NBMAX = 54)
```

<img src="README_files/figure-html/examples-2.png" width="576" />

```r
# Custom Angle and Hydrophobicity
draw_helical_wheel("ACDEFGHIKLMNPQRSTVWY",
                   FactC = 0.05, FONT1 = 3, FONT2 = 5,
                   CEXFT = 1, CEXTEXT = 0.8,
                   Ang = 30, Mom = 1,
                   FlFH = 0, ANGT = 100,
                   NBMIN = 18, NBM2 = 36, NBMAX = 54)
```

<img src="README_files/figure-html/examples-3.png" width="576" />
