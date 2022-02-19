---
author : "L. Cauquil"
title: "Modèle linéaire mixte avec les packages purrr et broom.mixed"
date: "26-01-2022"
output:
  html_document: 
    code_folding: show
    toc: yes
    toc_float: yes
    keep_md: TRUE
  pdf_document:
    toc: yes
  word_document:
    toc: yes
editor_options: 
  chunk_output_type: inline
---



## **Packages et importation des données**

### **Packages**


```r
## Manipulations des données
library(tibble)
library(dplyr)
library(tidyr)

## Mise en forme des résultats stat
library(purrr)
library(broom)
library(broom.mixed)

## Objet phyloseq
library(phyloseq)

## Fonctions stat
library(car)
library(multcomp)
library(lmerTest)
library(emmeans)

## Présentation des tables
library(DT)

## Visualisations
library(ggplot2)
```

### **Importation des données**

On utilise l'objet phyloseq `tab_phylum`


```r
load("data/data_phylum.RData")
tab_Phylum
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 9 taxa and 40 samples ]
## sample_data() Sample Data:       [ 40 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 9 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 9 tips and 8 internal nodes ]
```

Il y a 9 phyla présents dans l'écosystème.


## **Construction des objets**

Différents objets créés:

 - `RecapPhylum`: data.frame avec les moyennes, sd et SEM des abondances relatives de chaque phylum sans tenir compte des groupes
 - `top_Phylum`: vecteur de Phyla qui ont au moins une moyenne des abondances relatives au sein d'un groupe Segment_Age > à 0,005 
 - Rajout de la variable `Firm_bact_ratio` = ration Firmicutes/Bacteroidota
 - `df_stat`: data.frame utilisé pour appliqué un modèle linéaire mixte sur les abondances relatives des phyla sélectionnés + `Firm_bact_ratio`

Les analyses sont faites sur les abondances relatives


```r
## Transforme en abondance relative
relatabunphy <- transform_sample_counts(tab_Phylum, function(OTU) OTU/sum(OTU))
```

La fonction `psmelt()` fusionne les tables `otu_table`, `sample_otu` et `tax_table` d'un objet phyloseq et crée un seul data.frame au format long.  
Les phyla sont regroupés dans une seule colonne Phylum !


```r
dfrap <- psmelt(relatabunphy) 
head(dfrap)
```

```
##          OTU                      Sample Abundance       GPS Piglet Age Sexe
## 39 Cluster_1 GPS160115_ATGAAC-JLGMN_L001 0.9995071 GPS160115     14 D35    F
## 38 Cluster_1 GPS161671_CTCTAC-JLGMN_L001 0.9993956 GPS161671     23 D35    M
## 19 Cluster_1 GPS158690_CCTTGA-JLGMN_L001 0.9984081 GPS158690      5 D35    M
## 22 Cluster_1 GPS158965_CACCCA-JLGMN_L001 0.9983199 GPS158965      9 D21    F
## 18 Cluster_1 GPS161675_ACCGTG-JLGMN_L001 0.9981725 GPS161675     24 D35    M
## 26 Cluster_1 GPS160464_CCCAAA-JLGMN_L001 0.9978882 GPS160464     20 D21    F
##    Segment Segment_Age Concentration.ADN                    SampleID
## 39 Jejunum Jejunum_D35             7,736 GPS160115_ATGAAC-JLGMN_L001
## 38 Jejunum Jejunum_D35              9,75 GPS161671_CTCTAC-JLGMN_L001
## 19 Jejunum Jejunum_D35             75,69 GPS158690_CCTTGA-JLGMN_L001
## 22 Jejunum Jejunum_D21             25,34 GPS158965_CACCCA-JLGMN_L001
## 18 Jejunum Jejunum_D35             5,389 GPS161675_ACCGTG-JLGMN_L001
## 26 Jejunum Jejunum_D21             18,19 GPS160464_CCCAAA-JLGMN_L001
##    Slaughter_date  Kingdom     Phylum
## 39             D3 Bacteria Firmicutes
## 38             D5 Bacteria Firmicutes
## 19             D1 Bacteria Firmicutes
## 22             D2 Bacteria Firmicutes
## 18             D5 Bacteria Firmicutes
## 26             D4 Bacteria Firmicutes
```

### **Tableau général avec mean, sd, SEM des abondances relatives par phylum**

Construction de la table: tableau général sans tenir des groupes  


```r
RecapPhylum <- dfrap %>%  
  dplyr::select(Phylum,Abundance) %>% 
  group_by(Phylum) %>% 
  summarise(data.frame(mean = mean(Abundance),
                       sd = sd(Abundance), 
                       sem = sd(Abundance) / sqrt(length(Abundance))))
RecapPhylum[,2:4] <- round(RecapPhylum[,2:4],4)*100

datatable(RecapPhylum)
```

```{=html}
<div id="htmlwidget-b4d8925adc6cf284c63e" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-b4d8925adc6cf284c63e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9"],["Actinobacteriota","Bacteroidota","Campylobacterota","Desulfobacterota","Firmicutes","Fusobacteriota","Patescibacteria","Proteobacteria","Spirochaetota"],[0.43,11.58,0.11,0.55,82.87,0.81,0.06,3.56,0.03],[0.79,14.56,0.27,0.88,18.7,2.55,0.3,8.27,0.09],[0.12,2.3,0.04,0.14,2.96,0.4,0.05,1.31,0.01]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Phylum<\/th>\n      <th>mean<\/th>\n      <th>sd<\/th>\n      <th>sem<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


### **Sélection des phyla avec au moins une moyenne d'abondance relative > à 0,005 au sein d'un groupe Segment_Age**
  
Sélection des phyla  


```r
top_Phylum <- dfrap %>%  
  dplyr::select(Phylum, Abundance, Segment_Age) %>% 
  group_by(Phylum,Segment_Age) %>% 
  summarise_all(list(mean = mean)) %>% 
  dplyr::arrange(desc(mean)) %>% 
  dplyr::filter(mean >= 0.005) %>% 
  dplyr::arrange(Phylum, Segment_Age) %>% 
  distinct(Phylum) %>% 
  as_vector()
top_Phylum
```

```
##          Phylum1          Phylum2          Phylum3          Phylum4 
## Actinobacteriota     Bacteroidota Desulfobacterota       Firmicutes 
##          Phylum5          Phylum6 
##   Fusobacteriota   Proteobacteria 
## 9 Levels: Actinobacteriota Bacteroidota Campylobacterota ... Spirochaetota
```

Selection des phyla dans la table générale


```r
dfrap <- dfrap %>% 
  dplyr::filter(Phylum %in% top_Phylum)
```

### **Rajout de la variable `Firm_bact_ratio` aux phyla**


```r
dfrap %>% 
  dplyr::select(Phylum, Abundance, Segment_Age, Piglet, Slaughter_date, Segment, Age) %>% 
#  arrange(Phylum) %>% 
  pivot_wider(names_from = Phylum,
              values_from = Abundance) %>% 
  mutate(Firm_bact_ratio = Firmicutes/Bacteroidota) %>% 
  pivot_longer(cols = where(is.numeric),
               names_to = "phylum", 
               values_to = "abundance") -> df_stat
```

Transformation des abundances en racine de racine


```r
df_stat$abundance <- df_stat$abundance^0.25
```

## **Statistiques**  

**Paramètres des contrastes**


```r
options(contrasts = c("contr.sum", "contr.poly"))
```

**Modèle linéaire mixte**

Le modèle utilise 

- la fonction `lmer` du package `lmerTest` 
- l'effet piglet en variable aléatoire (1|Piglet) 

Mais Piglet est niché dans slaughter_date

- l'addition de Slaughter_date sur (1|Slaughter_date/Piglet)

NB l'écriture (1|Slaughter_date/Piglet) est équivalente à (1|Piglet) +(1|Slaughter_date) car les porcs ont des identifiants différents entre les date d'abattage
plus d'explication ici https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified  

https://github.com/lcauquil/tuto_mix_model/blob/master/Modele_mixte_guidelines_2019.pdf pour les explications REML (REML pour restricted maximum likelihood) et ML (maximum likelihood)

### **Modèle utilisé**

**lmer(abundance ~ Segment * Age + (1|Piglet)+ (1|Slaughter_date), data = as.data.frame(data))**

Application du modèle sur le phylum "Bacteroidota"


```r
df_stat %>% 
  filter(phylum == "Bacteroidota") -> lm_bacte

lmer(abundance ~ Segment * Age + (1|Piglet)+ (1|Slaughter_date), data = lm_bacte)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: abundance ~ Segment * Age + (1 | Piglet) + (1 | Slaughter_date)
##    Data: lm_bacte
## REML criterion at convergence: -31.5725
## Random effects:
##  Groups         Name        Std.Dev.
##  Piglet         (Intercept) 0.14663 
##  Slaughter_date (Intercept) 0.00000 
##  Residual                   0.06558 
## Number of obs: 40, groups:  Piglet, 24; Slaughter_date, 6
## Fixed Effects:
##   (Intercept)       Segment1           Age1  Segment1:Age1  
##      0.463270      -0.266397       0.026846      -0.004721  
## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings
```

### **Fonctions**

On crées plusieurs fonctions à appliquer pour chacun des phylum.  
Fonctions pour :

 - appliquer le modèle
 - effectuer le test shapiro
 - calcluer les p-values
 - récupérer les lettres

Pour appliquer chaque fonction on utilise la fonction `map()` du package purrr. Les fonctions `map()` et dérivées visent à remplacer les boucles et les fonctions de la famille des apply.  

Seule contrainte des fonctions `map()`, elles s'appliquent uniquement à des listes. Il faut donc convertir les data.frames en liste, mais c'est facile avec la fonction `split()`.

**Exemple d'utilisation de la fonction `map()`**

La fonction `map()` à généralement 2 arguments: une liste et une fonction


```r
x <- list(a = 1:10, b = 5:110)
x
```

```
## $a
##  [1]  1  2  3  4  5  6  7  8  9 10
## 
## $b
##   [1]   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22
##  [19]  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40
##  [37]  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58
##  [55]  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76
##  [73]  77  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94
##  [91]  95  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110
```

```r
map(x, mean)
```

```
## $a
## [1] 5.5
## 
## $b
## [1] 57.5
```

```r
## sortie sous forme de data.frame
map_df(x, mean)
```

```
## # A tibble: 1 × 2
##       a     b
##   <dbl> <dbl>
## 1   5.5  57.5
```

```r
rm(x)
```

**Fonction split()**

La fonction `split()` sépare un data.frame en liste à partir des différents niveaux d'une variable
 

```r
split(df_stat, df_stat$phylum)
```

```
## $Actinobacteriota
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum           abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>                <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Actinobacteriota     0.115
##  2 Jejunum_D35 23     D5             Jejunum D35   Actinobacteriota     0.109
##  3 Jejunum_D35 5      D1             Jejunum D35   Actinobacteriota     0.142
##  4 Jejunum_D21 9      D2             Jejunum D21   Actinobacteriota     0.136
##  5 Jejunum_D35 24     D5             Jejunum D35   Actinobacteriota     0.145
##  6 Jejunum_D21 20     D4             Jejunum D21   Actinobacteriota     0.120
##  7 Jejunum_D21 27     D6             Jejunum D21   Actinobacteriota     0.203
##  8 Jejunum_D21 26     D6             Jejunum D21   Actinobacteriota     0.178
##  9 Jejunum_D21 18     D4             Jejunum D21   Actinobacteriota     0.149
## 10 Jejunum_D35 7      D1             Jejunum D35   Actinobacteriota     0.125
## # … with 30 more rows
## 
## $Bacteroidota
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum       abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>            <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Bacteroidota    0.0966
##  2 Jejunum_D35 23     D5             Jejunum D35   Bacteroidota    0.120 
##  3 Jejunum_D35 5      D1             Jejunum D35   Bacteroidota    0.129 
##  4 Jejunum_D21 9      D2             Jejunum D21   Bacteroidota    0.139 
##  5 Jejunum_D35 24     D5             Jejunum D35   Bacteroidota    0.150 
##  6 Jejunum_D21 20     D4             Jejunum D21   Bacteroidota    0.131 
##  7 Jejunum_D21 27     D6             Jejunum D21   Bacteroidota    0.0855
##  8 Jejunum_D21 26     D6             Jejunum D21   Bacteroidota    0.141 
##  9 Jejunum_D21 18     D4             Jejunum D21   Bacteroidota    0.120 
## 10 Jejunum_D35 7      D1             Jejunum D35   Bacteroidota    0.121 
## # … with 30 more rows
## 
## $Desulfobacterota
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum           abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>                <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Desulfobacterota         0
##  2 Jejunum_D35 23     D5             Jejunum D35   Desulfobacterota         0
##  3 Jejunum_D35 5      D1             Jejunum D35   Desulfobacterota         0
##  4 Jejunum_D21 9      D2             Jejunum D21   Desulfobacterota         0
##  5 Jejunum_D35 24     D5             Jejunum D35   Desulfobacterota         0
##  6 Jejunum_D21 20     D4             Jejunum D21   Desulfobacterota         0
##  7 Jejunum_D21 27     D6             Jejunum D21   Desulfobacterota         0
##  8 Jejunum_D21 26     D6             Jejunum D21   Desulfobacterota         0
##  9 Jejunum_D21 18     D4             Jejunum D21   Desulfobacterota         0
## 10 Jejunum_D35 7      D1             Jejunum D35   Desulfobacterota         0
## # … with 30 more rows
## 
## $Firm_bact_ratio
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum          abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>               <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Firm_bact_ratio     10.4 
##  2 Jejunum_D35 23     D5             Jejunum D35   Firm_bact_ratio      8.31
##  3 Jejunum_D35 5      D1             Jejunum D35   Firm_bact_ratio      7.73
##  4 Jejunum_D21 9      D2             Jejunum D21   Firm_bact_ratio      7.20
##  5 Jejunum_D35 24     D5             Jejunum D35   Firm_bact_ratio      6.65
##  6 Jejunum_D21 20     D4             Jejunum D21   Firm_bact_ratio      7.63
##  7 Jejunum_D21 27     D6             Jejunum D21   Firm_bact_ratio     11.7 
##  8 Jejunum_D21 26     D6             Jejunum D21   Firm_bact_ratio      7.10
##  9 Jejunum_D21 18     D4             Jejunum D21   Firm_bact_ratio      8.35
## 10 Jejunum_D35 7      D1             Jejunum D35   Firm_bact_ratio      8.29
## # … with 30 more rows
## 
## $Firmicutes
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum     abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>          <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Firmicutes     1.00 
##  2 Jejunum_D35 23     D5             Jejunum D35   Firmicutes     1.00 
##  3 Jejunum_D35 5      D1             Jejunum D35   Firmicutes     1.00 
##  4 Jejunum_D21 9      D2             Jejunum D21   Firmicutes     1.00 
##  5 Jejunum_D35 24     D5             Jejunum D35   Firmicutes     1.00 
##  6 Jejunum_D21 20     D4             Jejunum D21   Firmicutes     0.999
##  7 Jejunum_D21 27     D6             Jejunum D21   Firmicutes     0.999
##  8 Jejunum_D21 26     D6             Jejunum D21   Firmicutes     0.999
##  9 Jejunum_D21 18     D4             Jejunum D21   Firmicutes     0.999
## 10 Jejunum_D35 7      D1             Jejunum D35   Firmicutes     0.999
## # … with 30 more rows
## 
## $Fusobacteriota
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum         abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>              <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Fusobacteriota    0     
##  2 Jejunum_D35 23     D5             Jejunum D35   Fusobacteriota    0     
##  3 Jejunum_D35 5      D1             Jejunum D35   Fusobacteriota    0     
##  4 Jejunum_D21 9      D2             Jejunum D21   Fusobacteriota    0     
##  5 Jejunum_D35 24     D5             Jejunum D35   Fusobacteriota    0     
##  6 Jejunum_D21 20     D4             Jejunum D21   Fusobacteriota    0.120 
##  7 Jejunum_D21 27     D6             Jejunum D21   Fusobacteriota    0.0855
##  8 Jejunum_D21 26     D6             Jejunum D21   Fusobacteriota    0     
##  9 Jejunum_D21 18     D4             Jejunum D21   Fusobacteriota    0.0846
## 10 Jejunum_D35 7      D1             Jejunum D35   Fusobacteriota    0.0770
## # … with 30 more rows
## 
## $Proteobacteria
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum         abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>              <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Proteobacteria     0.123
##  2 Jejunum_D35 23     D5             Jejunum D35   Proteobacteria     0.126
##  3 Jejunum_D35 5      D1             Jejunum D35   Proteobacteria     0.172
##  4 Jejunum_D21 9      D2             Jejunum D21   Proteobacteria     0.176
##  5 Jejunum_D35 24     D5             Jejunum D35   Proteobacteria     0.172
##  6 Jejunum_D21 20     D4             Jejunum D21   Proteobacteria     0.194
##  7 Jejunum_D21 27     D6             Jejunum D21   Proteobacteria     0.141
##  8 Jejunum_D21 26     D6             Jejunum D21   Proteobacteria     0.167
##  9 Jejunum_D21 18     D4             Jejunum D21   Proteobacteria     0.199
## 10 Jejunum_D35 7      D1             Jejunum D35   Proteobacteria     0.218
## # … with 30 more rows
```

**Fonction du modèle**


```r
fit_model <- function(data) {
  mod_REML <- lmer(abundance ~ Segment * Age + (1|Piglet)+ (1|Slaughter_date), data = as.data.frame(data))
	mod_REML_update <- update(mod_REML, REML = F)
}
```

**Application de la fonction `fit_model` sur les Firmicutes**  


```r
df_stat %>% 
  dplyr::filter(phylum == "Firmicutes") %>% 
  fit_model() %>% 
  summary()
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's
##   method [lmerModLmerTest]
## Formula: abundance ~ Segment * Age + (1 | Piglet) + (1 | Slaughter_date)
##    Data: as.data.frame(data)
## 
##      AIC      BIC   logLik deviance df.resid 
##   -108.8    -97.0     61.4   -122.8       33 
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -1.79409 -0.24396  0.06952  0.20530  1.87531 
## 
## Random effects:
##  Groups         Name        Variance  Std.Dev.
##  Piglet         (Intercept) 0.0045381 0.06737 
##  Slaughter_date (Intercept) 0.0000000 0.00000 
##  Residual                   0.0005582 0.02363 
## Number of obs: 40, groups:  Piglet, 24; Slaughter_date, 6
## 
## Fixed effects:
##                Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)    0.935449   0.014454 23.079580  64.720  < 2e-16 ***
## Segment1       0.034004   0.004453 14.911031   7.637 1.58e-06 ***
## Age1          -0.016489   0.014454 23.079580  -1.141    0.266    
## Segment1:Age1  0.004979   0.004453 14.911031   1.118    0.281    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) Sgmnt1 Age1  
## Segment1    -0.127              
## Age1        -0.034  0.112       
## Segmnt1:Ag1  0.112 -0.363 -0.127
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

**Application de la fonction `fit_model` sur tous les phyla**  


```r
df_stat %>% 
  split(df_stat$phylum) %>% 
  map(fit_model) -> result
```

**Resultat pour les Firmicutes**


```r
result$Firmicutes %>% summary()
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's
##   method [lmerModLmerTest]
## Formula: abundance ~ Segment * Age + (1 | Piglet) + (1 | Slaughter_date)
##    Data: as.data.frame(data)
## 
##      AIC      BIC   logLik deviance df.resid 
##   -108.8    -97.0     61.4   -122.8       33 
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -1.79409 -0.24396  0.06952  0.20530  1.87531 
## 
## Random effects:
##  Groups         Name        Variance  Std.Dev.
##  Piglet         (Intercept) 0.0045381 0.06737 
##  Slaughter_date (Intercept) 0.0000000 0.00000 
##  Residual                   0.0005582 0.02363 
## Number of obs: 40, groups:  Piglet, 24; Slaughter_date, 6
## 
## Fixed effects:
##                Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)    0.935449   0.014454 23.079580  64.720  < 2e-16 ***
## Segment1       0.034004   0.004453 14.911031   7.637 1.58e-06 ***
## Age1          -0.016489   0.014454 23.079580  -1.141    0.266    
## Segment1:Age1  0.004979   0.004453 14.911031   1.118    0.281    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) Sgmnt1 Age1  
## Segment1    -0.127              
## Age1        -0.034  0.112       
## Segmnt1:Ag1  0.112 -0.363 -0.127
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

### **Mise en forme des résultats**

Pour mettre en forme les résultats on peut utiliser la fonction `tidy()` des packages broom et broom.mixed qui récupèrent les résultats des tests et les formatent en data.frame, toujours en utilisant la fonction `map()`

La fonction `tidy()` s'applique à de nombreux tests:

[Listes des méthodes prises en charge par broom](https://broom.tidymodels.org/articles/available-methods.html){target="blank"}

[Listes des méthodes prises en charge par broom.mixed](https://cran.r-project.org/web/packages/broom.mixed/vignettes/broom_mixed_intro.html){target="blank"}

Mise en forme des résultats pour les Firmicutes


```r
tidy(result$Firmicutes)
```

```
## # A tibble: 7 × 8
##   effect   group          term      estimate std.error statistic    df   p.value
##   <chr>    <chr>          <chr>        <dbl>     <dbl>     <dbl> <dbl>     <dbl>
## 1 fixed    <NA>           (Interce…  0.935     0.0145      64.7   23.1  1.33e-27
## 2 fixed    <NA>           Segment1   0.0340    0.00445      7.64  14.9  1.58e- 6
## 3 fixed    <NA>           Age1      -0.0165    0.0145      -1.14  23.1  2.66e- 1
## 4 fixed    <NA>           Segment1…  0.00498   0.00445      1.12  14.9  2.81e- 1
## 5 ran_pars Piglet         sd__(Int…  0.0674   NA           NA     NA   NA       
## 6 ran_pars Slaughter_date sd__(Int…  0        NA           NA     NA   NA       
## 7 ran_pars Residual       sd__Obse…  0.0236   NA           NA     NA   NA
```

Appliqué à tous les phyla


```r
result_model <- df_stat %>%
  split(df_stat$phylum) %>% 
  map(fit_model) %>% 
  map(broom.mixed::tidy, .id = "phylum")
result_model
```

```
## $Actinobacteriota
## # A tibble: 7 × 8
##   effect   group          term      estimate std.error statistic    df   p.value
##   <chr>    <chr>          <chr>        <dbl>     <dbl>     <dbl> <dbl>     <dbl>
## 1 fixed    <NA>           (Interce…  0.220     0.0145     15.1    21.7  5.34e-13
## 2 fixed    <NA>           Segment1  -0.0124    0.00974    -1.27   15.3  2.23e- 1
## 3 fixed    <NA>           Age1       0.0130    0.0145      0.896  21.7  3.80e- 1
## 4 fixed    <NA>           Segment1…  0.00247   0.00974     0.254  15.3  8.03e- 1
## 5 ran_pars Piglet         sd__(Int…  0.0527   NA          NA      NA   NA       
## 6 ran_pars Slaughter_date sd__(Int…  0        NA          NA      NA   NA       
## 7 ran_pars Residual       sd__Obse…  0.0541   NA          NA      NA   NA       
## 
## $Bacteroidota
## # A tibble: 7 × 8
##   effect   group          term      estimate std.error statistic    df   p.value
##   <chr>    <chr>          <chr>        <dbl>     <dbl>     <dbl> <dbl>     <dbl>
## 1 fixed    <NA>           (Interce…  0.463      0.0310    15.0    22.6  3.37e-13
## 2 fixed    <NA>           Segment1  -0.267      0.0115   -23.2    14.4  7.60e-13
## 3 fixed    <NA>           Age1       0.0269     0.0310     0.867  22.6  3.95e- 1
## 4 fixed    <NA>           Segment1… -0.00476    0.0115    -0.414  14.4  6.85e- 1
## 5 ran_pars Piglet         sd__(Int…  0.141     NA         NA      NA   NA       
## 6 ran_pars Slaughter_date sd__(Int…  0         NA         NA      NA   NA       
## 7 ran_pars Residual       sd__Obse…  0.0612    NA         NA      NA   NA       
## 
## $Desulfobacterota
## # A tibble: 7 × 8
##   effect   group          term       estimate std.error statistic    df  p.value
##   <chr>    <chr>          <chr>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>
## 1 fixed    <NA>           (Intercep…  0.180      0.0114    15.8    11.8  2.55e-9
## 2 fixed    <NA>           Segment1   -0.141      0.0109   -12.9    11.1  4.64e-8
## 3 fixed    <NA>           Age1        0.00305    0.0114     0.269  11.8  7.93e-1
## 4 fixed    <NA>           Segment1:… -0.00312    0.0109    -0.286  11.1  7.80e-1
## 5 ran_pars Piglet         sd__(Inte…  0.0157    NA         NA      NA   NA      
## 6 ran_pars Slaughter_date sd__(Inte…  0         NA         NA      NA   NA      
## 7 ran_pars Residual       sd__Obser…  0.0640    NA         NA      NA   NA      
## 
## $Firm_bact_ratio
## # A tibble: 7 × 8
##   effect   group          term       estimate std.error statistic    df  p.value
##   <chr>    <chr>          <chr>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>
## 1 fixed    <NA>           (Intercep…    4.00      0.390    10.2   11.6   3.71e-7
## 2 fixed    <NA>           Segment1      2.79      0.289     9.64   7.75  1.39e-5
## 3 fixed    <NA>           Age1         -0.226     0.390    -0.579 11.6   5.74e-1
## 4 fixed    <NA>           Segment1:…   -0.106     0.289    -0.365  7.75  7.25e-1
## 5 ran_pars Piglet         sd__(Inte…    1.28     NA        NA     NA    NA      
## 6 ran_pars Slaughter_date sd__(Inte…    0        NA        NA     NA    NA      
## 7 ran_pars Residual       sd__Obser…    1.63     NA        NA     NA    NA      
## 
## $Firmicutes
## # A tibble: 7 × 8
##   effect   group          term      estimate std.error statistic    df   p.value
##   <chr>    <chr>          <chr>        <dbl>     <dbl>     <dbl> <dbl>     <dbl>
## 1 fixed    <NA>           (Interce…  0.935     0.0145      64.7   23.1  1.33e-27
## 2 fixed    <NA>           Segment1   0.0340    0.00445      7.64  14.9  1.58e- 6
## 3 fixed    <NA>           Age1      -0.0165    0.0145      -1.14  23.1  2.66e- 1
## 4 fixed    <NA>           Segment1…  0.00498   0.00445      1.12  14.9  2.81e- 1
## 5 ran_pars Piglet         sd__(Int…  0.0674   NA           NA     NA   NA       
## 6 ran_pars Slaughter_date sd__(Int…  0        NA           NA     NA   NA       
## 7 ran_pars Residual       sd__Obse…  0.0236   NA           NA     NA   NA       
## 
## $Fusobacteriota
## # A tibble: 7 × 8
##   effect   group          term       estimate std.error statistic    df  p.value
##   <chr>    <chr>          <chr>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>
## 1 fixed    <NA>           (Intercep…   0.132     0.0272      4.85  22.5  7.18e-5
## 2 fixed    <NA>           Segment1    -0.0262    0.0157     -1.67  15.1  1.15e-1
## 3 fixed    <NA>           Age1         0.0854    0.0272      3.14  22.5  4.68e-3
## 4 fixed    <NA>           Segment1:…  -0.0339    0.0157     -2.17  15.1  4.65e-2
## 5 ran_pars Piglet         sd__(Inte…   0.109    NA          NA     NA   NA      
## 6 ran_pars Slaughter_date sd__(Inte…   0        NA          NA     NA   NA      
## 7 ran_pars Residual       sd__Obser…   0.0857   NA          NA     NA   NA      
## 
## $Proteobacteria
## # A tibble: 7 × 8
##   effect   group          term       estimate std.error statistic    df  p.value
##   <chr>    <chr>          <chr>         <dbl>     <dbl>     <dbl> <dbl>    <dbl>
## 1 fixed    <NA>           (Intercep… 0.315       0.0366    8.61    5.37  2.41e-4
## 2 fixed    <NA>           Segment1   0.000497    0.0211    0.0236 12.2   9.82e-1
## 3 fixed    <NA>           Age1       0.0165      0.0366    0.451   5.37  6.69e-1
## 4 fixed    <NA>           Segment1:… 0.000562    0.0211    0.0266 12.2   9.79e-1
## 5 ran_pars Piglet         sd__(Inte… 0.123      NA        NA      NA    NA      
## 6 ran_pars Slaughter_date sd__(Inte… 0.0393     NA        NA      NA    NA      
## 7 ran_pars Residual       sd__Obser… 0.116      NA        NA      NA    NA
```

Pour regrouper les résultats de chaque phylum dans une même table on utilise la function `map_df()` qui crée un directement un data.frame.


```r
result_model <- df_stat %>%
  split(df_stat$phylum) %>% 
  map(fit_model) %>% 
  map_df(broom.mixed::tidy, .id = "phylum")
result_model
```

```
## # A tibble: 49 × 9
##    phylum        effect group term  estimate std.error statistic    df   p.value
##    <chr>         <chr>  <chr> <chr>    <dbl>     <dbl>     <dbl> <dbl>     <dbl>
##  1 Actinobacter… fixed  <NA>  (Int…  0.220     0.0145     15.1    21.7  5.34e-13
##  2 Actinobacter… fixed  <NA>  Segm… -0.0124    0.00974    -1.27   15.3  2.23e- 1
##  3 Actinobacter… fixed  <NA>  Age1   0.0130    0.0145      0.896  21.7  3.80e- 1
##  4 Actinobacter… fixed  <NA>  Segm…  0.00247   0.00974     0.254  15.3  8.03e- 1
##  5 Actinobacter… ran_p… Pigl… sd__…  0.0527   NA          NA      NA   NA       
##  6 Actinobacter… ran_p… Slau… sd__…  0        NA          NA      NA   NA       
##  7 Actinobacter… ran_p… Resi… sd__…  0.0541   NA          NA      NA   NA       
##  8 Bacteroidota  fixed  <NA>  (Int…  0.463     0.0310     15.0    22.6  3.37e-13
##  9 Bacteroidota  fixed  <NA>  Segm… -0.267     0.0115    -23.2    14.4  7.60e-13
## 10 Bacteroidota  fixed  <NA>  Age1   0.0269    0.0310      0.867  22.6  3.95e- 1
## # … with 39 more rows
```

### **Même méthode pour les p.value, shapiro et lettres**

**p.value**


```r
## fonction
p_val <- function(data)
{
  Anova(data, type = "III")
}

## résultats
result_pval <- df_stat %>%
  split(df_stat$phylum) %>% 
  map(fit_model) %>%
  map(p_val) %>% 
  map_dfr(tidy, .id = "phylum")
result_pval
```

```
## # A tibble: 28 × 5
##    phylum           term        statistic    df   p.value
##    <chr>            <chr>           <dbl> <dbl>     <dbl>
##  1 Actinobacteriota (Intercept)  229.         1 8.68e- 52
##  2 Actinobacteriota Segment        1.62       1 2.04e-  1
##  3 Actinobacteriota Age            0.803      1 3.70e-  1
##  4 Actinobacteriota Segment:Age    0.0645     1 7.99e-  1
##  5 Bacteroidota     (Intercept)  224.         1 1.52e- 50
##  6 Bacteroidota     Segment      540.         1 2.23e-119
##  7 Bacteroidota     Age            0.752      1 3.86e-  1
##  8 Bacteroidota     Segment:Age    0.172      1 6.79e-  1
##  9 Desulfobacterota (Intercept)  250.         1 2.48e- 56
## 10 Desulfobacterota Segment      168.         1 2.38e- 38
## # … with 18 more rows
```

**Ajout des p-adjusted**


```r
result_pval %>% 
  #dplyr::filter(effect == "fixed") %>% 
  dplyr::select(phylum, term, p.value) %>% 
  pivot_wider(names_from = term,
              values_from = p.value) %>% 
  mutate(p_adj_Segment = p.adjust(Segment, method = "BH"),
         p_adj_Age = p.adjust(Age, method = "BH"),
         p_adj_Age_Segment = p.adjust(`Segment:Age`, method = "BH")) -> df_pval
df_pval
```

```
## # A tibble: 7 × 8
##   phylum   `(Intercept)`   Segment     Age `Segment:Age` p_adj_Segment p_adj_Age
##   <chr>            <dbl>     <dbl>   <dbl>         <dbl>         <dbl>     <dbl>
## 1 Actinob…      8.68e-52 2.04e-  1 0.370          0.799      2.38e-  1    0.675 
## 2 Bactero…      1.52e-50 2.23e-119 0.386          0.679      1.56e-118    0.675 
## 3 Desulfo…      2.48e-56 2.38e- 38 0.788          0.775      8.33e- 38    0.788 
## 4 Firm_ba…      1.31e-24 5.42e- 22 0.563          0.715      1.27e- 21    0.760 
## 5 Firmicu…      0        2.22e- 14 0.254          0.264      3.89e- 14    0.675 
## 6 Fusobac…      1.24e- 6 9.45e-  2 0.00169        0.0302     1.32e-  1    0.0118
## 7 Proteob…      7.01e-18 9.81e-  1 0.652          0.979      9.81e-  1    0.760 
## # … with 1 more variable: p_adj_Age_Segment <dbl>
```

**Shapiro**


```r
## fonction
p_shap <- function(data)
{
  tmp <- summary(data)
  tmp2 <- shapiro.test(tmp$residuals)
}

## résultats
df_shap <- df_stat %>%
  split(df_stat$phylum) %>% 
  map(fit_model) %>%
  map(p_shap) %>% 
  map_dfr(tidy, .id = "phylum") %>% 
  dplyr::select(-method) %>% 
  rename(W = statistic,
         W_pval = p.value)
df_shap
```

```
## # A tibble: 7 × 3
##   phylum               W   W_pval
##   <chr>            <dbl>    <dbl>
## 1 Actinobacteriota 0.953 0.0948  
## 2 Bacteroidota     0.978 0.625   
## 3 Desulfobacterota 0.874 0.000363
## 4 Firm_bact_ratio  0.891 0.00103 
## 5 Firmicutes       0.941 0.0382  
## 6 Fusobacteriota   0.958 0.139   
## 7 Proteobacteria   0.923 0.00942
```

**Lettres**


```r
## fonction
p_letters <- function(data)
{
  tmp <- emmeans(data, pairwise ~ Segment * Age)
  tmp2 <- cld(tmp$emmeans, alpha = 0.05, Letters = letters, adjust ="tukey")
}

## résultats
result_letters <- df_stat %>%
  split(df_stat$phylum) %>% 
  map(fit_model) %>%
  map(p_letters) %>% 
  map_dfr(tidy, .id = "phylum")
result_letters
```

```
## # A tibble: 28 × 9
##    phylum       Segment Age   estimate std.error    df conf.low conf.high .group
##    <chr>        <chr>   <chr>    <dbl>     <dbl> <dbl>    <dbl>     <dbl> <chr> 
##  1 Actinobacte… Jejunum D35     0.192     0.0229  6.88  0.115      0.268  " a"  
##  2 Actinobacte… Colon   D35     0.222     0.0354 15.0   0.121      0.322  " a"  
##  3 Actinobacte… Jejunum D21     0.223     0.0229  6.74  0.146      0.300  " a"  
##  4 Actinobacte… Colon   D21     0.243     0.0238  7.58  0.165      0.320  " a"  
##  5 Bacteroidota Jejunum D35     0.175     0.0464  5.22  0.00242    0.347  " a " 
##  6 Bacteroidota Jejunum D21     0.219     0.0464  5.17  0.0460     0.392  " a " 
##  7 Bacteroidota Colon   D35     0.698     0.0562  9.37  0.526      0.871  "  b" 
##  8 Bacteroidota Colon   D21     0.762     0.0470  5.44  0.590      0.933  "  b" 
##  9 Desulfobact… Jejunum D21     0.0384    0.0200 10.6  -0.0216     0.0984 " a " 
## 10 Desulfobact… Jejunum D35     0.0386    0.0200 10.7  -0.0213     0.0985 " a " 
## # … with 18 more rows
```

**Mise en forme des résultats des lettres en colonnes**


```r
df_letters <- result_letters %>%
  dplyr::select(phylum, Segment, Age, .group) %>% 
  pivot_wider(names_from = c(Segment, Age), 
              values_from = .group)
df_letters
```

```
## # A tibble: 7 × 5
##   phylum           Jejunum_D35 Colon_D35 Jejunum_D21 Colon_D21
##   <chr>            <chr>       <chr>     <chr>       <chr>    
## 1 Actinobacteriota " a"        " a"      " a"        " a"     
## 2 Bacteroidota     " a "       "  b"     " a "       "  b"    
## 3 Desulfobacterota " a "       "  b"     " a "       "  b"    
## 4 Firm_bact_ratio  "  b"       " a "     "  b"       " a "    
## 5 Firmicutes       "  b d"     " a c "   "   cd"     " ab  "  
## 6 Fusobacteriota   " a "       " a "     " a "       "  b"    
## 7 Proteobacteria   " a"        " a"      " a"        " a"
```

### **Merge all results**


```r
left_join(df_shap, df_pval) %>% 
  left_join(df_letters) -> full_result
full_result
```

```
## # A tibble: 7 × 14
##   phylum               W   W_pval `(Intercept)`   Segment     Age `Segment:Age`
##   <chr>            <dbl>    <dbl>         <dbl>     <dbl>   <dbl>         <dbl>
## 1 Actinobacteriota 0.953 0.0948        8.68e-52 2.04e-  1 0.370          0.799 
## 2 Bacteroidota     0.978 0.625         1.52e-50 2.23e-119 0.386          0.679 
## 3 Desulfobacterota 0.874 0.000363      2.48e-56 2.38e- 38 0.788          0.775 
## 4 Firm_bact_ratio  0.891 0.00103       1.31e-24 5.42e- 22 0.563          0.715 
## 5 Firmicutes       0.941 0.0382        0        2.22e- 14 0.254          0.264 
## 6 Fusobacteriota   0.958 0.139         1.24e- 6 9.45e-  2 0.00169        0.0302
## 7 Proteobacteria   0.923 0.00942       7.01e-18 9.81e-  1 0.652          0.979 
## # … with 7 more variables: p_adj_Segment <dbl>, p_adj_Age <dbl>,
## #   p_adj_Age_Segment <dbl>, Jejunum_D35 <chr>, Colon_D35 <chr>,
## #   Jejunum_D21 <chr>, Colon_D21 <chr>
```

```r
datatable(full_result)
```

```{=html}
<div id="htmlwidget-d9cd19af78ac36d9f2c3" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-d9cd19af78ac36d9f2c3">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7"],["Actinobacteriota","Bacteroidota","Desulfobacterota","Firm_bact_ratio","Firmicutes","Fusobacteriota","Proteobacteria"],[0.952816431913969,0.978257055429324,0.87395600563975,0.890610644660262,0.941274799188378,0.957652847375552,0.922773968771566],[0.0948139591776235,0.625039810354336,0.000363131036589038,0.00103072631667422,0.0381823346558611,0.139160797423286,0.00942417147296379],[8.67902386441662e-52,1.51879223584256e-50,2.47639052015161e-56,1.31007350652692e-24,0,1.24096039370785e-06,7.00563748427384e-18],[0.20368154415844,2.22884598375759e-119,2.37893614847488e-38,5.4229034510239e-22,2.22262996197009e-14,0.0944989730731643,0.981193815184475],[0.370329339067513,0.385824079973006,0.788202948387599,0.562779217750694,0.253949991201482,0.00169223904251737,0.651682555154706],[0.799475377478624,0.678548595051256,0.774813553541423,0.715141158473788,0.263501323054347,0.0301846754262807,0.97875681540932],[0.237628468184846,1.56019218863031e-118,8.32627651966209e-38,1.26534413857224e-21,3.88960243344766e-14,0.13229856230243,0.981193815184475],[0.67519213995276,0.67519213995276,0.788202948387599,0.760296314347157,0.67519213995276,0.0118456732976216,0.760296314347157],[0.932721273725062,0.932721273725062,0.932721273725062,0.932721273725062,0.922254630690214,0.211292727983965,0.97875681540932],[" a"," a "," a ","  b","  b d"," a "," a"],[" a","  b","  b"," a "," a c "," a "," a"],[" a"," a "," a ","  b","   cd"," a "," a"],[" a","  b","  b"," a "," ab  ","  b"," a"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>phylum<\/th>\n      <th>W<\/th>\n      <th>W_pval<\/th>\n      <th>(Intercept)<\/th>\n      <th>Segment<\/th>\n      <th>Age<\/th>\n      <th>Segment:Age<\/th>\n      <th>p_adj_Segment<\/th>\n      <th>p_adj_Age<\/th>\n      <th>p_adj_Age_Segment<\/th>\n      <th>Jejunum_D35<\/th>\n      <th>Colon_D35<\/th>\n      <th>Jejunum_D21<\/th>\n      <th>Colon_D21<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


**Save in excel sheet**


```r
full_result
```

## **Remarque**

**Attention à l'ordre des lignes quand on fusionne des data.frame, éviter les `cbind()`**  

Utiliser la fonction `merge()` ou les fonctions `_join()` de dplyr qui garantissent le respect de l'ordre des lignes. La fusion s'effectue à partir d'une colonne commune aux 2 tables.

Exemple fusion des taxon et de la table d'abondance


```r
TA <- data.frame(otu_table(relatabunphy))
TAXA <- tax_table(tab_Phylum)[,"Phylum"]
```

**merge()**


```r
## merge()
temp_merge <- merge(TAXA , TA, by = "row.names")
```

**left_join()**

N'accepte que des data.frames
On utilise la fonction rownames_to_column() du package tibble qui transforme les rownames en colonne avec "rowname" en en-tête par défaut.
left_join peut repèrer automatiquement la colonne en commun: "rowname"  


```r
## left_join()
temp_join <- left_join(rownames_to_column(data.frame(TAXA)), rownames_to_column(TA))
```


```r
## par défaut merge() tri sur la colonne en commun
temp_merge[1:5,1:5]
```

```
##     Row.names           Phylum GPS158690_CCTTGA.JLGMN_L001
## 1   Cluster_1       Firmicutes                0.9984080886
## 2 Cluster_189 Campylobacterota                0.0000000000
## 3 Cluster_294  Patescibacteria                0.0000000000
## 4  Cluster_31     Bacteroidota                0.0002796601
## 5  Cluster_57   Fusobacteriota                0.0000000000
##   GPS158692_TCGTTC.JLGMN_L001 GPS158694.PCR450.6E_GTTTCT.JLGMN_L001
## 1                7.877842e-01                            0.88554166
## 2                9.093665e-05                            0.00000000
## 3                0.000000e+00                            0.00000000
## 4                1.965444e-01                            0.08483679
## 5                0.000000e+00                            0.00000000
```

```r
temp_join[1:5,1:5]
```

```
##       rowname           Phylum GPS158690_CCTTGA.JLGMN_L001
## 1  Cluster_31     Bacteroidota                0.0002796601
## 2 Cluster_189 Campylobacterota                0.0000000000
## 3   Cluster_1       Firmicutes                0.9984080886
## 4  Cluster_97 Desulfobacterota                0.0000000000
## 5  Cluster_86 Actinobacteriota                0.0004087340
##   GPS158692_TCGTTC.JLGMN_L001 GPS158694.PCR450.6E_GTTTCT.JLGMN_L001
## 1                1.965444e-01                            0.08483679
## 2                9.093665e-05                            0.00000000
## 3                7.877842e-01                            0.88554166
## 4                3.485905e-03                            0.00268524
## 5                8.487420e-04                            0.01476882
```

**Session Info**


```r
sessionInfo()
```

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.3 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
## 
## locale:
##  [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
##  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8   
##  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] ggplot2_3.3.5     DT_0.20           emmeans_1.7.2     lmerTest_3.1-3   
##  [5] lme4_1.1-28       Matrix_1.2-18     multcomp_1.4-18   TH.data_1.1-0    
##  [9] MASS_7.3-51.5     survival_3.1-8    mvtnorm_1.1-3     car_3.0-12       
## [13] carData_3.0-5     phyloseq_1.30.0   broom.mixed_0.2.7 broom_0.7.12     
## [17] purrr_0.3.4       tidyr_1.2.0       dplyr_1.0.8       tibble_3.1.6     
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-144        pbkrtest_0.5.1      numDeriv_2016.8-1.1
##  [4] tools_3.6.3         backports_1.4.1     bslib_0.3.1        
##  [7] utf8_1.2.2          R6_2.5.1            vegan_2.5-7        
## [10] DBI_1.1.2           BiocGenerics_0.32.0 mgcv_1.8-31        
## [13] colorspace_2.0-2    permute_0.9-7       ade4_1.7-18        
## [16] withr_2.4.3         tidyselect_1.1.1    compiler_3.6.3     
## [19] cli_3.2.0           Biobase_2.46.0      sandwich_3.0-1     
## [22] sass_0.4.0          scales_1.1.1        multcompView_0.1-8 
## [25] stringr_1.4.0       digest_0.6.29       minqa_1.2.4        
## [28] rmarkdown_2.11      XVector_0.26.0      pkgconfig_2.0.3    
## [31] htmltools_0.5.2     fastmap_1.1.0       htmlwidgets_1.5.4  
## [34] rlang_1.0.1         rstudioapi_0.13     jquerylib_0.1.4    
## [37] generics_0.1.2      zoo_1.8-9           jsonlite_1.7.3     
## [40] crosstalk_1.2.0     magrittr_2.0.2      biomformat_1.14.0  
## [43] Rcpp_1.0.8          munsell_0.5.0       S4Vectors_0.24.4   
## [46] Rhdf5lib_1.8.0      fansi_1.0.2         ape_5.6-1          
## [49] abind_1.4-5         lifecycle_1.0.1     stringi_1.7.6      
## [52] yaml_2.2.2          zlibbioc_1.32.0     rhdf5_2.30.1       
## [55] plyr_1.8.6          grid_3.6.3          parallel_3.6.3     
## [58] crayon_1.5.0        lattice_0.20-40     Biostrings_2.54.0  
## [61] splines_3.6.3       multtest_2.42.0     knitr_1.37         
## [64] pillar_1.7.0        igraph_1.2.11       boot_1.3-24        
## [67] estimability_1.3    reshape2_1.4.4      codetools_0.2-16   
## [70] stats4_3.6.3        glue_1.6.1          evaluate_0.14      
## [73] data.table_1.14.2   vctrs_0.3.8         nloptr_2.0.0       
## [76] foreach_1.5.2       gtable_0.3.0        assertthat_0.2.1   
## [79] xfun_0.29           xtable_1.8-4        coda_0.19-4        
## [82] iterators_1.0.14    IRanges_2.20.2      cluster_2.1.0      
## [85] ellipsis_0.3.2
```


