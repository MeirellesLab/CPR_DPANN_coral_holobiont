Holobiont project KW
================
Arthur & Leticia
12/10/2019

# Kruskal-Wallis Univariate test

## Libraries we used

``` r
library(janitor)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(broom)
```

## Importing and tidying the data

``` r
data <- read.csv("input/detailed_abundance_phyla_157_corais_matrix_porcentagem_prot_divided_17abr2019.csv",
                 check.names = FALSE) %>% 
  select(-1)

levels(data$health_status)[levels(data$health_status) == "iminent bleaching"] <- "bleached"
levels(data$health_status)[levels(data$health_status) == "signs of bleaching"] <- "bleached"
levels(data$health_status)[levels(data$health_status) == "bleaching event"] <- "bleached" 
levels(data$health_status)[levels(data$health_status) == "bleaching event "] <- "bleached"
levels(data$health_status)[levels(data$health_status) == "Cyanobacterial patches"] <- "Black Band Disease"


data_num <- data %>% select_if(is.numeric)

abundantes <- data_num %>% gather() %>% group_by(key) %>% 
  summarise(mean_abund = mean(value)) %>% 
  filter(mean_abund > 1/100) %>% pull(key)


data <- data %>% 
  mutate(
    seq_method = factor(
      case_when(
        .$seq_method %in% c("Illumina Genome Analyzer IIx", "Illumina HiSeq 1000",
                            "Illumina HiSeq 1000 WGS", "Illumina HiSeq 2000",
                            "Illumina HiSeq 2500", "Illumina HiSeq 2500 WGS",
                            "WGS illumina") ~ "Illumina",
        TRUE ~ "454"
      )
    )
  )

data <- data %>% filter(genus != "Fungia",
                        genus != "Pseudodiploria", 
                        !(.$genus == "Porites" & .$authors == "Thurber et al. 2009"),
                        !(.$genus == "Acropora" & .$region == "mesocosmos" & .$authors == "not published"),
                        !(.$genus == "Pocillopora" & .$region == "mesocosmos" & .$authors =="not published"),
                        !(.$genus == "Seriatopora" &.$region == "mesocosmos" & .$authors =="not published"),
                        !(.$id == "mgm4694757_prinseq_good_Crf1")) %>% 
  mutate(genus = factor(genus),
         authors = factor(authors),
         growth_form = factor(growth_form))

data_num <- data %>% select_if(is.numeric)
```

## Running the test

``` r
pvalue <- matrix(0, nrow=ncol(data_num), ncol=length(2:8)) %>% as.data.frame()
statistic <- pvalue
d_f <- pvalue

for(i in 1:ncol(data_num)){
  
  for(k in 2:8){
    
  le_test <- kruskal.test(data_num[, i] ~ data[, k])
  
  pvalue[i, k] <- le_test$p.value
  statistic[i, k] <- le_test$statistic
  }
  
}
```

## Tidying results

``` r
df_kw_individual <- bind_cols(
  pvalue,
  statistic,
  as.data.frame(colnames(data_num))
) %>%
  select(-c(V1, V11))

colnames(df_kw_individual) <- c(paste0("pvalue_", colnames(data)[2:8]),
                                paste0(paste0("statistic_", colnames(data)[2:8])),
                                "taxa") 

df_kw_individual <- df_kw_individual %>% 
  mutate(taxa = str_replace_all(taxa, "\\.", " ") %>% 
           str_replace("Candidatus", "Ca.") %>% 
           str_replace("candidate", "Ca.") %>% 
           str_replace("Unclassified class", "Un. class"))

final_table <- df_kw_individual %>% select(contains("pvalue")) %>% gather() %>% 
  mutate(significant = case_when(.$value >= 0.05 ~ 0,
                                 TRUE            ~ 1)
         ) %>% 
  group_by(key) %>% 
  summarise(Significance_ratio = paste0(
    round(sum(significant)*100/ncol(data_num),2), " %")
    ) %>% 
  arrange(desc(Significance_ratio)) %>% 
  mutate(Variable = str_replace(key, "pvalue_","") %>% 
           str_replace("_"," ")) %>% 
  select(Variable, Significance_ratio)

write_csv(x = final_table, path = "outputs/KW/univariate_test.csv")
```

## Ploting results

``` r
df_kw_individual %>% select(contains("pvalue"), taxa) %>% gather(key, value, -taxa) %>% 
  mutate(
    key = str_replace(key, "pvalue_","") %>% str_replace("_", " "),
    significant = case_when(.$value >= 0.01 ~ 0,
                            TRUE            ~ 1)
  ) %>%
  ggplot(aes(x = reorder(key, significant), y = reorder(taxa, significant),
             fill = factor(significant)))+
  geom_tile(colour = "grey90")+
  scale_fill_manual(values = c("1" = "gray70", "0" = "firebrick"),
                    labels = c("No", "Yes"))+
  labs(x = NULL, y = NULL, fill = "Significant with\nP-value = 0.01")+
  theme_minimal()+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90,hjust = 1, size = 14))+coord_flip()+
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(face="bold", size=14),
        #panel.grid.major.y = element_line(colour = "black")
        )
```

![](KW_test_files/figure-gfm/ploting%20signif%20.01-1.png)<!-- -->

``` r
ggsave("outputs/KW/heatmap_KW_sig_01.png", width = 160, height = 80, units = "mm", scale = 4)
ggsave("outputs/KW/heatmap_KW_sig_01.pdf", width = 160, height = 80, units = "mm", scale = 4)
ggsave("outputs/KW/heatmap_KW_sig_01.svg", width = 160, height = 80, units = "mm", scale = 4)
```

``` r
df_kw_individual %>% select(contains("pvalue"), taxa) %>% gather(key, value, -taxa) %>% 
  mutate(
    key = str_replace(key, "pvalue_","") %>% str_replace("_", " "),
    significant = case_when(.$value >= 0.05 ~ 0,
                            TRUE            ~ 1)
  ) %>%
  ggplot(aes(x = reorder(key, significant), y = reorder(taxa, significant),
             fill = factor(significant)))+
  geom_tile(colour = "grey90")+
  scale_fill_manual(values = c("0" = "#FA8716", "1" = "#3B2005"),
                    labels = c("No", "Yes"))+
  labs(x = NULL, y = NULL, fill = "Significant with\nP-value = 0.05")+
  theme_minimal()+
  theme(legend.position = "top", axis.text.x = element_text(angle = 90,hjust = 1, size = 14,
                                                            vjust = 0.5))+coord_flip()+
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.y = element_text(size = 12))+
  theme(axis.title.y = element_text(face="bold", size=14),
        #panel.grid.major.y = element_line(colour = "black")
        )+
  annotation_custom(tableGrob(final_table, rows=NULL), xmin = 3.89, ymin = 125)
```

![](KW_test_files/figure-gfm/ploting%20signif%20.05-1.png)<!-- -->

``` r
ggsave("outputs/KW/heatmap_KW_sig_05.png", width = 160, height = 80, units = "mm", scale = 4)
ggsave("outputs/KW/heatmap_KW_sig_05.pdf", width = 160, height = 80, units = "mm", scale = 4)
ggsave("outputs/KW/heatmap_KW_sig_05.svg", width = 160, height = 80, units = "mm", scale = 4)
```

# Kruskal-Wallis MuLtivariate test

## Separating rare and abundant subjects

``` r
multkw <- function(group,y,simplify=FALSE){
  
  ### sort data by group ###
  o<-order(group)
  group<-group[o]
  y<-as.matrix(y[o,])
  n<-length(group)
  p<-dim(y)[2]
  if (dim(y)[1] != n)
    return("number of oberservations not equal to length of group")
  
  groupls<-unique(group)
  g<-length(groupls) #number of groups#
  groupind<-sapply(groupls,"==",group) #group indicator#
  ni<-colSums(groupind) #num of subj of each group#
  r<-apply(y,2,rank) #corresponding rank variable#
  
  ### calculation of statistic ###
  r.ik<-t(groupind)%*%r*(1/ni) #gxp, mean rank of kth variate in ith
  group#
  m<- (n+1)/2 #expected value of rik#
  u.ik<-t(r.ik-m)
  U<-as.vector(u.ik)
  V<-1/(n-1)*t(r-m)%*%(r-m) #pooled within-group cov matrix
  Vstar<-bdiag(lapply(1/ni,"*",V))
  W2<-as.numeric(t(U)%*%solve(Vstar)%*%U)
  
  ### return stat and p-value ###
  returnlist<-list(statistic=W2,d.f.=p*(g-1),
                   p.value=pchisq(W2,p*(g-1),lower.tail=F))
  if (simplify==TRUE) return (W2)
  else return (returnlist)
}
```

``` r
abundantes <- data_num %>% gather() %>% group_by(key) %>% 
  summarise(mean_abund = mean(value)) %>% 
  filter(mean_abund > 1/100) %>% pull(key)

raro <- data_num %>% gather() %>% group_by(key) %>% 
  summarise(mean_abund = mean(value)) %>% 
  filter(mean_abund < 1/100) %>% pull(key)
```

## Tidying

``` r
mkw_abun_genus <- multkw(y = select(data_num, abundantes), group = data$genus) %>% 
  unlist() %>% t() %>% as.data.frame()
mkw_abun_health_status <- multkw(y = select(data_num, abundantes), group = data$health_status) %>% 
  unlist() %>% t() %>% as.data.frame()
mkw_abun_growth_form <- multkw(y = select(data_num, abundantes), group = data$growth_form) %>% 
  unlist() %>% t() %>% as.data.frame()
mkw_abun_seq_method <- multkw(y = select(data_num, abundantes), group = data$seq_method) %>% 
  unlist() %>% t() %>% as.data.frame()
mkw_abun_region <- multkw(y = select(data_num, abundantes), group = data$region) %>% 
  unlist() %>% t() %>% as.data.frame()
mkw_abun_authors <- multkw(y = select(data_num, abundantes), group = data$authors) %>% 
  unlist() %>% t() %>% as.data.frame()

mkw_abun_geral <- bind_rows(mkw_abun_genus, mkw_abun_health_status,
                            mkw_abun_growth_form,mkw_abun_seq_method, 
                            mkw_abun_region, mkw_abun_authors)
mkw_abun_geral$var <- c("genus", "health_status", "growth_form",
                        "seq_method", "region", "authors")
```

## Ploting

``` r
mkw_plot <- mkw_abun_geral %>% 
  mutate(var = fct_reorder(var, -statistic)) %>% 
  ggplot(aes(x = var, y = statistic))+
  labs(x = NULL, y = "Statistic")+
  geom_col(position = "dodge", fill = "#B39266")+
  scale_x_discrete(labels = c("seq_method" = "Sequency Method", "growth_form" = "Growth form",
                              "health_status" = "Health status",
                              "region" = "Region", "genus" = "Genus", "authors" = "Authors"))+
  theme_minimal()+
  theme(axis.title.x = element_text(face="bold", size=14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))+
  theme(axis.title.y = element_text(face="bold", size=14))

ggsave(mkw_plot, file="outputs/KW/MKW_sep_abundant_rare.jpg", width = 160, height = 80, units = "mm", scale = 4)
ggsave(mkw_plot, file="outputs/KW/MKW_sep_abundant_rare.png", width = 160, height = 80, units = "mm", scale = 4)
ggsave(mkw_plot, file="outputs/KW/MKW_sep_abundant_rare.pdf", width = 160, height = 80, units = "mm", scale = 4)
ggsave(mkw_plot, file="outputs/KW/MKW_sep_abundant_rare.svg", width = 160, height = 80, units = "mm", scale = 4)

mkw_plot
```

![](KW_test_files/figure-gfm/ploting%20-1.png)<!-- -->
