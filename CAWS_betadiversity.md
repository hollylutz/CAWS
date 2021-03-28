
### Load them libraries...


```R
library(phyloseq)
library(ggplot2)
require(scales)
library(vegan)
library(readr)
library(data.table)
library(dplyr)
library(stats)
library(ggpubr)
library(dunn.test)
```

#### Read in Sho original phyloseq object


```R
#Non-transformed OTU counts
ps.nt = readRDS("/home/holutz/projects/CAWS/phyloseq.rds") 
```


```R
#New ps object with 2019 water data
ps2.nt = readRDS("/home/holutz/projects/CAWS/phyloseq2.rds")
#Transform ASV counts to relative abundance
ps2.t = ps2.nt %>%
     transform_sample_counts(function(x) x / sum(x))
```

#### Introduce some modifications to original phyloseq object


```R
#Import phylogenetic information
tre = read_tree("/home/holutz/projects/CAWS/study_tree.tree")

#Import modified metadata w numeric variables as factors
map = import_qiime_sample_data("/home/holutz/projects/CAWS/Metadata_2019_corrected_water2.txt")
map$year_coll.factor = as.factor(map$year_coll.factor)
map$SITE_CODE.factor = as.factor(map$SITE_CODE.factor)
#map = sample_data(ps2.nt)

#breakdown ps.nt into components to create new phyloseq object
otu = otu_table(ps2.nt)
tax = tax_table(ps2.nt)

#create new phyloseq object by combining all of the above
ps.nt = merge_phyloseq(tre,map,otu,tax)

#Trasform ASVs for relative abundance within each library
ps.t = ps.nt %>%
    transform_sample_counts(function(x) x / sum(x))
                            
#Save new non-transformed and transformed phyloseq objects for subsequent analysis
saveRDS(ps.nt, "/home/holutz/projects/CAWS/CAWS_nt.2.rds")
saveRDS(ps.t, "/home/holutz/projects/CAWS/CAWS_t.2.rds")
```

### BEGIN ANALYSES HERE


```R
#Read in phyloseq objects (ps)
ps.nt = readRDS("/home/holutz/projects/CAWS/CAWS_nt.2.rds")
ps.t = readRDS("/home/holutz/projects/CAWS/CAWS_t.2.rds")
#Display ps information
ps.t
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 50223 taxa and 832 samples ]
    sample_data() Sample Data:       [ 832 samples by 42 sample variables ]
    tax_table()   Taxonomy Table:    [ 50223 taxa by 7 taxonomic ranks ]
    phy_tree()    Phylogenetic Tree: [ 50223 tips and 50201 internal nodes ]



```R
#Explore metadata as a dataframe
ps.df = data.frame(sample_data(ps.t))
colnames(ps.df)
head(ps.df$SITE_CODE.factor)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'X.SampleID'</li><li>'SITE_CODE'</li><li>'SITE_CODE.factor'</li><li>'Sample_type'</li><li>'SourceSink'</li><li>'Env'</li><li>'kit'</li><li>'Disinfection_status'</li><li>'date_collected'</li><li>'year_coll.factor'</li><li>'year_collected'</li><li>'monthyear'</li><li>'season_collected'</li><li>'month_collected'</li><li>'Herbert'</li><li>'Weather'</li><li>'Location'</li><li>'Land_use'</li><li>'run_prefix'</li><li>'DO'</li><li>'Temp'</li><li>'pH'</li><li>'NO2_NO3'</li><li>'NH3_N'</li><li>'TOT_PHOS'</li><li>'SO4'</li><li>'TDS'</li><li>'SS'</li><li>'ALK'</li><li>'Cl'</li><li>'Fluorine'</li><li>'TOC'</li><li>'FEC_COL'</li><li>'Chlorophyll'</li><li>'conductivity'</li><li>'Turbidity'</li><li>'Ref_to_WRP'</li><li>'Miles_from_WRP'</li><li>'Region'</li><li>'Region_type'</li><li>'Disinfection'</li><li>'Description'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>SiteCode_56</li><li>SiteCode_56</li><li>SiteCode_59</li><li>SiteCode_59</li><li>SiteCode_99</li><li>SiteCode_52</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'SiteCode_100'</li><li>'SiteCode_108'</li><li>'SiteCode_112'</li><li>'SiteCode_36'</li><li>'SiteCode_43'</li><li>'SiteCode_52'</li><li>'SiteCode_55'</li><li>'SiteCode_56'</li><li>'SiteCode_57'</li><li>'SiteCode_59'</li><li>'SiteCode_73'</li><li>'SiteCode_75'</li><li>'SiteCode_76'</li><li>'SiteCode_86'</li><li>'SiteCode_96'</li><li>'SiteCode_97'</li><li>'SiteCode_99'</li></ol>
</details>


### Calculate betadiversity metrics for PCoA

#### Weighted and Unweighted UniFrac for all samples combined


```R
#Wunif
#wunif.dist = distance(ps.t, method = "wunifrac", type = "samples")
#wunif.dist.log = ordinate(ps.t, method = "PCoA", distance = "wunifrac")
#wunif.dist.evals = wunif.dist.log$values$Eigenvalues

#Save estimates so you don't have to rerun every time
#saveRDS(wunif.dist, "/home/holutz/projects/CAWS/wunifrac/wunif.ps2.dist.rds")
#saveRDS(wunif.dist.log, "/home/holutz/projects/CAWS/wunifrac/wunif.ps2.dist.log.rds")
#saveRDS(wunif.dist.evals, "/home/holutz/projects/CAWS/wunifrac/wunif.ps2.dist.evals.rds")
```


```R
#Unif

#unif.dist = distance(ps.t, method = "unifrac", type = "samples")
#unif.dist.log = ordinate(ps.t, method = "PCoA", distance = "unifrac")
#unif.dist.evals = unif.dist.log$values$Eigenvalues

#Save estimates so you don't have to rerun every time
#saveRDS(unif.dist, "/home/holutz/projects/CAWS/unifrac/unif.ps2.dist.rds")
#saveRDS(unif.dist.log, "/home/holutz/projects/CAWS/unifrac/unif.ps2.dist.log.rds")
#saveRDS(unif.dist.evals, "/home/holutz/projects/CAWS/unifrac/unif.ps2.dist.evals.rds")
```


```R
#Read in saved distance metrics for plotting and statistical analysis
#Weighted UniFrac
wunif.dist = readRDS("/home/holutz/projects/CAWS/wunifrac/wunif.ps2.dist.rds")
wunif.dist.log = readRDS("/home/holutz/projects/CAWS/wunifrac/wunif.ps2.dist.log.rds")
wunif.dist.evals = readRDS("/home/holutz/projects/CAWS/wunifrac/wunif.ps2.dist.evals.rds")

#Unweighted UniFrac
unif.dist = readRDS("/home/holutz/projects/CAWS/unifrac/unif.ps2.dist.rds")
unif.dist.log = readRDS("/home/holutz/projects/CAWS/unifrac/unif.ps2.dist.log.rds")
unif.dist.evals = readRDS("/home/holutz/projects/CAWS/unifrac/unif.ps2.dist.evals.rds")
```

#### Plot PCoA


```R
library(viridis)

p = plot_ordination(ps.t, wunif.dist.log, color = "year_coll.factor") + #, shape = "SampleType", label="X.SampleID") + 
    ggtitle("Weighted UniFrac by SITE_CODE and year_collected") +
    geom_point(size=1, alpha=0.8) + 
    scale_colour_viridis(option="viridis", discrete=TRUE, direction=1) +
    #stat_ellipse(level=0.95, geom="polygon",alpha = .2, aes(fill = year_collected), linetype=0) +
    #scale_fill_viridis(option="viridis", discrete=FALSE, direction=1)+
    coord_fixed(sqrt(wunif.dist.evals[2] / wunif.dist.evals[1]))

p = p + guides(colour = guide_legend(title = "year_collected", ncol = 2, keywidth = 1, keyheight = 1))

p = p + theme_minimal() + 
        theme(panel.grid.major = element_blank(),
              axis.title.y=element_text(margin=margin(0,20,0,0)),
              text=element_text(size=8, color="black",family="Arial"),
              panel.grid.minor = element_blank(),
              #axis.ticks = element_blank(),
              legend.key.size = unit(.1, "in"),
              legend.spacing.x = unit(.01,"in"),
              legend.position = "bottom",
              plot.title = element_text(hjust = 0.5))

p$layers <- p$layers[-1]

p

p2 = p + facet_wrap(~SITE_CODE.factor)

p2

png('/home/holutz/projects/CAWS/wunifrac/wunif_by_year_v2.png', width=8, height=12, units='in', res=500)
plot(p)
dev.off()

png('/home/holutz/projects/CAWS/wunifrac/wunif_by_sitecode.facet_v2.png', width=8, height=8, units='in', res=500)
plot(p2)
dev.off()
```


![png](output_16_0.png)



<strong>png:</strong> 2



<strong>png:</strong> 2



![png](output_16_3.png)



```R
sample_tab = data.frame(sample_data(ps.t))
head(sample_tab)
```


<table>
<caption>A data.frame: 6 × 93</caption>
<thead>
	<tr><th></th><th scope=col>X.SampleID</th><th scope=col>SITE_CODE</th><th scope=col>SITE_CODE.factor</th><th scope=col>Sample_type</th><th scope=col>SourceSink</th><th scope=col>Env</th><th scope=col>kit</th><th scope=col>Disinfection_status</th><th scope=col>date_collected</th><th scope=col>year_collected</th><th scope=col>⋯</th><th scope=col>Total_Length</th><th scope=col>Standard_Length</th><th scope=col>Weight</th><th scope=col>Sex</th><th scope=col>Final.Hg</th><th scope=col>FEC_COL_pre</th><th scope=col>Hg</th><th scope=col>Description</th><th scope=col>year_coll.factor.1</th><th scope=col>SITE_CODE.factor.1</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>⋯</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;lgl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>3773</th><td>3773</td><td>96</td><td>SC96</td><td>Water</td><td>NA</td><td>NA</td><td>NA</td><td>Post</td><td>9/9/19 </td><td>2019</td><td>⋯</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td></td><td>2019</td><td>96</td></tr>
	<tr><th scope=row>3710</th><td>3710</td><td>99</td><td>SC99</td><td>Water</td><td>NA</td><td>NA</td><td>NA</td><td>Post</td><td>6/17/19</td><td>2019</td><td>⋯</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td></td><td>2019</td><td>99</td></tr>
	<tr><th scope=row>3581</th><td>3581</td><td>99</td><td>SC99</td><td>Water</td><td>NA</td><td>NA</td><td>NA</td><td>Post</td><td>3/18/19</td><td>2019</td><td>⋯</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td></td><td>2019</td><td>99</td></tr>
	<tr><th scope=row>3703</th><td>3703</td><td>96</td><td>SC96</td><td>Water</td><td>NA</td><td>NA</td><td>NA</td><td>Post</td><td>6/10/19</td><td>2019</td><td>⋯</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td></td><td>2019</td><td>96</td></tr>
	<tr><th scope=row>3621</th><td>3621</td><td>34</td><td>SC34</td><td>Water</td><td>NA</td><td>NA</td><td>NA</td><td>Post</td><td>5/1/19 </td><td>2019</td><td>⋯</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td></td><td>2019</td><td>34</td></tr>
	<tr><th scope=row>3781</th><td>3781</td><td>86</td><td>SC86</td><td>Water</td><td>NA</td><td>NA</td><td>NA</td><td>Post</td><td>9/23/19</td><td>2019</td><td>⋯</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td><td></td><td>2019</td><td>86</td></tr>
</tbody>
</table>




```R
#Estimate betadispersion - weighted unifrac
beta.disp.wuf = betadisper(wunif.dist,group = sample_data(ps.t)$SITE_CODE.factor)
beta.disp.wuf = beta.disp.wuf$distances
beta.disp.wuf = data.frame(d = beta.disp.wuf, X.SampleID = names(beta.disp.wuf)) #The "X." in SampleID is introduced because of the # sign in the metadata
beta.disp.wuf2 = left_join(beta.disp.wuf, sample_tab, by="X.SampleID")
```


```R
kruskal.test(d~SITE_CODE.factor, data=beta.disp.wuf2)

krusk.wuf = dunn.test(beta.disp.wuf2$d, beta.disp.wuf2$SITE_CODE.factor, method = "Bonferroni")
write.table(krusk.wuf, "/home/holutz/projects/CAWS/unifrac/kw.wuf.bdisp.txt", sep="\t")
```


    
    	Kruskal-Wallis rank sum test
    
    data:  d by SITE_CODE.factor
    Kruskal-Wallis chi-squared = 155.41, df = 17, p-value < 2.2e-16



      Kruskal-Wallis rank sum test
    
    data: x and group
    Kruskal-Wallis chi-squared = 155.412, df = 17, p-value = 0
    
    
                               Comparison of x by group                            
                                     (Bonferroni)                                  
    Col Mean-|
    Row Mean |      SC100      SC108      SC112       SC34       SC36       SC43
    ---------+------------------------------------------------------------------
       SC108 |  -1.666786
             |     1.0000
             |
       SC112 |  -0.359618   1.290544
             |     1.0000     1.0000
             |
        SC34 |   1.641248   2.052470   1.728225
             |     1.0000     1.0000     1.0000
             |
        SC36 |   0.358727   2.022760   0.713962  -1.554910
             |     1.0000     1.0000     1.0000     1.0000
             |
        SC43 |   1.229685   2.426420   1.480074  -1.194014   0.974632
             |     1.0000     1.0000     1.0000     1.0000     1.0000
             |
        SC52 |   0.197260   1.460813   0.467707  -1.543085  -0.071576  -0.880611
             |     1.0000     1.0000     1.0000     1.0000     1.0000     1.0000
             |
        SC55 |   2.860253   4.053865   3.100526  -0.674208   2.604443   1.316108
             |     0.3238    0.0039*     0.1478     1.0000     0.7040     1.0000
             |
        SC56 |   4.382983   5.991792   4.680873  -0.614972   4.029352   1.847621
             |    0.0009*    0.0000*    0.0002*     1.0000    0.0043*     1.0000
             |
        SC57 |   0.456148   2.166004   0.820129  -1.537294   0.086739  -0.928644
             |     1.0000     1.0000     1.0000     1.0000     1.0000     1.0000
             |
        SC59 |   3.893582   5.510537   4.198994  -0.726967   3.539042   1.507798
             |    0.0076*    0.0000*    0.0021*     1.0000     0.0307     1.0000
             |
        SC73 |   3.418250   4.997191   3.723103  -0.808835   3.073632   1.234763
             |     0.0482    0.0000*    0.0151*     1.0000     0.1618     1.0000
             |
        SC75 |   3.395522   4.075841   3.535195   0.330516   3.251642   2.508868
             |     0.0524    0.0035*     0.0312     1.0000     0.0878     0.9266
             |
        SC76 |   4.260968   5.864053   4.559544  -0.637143   3.908962   1.774656
             |    0.0016*    0.0000*    0.0004*     1.0000    0.0071*     1.0000
             |
        SC86 |   0.143133   1.845712   0.509042  -1.609646  -0.224020  -1.143538
             |     1.0000     1.0000     1.0000     1.0000     1.0000     1.0000
             |
        SC96 |   3.291517   4.867890   3.597385  -0.836340   2.947740   1.150449
             |     0.0762    0.0001*    0.0246*     1.0000     0.2449     1.0000
             |
        SC97 |   1.113369   2.347358   1.372095  -1.247229   0.850201  -0.125580
             |     1.0000     1.0000     1.0000     1.0000     1.0000     1.0000
             |
        SC99 |  -3.152277  -1.422913  -2.746450  -2.410629  -3.519859  -3.499472
             |     0.1239     1.0000     0.4609     1.0000     0.0330     0.0357
    Col Mean-|
    Row Mean |       SC52       SC55       SC56       SC57       SC59       SC73
    ---------+------------------------------------------------------------------
        SC55 |   2.247193
             |     1.0000
             |
        SC56 |   3.044817   0.237750
             |     0.1781     1.0000
             |
        SC57 |   0.136770  -2.584012  -4.072862
             |     1.0000     0.7471    0.0036*
             |
        SC59 |   2.686746  -0.105101  -0.493464   3.565889
             |     0.5520     1.0000     1.0000     0.0277
             |
        SC73 |   2.377314  -0.353608  -0.829435   3.081331  -0.352950
             |     1.0000     1.0000     1.0000     0.1576     1.0000
             |
        SC75 |   3.095931   1.678374   1.675248   3.233305   1.863740   1.993393
             |     0.1501     1.0000     1.0000     0.0936     1.0000     1.0000
             |
        SC76 |   2.964023   0.168697  -0.098219   3.947399   0.392287   0.729956
             |     0.2323     1.0000     1.0000    0.0060*     1.0000     1.0000
             |
        SC86 |  -0.093538  -2.794216  -4.351625  -0.319973  -3.848977  -3.360031
             |     1.0000     0.3980    0.0010*     1.0000    0.0091*     0.0596
             |
        SC96 |   2.286516  -0.435949  -0.943182   2.951219  -0.468344  -0.112930
             |     1.0000     1.0000     1.0000     0.2421     1.0000     1.0000
             |
        SC97 |   0.770386  -1.470808  -2.064481   0.801174  -1.713301  -1.427120
             |     1.0000     1.0000     1.0000     1.0000     1.0000     1.0000
             |
        SC99 |  -2.574273  -5.158816  -7.606402  -3.712000  -7.111350  -6.532460
             |     0.7684    0.0000*    0.0000*    0.0157*    0.0000*    0.0000*
    Col Mean-|
    Row Mean |       SC75       SC76       SC86       SC96       SC97
    ---------+-------------------------------------------------------
        SC76 |  -1.711383
             |     1.0000
             |
        SC86 |  -3.351925  -4.225859
             |     0.0614    0.0018*
             |
        SC96 |  -2.038656  -0.843439   3.230099
             |     1.0000     1.0000     0.0947
             |
        SC97 |  -2.609705  -1.988270   1.023416  -1.339714
             |     0.6932     1.0000     1.0000     1.0000
             |
        SC99 |  -4.686351  -7.465874  -3.371398  -6.394791  -3.453384
             |    0.0002*    0.0000*     0.0572    0.0000*     0.0424
    
    alpha = 0.05
    Reject Ho if p <= alpha/2



```R
avg_obs <- aggregate(d ~ SITE_CODE.factor, data = beta.disp.wuf2, mean)
beta.disp.wuf2$SITE_CODE.factor <-factor(beta.disp.wuf2$SITE_CODE.factor, levels=avg_obs[order(avg_obs$d), "SITE_CODE.factor"])
```


```R
p = ggplot(beta.disp.wuf2, aes(x = year_coll.factor, y = d, fill = year_coll.factor)) + ggtitle("Weighted Unifrac betadispersion")

p = p + geom_boxplot(outlier.stroke = 0.2, outlier.shape=20) + 
        ylab("Distance of samples from group centroid ") + 
        xlab("Year") #+ 
        #scale_y_continuous(limits=c(0,0.2))
        #coord_flip(ylim = c(0,1))

p = p + scale_fill_viridis(discrete=TRUE, alpha=0.8)

p = p + theme_bw() + theme(text=element_text(size=12, color="black",family = "Arial"),
                           legend.key = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           legend.position="none",
                           axis.text.x=element_text(angle=45,vjust=0.4))
                           #axis.title.x=element_blank(),
                           #axis.title.y=element_text(margin=margin(0,10,0,0)))

#p = p + geom_hline(yintercept = 0.1975)
p = p + facet_wrap(~SITE_CODE.factor)

p

png('/home/holutz/projects/CAWS/wunifrac/wunif.bdisp.year_by_site.png', width=8, height=8, units='in', res=500)
plot(p)
dev.off()
```


<strong>png:</strong> 2



![png](output_21_1.png)



```R
library(viridis)

p = plot_ordination(ps.t, unif.dist.log, color = "year_coll.factor") + #, shape = "SampleType", label="X.SampleID") + 
    ggtitle("Unweighted UniFrac by SITE_CODE and year_collected") +
    geom_point(size=3, alpha=0.8) + 
    scale_colour_viridis(option="viridis", discrete=TRUE, direction=1) +
    #stat_ellipse(level=0.95, geom="polygon",alpha = .2, aes(fill = year_collected), linetype=0) +
    #scale_fill_viridis(option="viridis", discrete=FALSE, direction=1)+
    coord_fixed(sqrt(unif.dist.evals[2] / unif.dist.evals[1]))

p = p + guides(colour = guide_legend(title = "year_collected", ncol = 2, keywidth = 1, keyheight = 1))

p = p + theme_minimal() + 
        theme(panel.grid.major = element_blank(),
              axis.title.y=element_text(margin=margin(0,20,0,0)),
              text=element_text(size=8, color="black",family="Arial"),
              panel.grid.minor = element_blank(),
              #axis.ticks = element_blank(),
              legend.key.size = unit(.1, "in"),
              legend.spacing.x = unit(.01,"in"),
              legend.position = "bottom",
              plot.title = element_text(hjust = 0.5))

p$layers <- p$layers[-1]

p

p2 = p + facet_wrap(~SITE_CODE.factor)

p2

png('/home/holutz/projects/CAWS/unifrac/unif_by_year_v2.png', width=8, height=8, units='in', res=500)
plot(p)
dev.off()

png('/home/holutz/projects/CAWS/unifrac/unif_by_sitecode.facet_v2.png', width=8, height=8, units='in', res=500)
plot(p2)
dev.off()
```


![png](output_22_0.png)



<strong>png:</strong> 2



<strong>png:</strong> 2



![png](output_22_3.png)



```R
#Estimate betadispersion - eighted unifrac
beta.disp.uf = betadisper(unif.dist,group = sample_data(ps.t)$SITE_CODE.factor)
beta.disp.uf = beta.disp.uf$distances
beta.disp.uf = data.frame(d = beta.disp.uf, X.SampleID = names(beta.disp.uf))
beta.disp.uf2 = left_join(beta.disp.uf, sample_tab, by="X.SampleID")
```


```R
kruskal.test(d~year_coll.factor, data=beta.disp.uf2)

krusk.uf = dunn.test(beta.disp.uf2$d, beta.disp.uf2$year_coll.factor, method = "Bonferroni")
write.table(krusk.uf, "/home/holutz/projects/CAWS/unifrac/kw.uf.bdisp.txt", sep="\t")
```


    
    	Kruskal-Wallis rank sum test
    
    data:  d by year_coll.factor
    Kruskal-Wallis chi-squared = 47.936, df = 6, p-value = 1.217e-08



      Kruskal-Wallis rank sum test
    
    data: x and group
    Kruskal-Wallis chi-squared = 47.9359, df = 6, p-value = 0
    
    
                               Comparison of x by group                            
                                     (Bonferroni)                                  
    Col Mean-|
    Row Mean |     YR2013     YR2014     YR2015     YR2016     YR2017     YR2018
    ---------+------------------------------------------------------------------
      YR2014 |  -1.378688
             |     1.0000
             |
      YR2015 |   0.151936   1.902153
             |     1.0000     0.6001
             |
      YR2016 |  -2.577716  -1.569486  -3.326786
             |     0.1044     1.0000    0.0092*
             |
      YR2017 |  -3.540567  -2.784956  -4.530136  -1.127399
             |    0.0042*     0.0562    0.0001*     1.0000
             |
      YR2018 |  -2.609084  -1.649962  -3.271226  -0.202871   0.832401
             |     0.0953     1.0000    0.0112*     1.0000     1.0000
             |
      YR2019 |  -4.603495  -4.078960  -5.572175  -2.602677  -1.642558  -2.259155
             |    0.0000*    0.0005*    0.0000*     0.0971     1.0000     0.2507
    
    alpha = 0.05
    Reject Ho if p <= alpha/2



```R
#This code will arrange your values by the mean betadispersion so the boxplots look nice; 
#Need to be careful if color codes matter, may want to manually color if so
avg_obs <- aggregate(d ~ SITE_CODE.factor, data = beta.disp.uf2, mean)
beta.disp.uf2$SITE_CODE.factor <-factor(beta.disp.uf2$SITE_CODE.factor, levels=avg_obs[order(avg_obs$d), "SITE_CODE.factor"])
```


```R
p = ggplot(beta.disp.uf2, aes(x = year_coll.factor, y = d, fill = year_coll.factor)) + ggtitle("Unweighted Unifrac betadispersion")

p = p + geom_boxplot(outlier.stroke = 0.2, outlier.shape=20) + 
        ylab("Distance of samples from group centroid ") + 
        xlab("Year") #+ 
        #scale_y_continuous(limits=c(0,0.2))
        #coord_flip(ylim = c(0,1))

p = p + scale_fill_viridis(discrete=TRUE, alpha=0.8)

p = p + theme_bw() + theme(text=element_text(size=12, color="black",family = "Arial"),
                           legend.key = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           legend.position="none",
                           axis.text.x=element_text(angle=45,vjust=0.4))
                           #axis.title.x=element_blank(),
                           #axis.title.y=element_text(margin=margin(0,10,0,0)))

#p = p + geom_hline(yintercept = 0.1975)
p = p + facet_wrap(~SITE_CODE.factor)

p

png('/home/holutz/projects/CAWS/unifrac/unif.bdisp.year_by_site.png', width=8, height=8, units='in', res=500)
plot(p)
dev.off()
```


<strong>png:</strong> 2



![png](output_26_1.png)


#### Plotting fecal coliform data for Fig.X


```R
fc = read.csv("/home/holutz/projects/CAWS/FC_data/FC_Figure_Data.csv")
fc$date.collected = as.Date(fc$date.collected, format = "%m/%d/%Y")
```


```R
head(fc)
tail(fc)
```


<table>
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>CAWS.</th><th scope=col>LIMS.</th><th scope=col>Site</th><th scope=col>date.collected</th><th scope=col>FEC.COL_CFU</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;date&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>221 </td><td>6859117</td><td>OB</td><td>2013-05-14</td><td> 5500</td></tr>
	<tr><th scope=row>2</th><td>#N/A</td><td>6876034</td><td>OB</td><td>2013-06-04</td><td>24000</td></tr>
	<tr><th scope=row>3</th><td>182 </td><td>6880799</td><td>OB</td><td>2013-06-11</td><td> 9800</td></tr>
	<tr><th scope=row>4</th><td>209 </td><td>6886600</td><td>OB</td><td>2013-06-18</td><td> 8800</td></tr>
	<tr><th scope=row>5</th><td>222 </td><td>6891959</td><td>OB</td><td>2013-06-25</td><td>21000</td></tr>
	<tr><th scope=row>6</th><td>164 </td><td>6897267</td><td>OB</td><td>2013-07-02</td><td>11000</td></tr>
</tbody>
</table>




<table>
<caption>A data.frame: 6 × 5</caption>
<thead>
	<tr><th></th><th scope=col>CAWS.</th><th scope=col>LIMS.</th><th scope=col>Site</th><th scope=col>date.collected</th><th scope=col>FEC.COL_CFU</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;date&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>276</th><td>3698</td><td>8485065</td><td>CA</td><td>2019-05-28</td><td>40</td></tr>
	<tr><th scope=row>277</th><td>3717</td><td>8504745</td><td>CA</td><td>2019-06-24</td><td> 5</td></tr>
	<tr><th scope=row>278</th><td>3737</td><td>8525070</td><td>CA</td><td>2019-07-22</td><td>90</td></tr>
	<tr><th scope=row>279</th><td>3757</td><td>8550720</td><td>CA</td><td>2019-08-26</td><td>30</td></tr>
	<tr><th scope=row>280</th><td>3807</td><td>8597145</td><td>CA</td><td>2019-10-28</td><td>60</td></tr>
	<tr><th scope=row>281</th><td>3832</td><td>8616661</td><td>CA</td><td>2019-11-25</td><td> 5</td></tr>
</tbody>
</table>




```R
site.labs <- c("Calumet", "O'Brien")
names(site.labs) <- c("CA", "OB")
h <- 400

p = ggplot(fc, aes(x = date.collected, y = FEC.COL_CFU, color = FEC.COL_CFU)) + ggtitle("Fecal Coliform Counts by Site")
p = p + geom_point(size = 2) + 
        ylab("Coliform Counts (CTU/100mL)") +
        scale_y_continuous(trans='log10') +
        geom_hline(yintercept=400,linetype="dashed") +
        scale_color_gradient(trans="log10",low="blue",high="red")

p = p + theme_bw() + theme(panel.grid.major = element_blank(),
              axis.title.y=element_text(margin=margin(0,20,0,0)),
              text=element_text(size=12, color="black",family="Arial"),
              panel.grid.minor = element_blank(),
              #axis.text.x = element_blank(),
              #axis.ticks = element_blank(),
              axis.text.x=element_text(angle=90,vjust=0.4),
              #legend.key.size = unit(.25, "mm"),
              #legend.spacing.x = unit(.2,"in"),
              legend.position = "right",
              axis.title.x=element_blank())

p = p + scale_x_date(breaks = seq(as.Date("2013-05-14"), as.Date("2019-11-25"), by="1 year"),labels=date_format("%Y"))

p = p + facet_wrap(~Site, nrow=2,labeller = labeller(Site = site.labs))

p

png('/home/holutz/projects/CAWS/FC_data/FC_plot_gradient.png', width=12, height=8, units='in', res=300)
plot(p)
dev.off()
```

    Warning message:
    “Removed 2 rows containing missing values (geom_point).”
    Warning message:
    “Removed 2 rows containing missing values (geom_point).”



<strong>png:</strong> 2



![png](output_30_2.png)

