#paper #results

## Infected mice undergo significant body weight loss

Mouse behavioral changes and body weight were evaluated to determine its health status. In particular, mouse body weight was measured every 24 hours for a period of seven days. Compared to naive mice, infected mice experienced a significant decrease in body weight soon after infection, which continued until day 7. Uninfected mice increased their body weight over time, which reflected their health status. The difference in body weight between uninfected and infected mice was assessed using the Mann-Whitney test. Significant differences were present from day 2 to day 7 after infection, indicating a long-lasting effect of the infection (Figure 2).

## _In vitro_ spleen stimulation with pneumococcal strain TIGR4 activates several genes related to different branches of the immune system

Transcriptomic data from spleens of infected and uninfected mice with or without homologous _in vitro _stimulation were analyzed. We performed an Independent Principal Component Analysis (IPCA) and its sparse version, sIPCA, both proposed by MixOmics package (Le Cao, 2016), (i) to observe the distribution of our data, (ii) to understand how stimulation at different time points affects the clustering of samples, (iii) to identify the genes responsible for the main variance among samples, and (iv) to find possible outlier samples.

The IPCA approach (Supplementary figure 1) yielded a better clusterization among experimental groups and time points compared to PCA (data not shown). An outlier control sample was detected and removed. The sparse version of the IPCA (sIPCA, Figure 3) applies soft-thresholding in the independent loading vectors in IPCA, performing feature selection. The graph shows the presence of two well-defined groups in the sIPC1: the stimulated and unstimulated samples. 

To better understand the genes that drive the formation of these clusters, the normalized expression values of the 50 genes selected by the first component of the sIPCA were divided into two heatmaps (Figure 3). Genes driving the unstimulated cluster included positive and negative regulators of the immune response, and they presented a decreased expression after stimulation The stimulation was mainly driven by genes related to cytokines, chemokines, and inflammation, all of them presenting an increased expression compared to unstimulated samples. 

## Stimulation of infected spleens highlights biological pathways of pneumococcus infection and new pathways related to stimulation

We then proceeded with the differential expression analysis using the _DESeq2 _package. To understand the biological alterations caused by the infection and the subsequent _in vitro _stimulation with inactivated pneumococcus, enrichment analysis was performed using three different comparisons. (i) Spleens from infected mice, (ii) stimulated spleens from non-infected mice, and (iii) stimulated spleens from infected mice, at different time points after infection, were all compared with control spleens. 

The number of differentially expressed genes for each condition at each time point is presented in Figure 4. As expected, the stimulation of infected samples led to a higher number of differentially expressed genes in comparison to only infected samples. 

The enrichment analysis was performed using the _Blood Transcription Modules (BTM)_ database and the _tmod_ package, and the results of the different groups are shown in Figure 5. In total, 87 modules were significantly enriched, only 3 of them being specifically activated in infected samples, while 40 modules were only activated after stimulation of previously infected samples. 

### Activation of extracellular matrix, cell adhesion, and Innate Immune response modules

Five modules were consistently activated at almost all time points in both unstimulated and _in vitro _stimulated groups. Related to monocytes, immune activation, TLR signaling, and cell cycle, these modules showed a different pattern after stimulation, presenting more down-regulated genes. Following the same direction, modules related to the extracellular region, monocytes, and cell cycle are especially enriched in down-regulated genes after stimulation (Figure 5). 

The “extracellular region cluster” module shows the downregulation of genes involved in the interaction with extracellular components, growth control, and the vascular endothelium/angiogenesis (HSPG2, GH1, ENG). Inside the same module, the monocyte chemoattractant CCL2 is also down-regulated, while CCL18, important for the recruitment of T lymphocytes but not monocytes, is up-regulated. 

_Moreover, the stimulation down-regulates genes responsible for the proliferation and differentiation of monocytes and macrophages like CSF2RA, CSF1R, and CSF3R, the latter one also important for adhesion. Inside the monocytes modules, other genes linked to the extracellular matrix and cell adhesion were also found down-regulated after stimulation (EMILIN2, CD36, CD93, CD302, HPSE, CLEC12A). _

_The downregulation of extracellular matrix and adhesion genes could be due to the process of in vitro stimulation, decreasing the cell adhesion to microtubes by shifting the cell priorities in response to the presence of a pathogen. However, the number of monocytes modules enriched after stimulation indicates that, indeed, the recruitment of monocytes seems negatively regulated._

### Activation of cell cycle, cytokines, and adaptive immune response modules

Unstimulated infected samples presented very few specific modules, such as the “recruitment of neutrophils” and “enriched in cell cycle”, while uninfected stimulated samples did not enrich any specific modules. On the other hand, the stimulation of infected samples led to the enrichment of many biological pathways not activated in the previous comparisons. Overall, the modules were related to antiviral response, antigen presentation, T cells,  B cells, and chemokines (Figure 5).  

On the first day, unstimulated samples presented the enrichment of T cell and cell cycle modules. After stimulation, these same modules are activated, together with many others related to T cells and cell cycle, in both cases enriched mainly by down-regulated genes. Despite the significant enrichment only at day 1, most of the activated genes in these modules are also present at other time points.

Among the down-regulated genes in T cell modules, there are cell-cycle genes and genes linked to cell adhesion, like VCAM1 and SIR3PG, while ITGA4, another adhesion-related gene, was up-regulated in unstimulated samples but presented no change after stimulation. Negative regulators of the T cell activity (LILRB4, LILRB3, SIT1) were also downregulated, while the few up-regulated genes were mainly related to T cell activation (CD3E, GRAP2, CDCA7, and LAT). 

On the other hand, the specific modules in late time points were mostly activated by up-regulated genes. We observe a stronger activation of the “cytokines − receptors cluster” module and the specific enrichment of pathways like leukocyte differentiation, signaling in T cells, enriched in B cells, among other modules.	

Cytokine modules are activated from the stimulation of baseline samples and looking inside these modules, indeed we see many up-regulated genes independently of the time point, especially those from the CCL family (CCL4, CCL18, CCL3, CCL3L3, and CCL3L1), IL1A, IL1RN, and TNF. On the other hand, many cytokine genes (CSF2, IL2RA, IL6, IL10, IL1B) present a slight up-regulation after stimulation of baseline samples, but an increase in the fold change in the late time points after infection.

Despite the high number of activated modules, the response to the stimuli after a previous infection does not show general up-regulation of the immune and inflammatory response. The second contact with the pathogen through the in vitro stimulation permitted us to appreciate biological processes which could not be detected in the primary infection, especially those related to antigen presentation, adaptive immunity, and cytokines. These processes are possibly related to a recall immune response starting within the first days after infection. 

## Cytokines Assay suggests a recall immune response composed of both innate and adaptive branches of the immune system
**Regarding the concentration of cytokines in spleen culture supernatants, the infection caused punctual changes. The main increases in the concentration occurred on days 4 and 7 after infection and were much more expressive in stimulated samples (Figure or Supplementary material XX). **

**At baseline, the stimulation process induced a significant increase in four cytokines concentrations, two of which presented a considerable difference between stimulated and control samples: KC and MIP-1a (Median of differences of 42,46 and 184,4 respectively). Despite cytokine changes between stimulated and control samples being noticeable in early stages, they increased considerably upon in vitro stimulation, at days 4 and 7 after the infection (Figure 6)**

In samples from infected animals, some cytokines presented a significant decrease on day 2 after infection, including GM-CSF, IL-10, IL12p40, IL12p70, MCP-1, MIP-1b, Rantes, TNFa, and IL-2, while a significant increase in concentration occurs only on day 7 for IL-17a in these samples. 

Comparison of stimulated samples from days 4 and 7 after infection with only infected samples showed a significant increase in almost all cytokines concentrations (Figure 6 and Suppl. XX), indicating the involvement of both innate and adaptive branches of the immune system. This increase was more accentuated on day 7, in which the difference in the median between the groups was 2597 pg/ml for IL-17a, 769 for GM-CSF, 374 for IFNgamma, and 182 for G-CSF.

	

## Gene expression and cytokines data integration indicates a recall immune response after stimulation, especially from day 4 after infection

To identify the genes correlated with the increase in the concentration of cytokines, especially at late time points, data integration was performed using the sparse version of Partial Least Squares (sPLS), from MixOmics package. As expected, the _in vitro _stimulated samples formed a different cluster compared to the non-stimulated samples, although there is a different behavior regarding time points in each cluster (Figure 10A).

In the non-stimulated cluster, we see a perturbation caused by infection, but some samples from day 7 cluster together with control samples from day 0, a pattern probably related to the resolution of the infection. On the other hand, _in vitro _stimulated samples presented a different pattern, samples from days 4 and 7 form a new cluster, driven by the increase in cytokine concentration and the expression of certain genes, including many cytokine genes like Il2, Il2ra, Csf2, and Il10 as observed in the correlation circle plot (Figure 10B).

The highest values of correlation were found for Cd69, Csf2, Il2ra, and Il10 (Figure 3C). Among the genes positively correlated with cytokine concentrations, three of them are members of the TNF Receptor Superfamily (Tnfrsf9, Tnfrsf4 and Tnfrsf4), two of them encode for protooncogene that acts as a serine/threonine protein kinase (Pim1 and Pim2) and one is a microRNA, Mir155hg. 

 some of them are negative regulators of the immune system, like Cd300a, Sirpa, and Cd5l, which was already described as a repressor of Th17 cell pathogenicity (Wang, C. et al., 2015). 

## Biomarkers of _in vitro _stimulation

We aimed to understand if feature selection could summarize the impact of a previous infection on stimulated samples, indicating possible biomarkers of this infection.** **We applied the DaMiR-seq package, which provides data normalization, feature selection, and classification, based on different machine learning techniques. Three groups were established based on the transcriptomics and cytokine data distribution, focusing on the stimulation of uninfected samples, samples from early time points after infection (days 1 and 2), and samples from late time points (days 4 and 7).

Eleven genes were chosen by applying a threshold of 0.5 to the scaled importance score identified by the DaMiR-seq package (supplementary materials). These 11 genes allowed a clear clusterization of the three groups (Figure 7). Three of these genes are involved in innate immunity and were selected by the package as representants of the stimulation of infected samples: Fpr1, Nlrp3, and Slpi. They presented an increased expression, independently of the time point, when compared with the stimulation of uninfected samples.

The stimulation of early time points led to the increase in the expression of other inflammatory genes like Serpinb2 and Chil1, compared to the stimulated control group, and these values started to decrease in the subsequent days. At late time points, three other important genes, related to cytokine activity, had their expression increased when compared to the other groups: Ccr4, Csf2, and Il2.

Despite the small number of samples being a limitation for this type of analysis, the feature selection summarizes the new immunological processes that arise after the stimulation of infected samples and suggests the use of the _in vitro _stimulation model to confirm the presence of a previous infection by measuring the expression of a few genes. 