# The pQBR mercury resistance plasmids: a model set of sympatric environmental mobile genetic elements: raw data

Associated with manuscript: ***The pQBR mercury resistance plasmids: a model set of sympatric environmental mobile genetic elements***

Manuscript authors: Victoria T. Orr, Ellie Harrison, Rosanna Wright, Damian Rivett , James P. J. Hall

Analysis scripts are available at https://github.com/jpjh/PQBR_PLASMIDS

---
### Conj_Fit

Raw data and rmarkdown code for analyis of plasmid conjugation rates and transposon mobilisation rates, as well as comparative fitness experiment. 

- Hybrid_data.csv is the raw data files for mobilisation rate analyis and comparative fitness analysis. Dataset includes 3 experiments conducted and includes multiple donor strains (strains with plasmid and target chromosomal tranposon) created independently. Headers:
  		- Label (labelling of donor strains for lab use), 
  		- Plasmid (plasmid used in replicate, "Con" refers to control replicates, "_L" refers to in-house strain. "pQBR57_L", "pQBR103_L", "pQBR55_L" are the strain analyses in Hall et al. 2021. "pQBR57" is "pQBR57R".
  		- Plas_group - groups that plasmids belong to 1:	pQBR Group II	2: pQBR Group III	3: pQBR Group I	4:	5: IncP-2 (pP19E3.1)	6: pQBR Group IV
  		- Strain_rep - replicate number of that particular donor strain
  		- Plas_rep - replicate number of that particular plasmid
  		- Plas - labelling of donor strains for each plasmid
  		- Quantity_S - quantiy of culture plated at start of experiment for population counts
  		- Dilution_S - dilution of culture plated at start of experiment for population counts
  		- Quantity_E - quantiy of culture plated at end of experiment for population counts
  		- Dilution_E - dilution of culture plated at end of experiment for population counts
  		- Quantity_P - quantiy of culture plated at end of experiment for plasmid conjugation counts
  		- Dilution_P - dilution of culture plated at end of experiment for plasmid conjugation counts
  		- Quantity_T - quantiy of culture plated at end of experiment for transposon mobilisation counts
  		- Dilution_T - dilution of culture plated at end of experiment for transposon mobilisation counts
  		- Start_DC - population counts at start of experiment for donor strain
  		- Start_RC - population counts at start of experiment for recipient strain
  		- End_DC - population counts at end of experiment for donor strain
  		- End_RC - population counts at end of experiment for recipient strain
  		- End_NC - total population counts at end of experiment
  		- Plas_TC - transconjugant counts of plasmid selective plates, includes both "just plasmid" and "plasmid + transposon" transconjugants
  		- Tran_TC - transconjugant counts of transpson selective plates, includes only "plasmid + transposon" transconjugants
  		- Data - Date experimental data was collected
  		- Dat_rep - Counting combinations of date and donor e.g. plasmid1 (donor1:date1) = 1, plasmid1 (donor2:date1) = 2, plasmid1 (donor1:date2) = 3, plasmid1 (donor2:date2) = 4
- VTO_conjugation_comfit.Rmd is R markdown file containing code to analyse data and create plots from manuscript.

### Growth_curves

Raw and processed data, and R markdown for analysis from plate reader growth curves.

- rawdata - folder containing plate reader raw data
- tmp - folder containing information for procesing raw data
- VTO_growth_curves.Rmd - R markdown file containing analysis and code for production of manuscript plots

### Heatmaps&Trees

- GroupI: PIRATE genes family analysis results and Rmd for producing plots in manuscript for Group I plasmids
- GroupIII: PIRATE genes family analysis results and Rmd for producing plots in manuscript for Group III plasmids
- GroupIV: PIRATE genes family analysis results and Rmd for producing plots in manuscript for Group IV plasmids

### PAES_heatmap

Plasmids sequences (*.fna), mash distances and R markdown code to producing plots in manuscript.

### PlasmidSize_GC

- nuccore.csv - PLSDB download containing plasmid data
- pqbr2.csv - plasmid data from pQBR collection
- size_gc.Rmd - R markdown for analyis and producing manuscript figures. 

