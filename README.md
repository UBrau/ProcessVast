# ProcessVast
Scripts and files to perform routine analyses of alternative splicing, based on [_vast-tools_](https://github.com/vastgroup/vast-tools) output 

## Summary
This is a tool to make sense of _vast-tools_ output by providing several standard analyses, including similarily of single samples, numbers of changing events between specified contrasts, correlation between multiple contrasts, and files of genes containing changing AS events for GO analysis, including specific backgrounds.

## Dependencies
R packages _optparse_ and _gplots_. An attempt is made to install them.

## Usage
`ProcessVast.R -h`

## Input
### Required
Necessary to run are, a completed vast-tools analysis, including `diff` for every contrast that will be compared.
- **INCLUSION... table** (primary vast-tools output)
- **vast-tools _diff_ results tables**. It is recommended to run `diff` with default statistical options except increase --size to 2000.
- **Sample table**: CSV file with (at least) columns _Sample_ and _Type_. _Sample_ must be identical to what is specified in the INCLUSION... table; all samples that are considered replicates have the same _Type_. See examples.
- **Contrast table**: CSV file with (at least) columns _Experimental_, _Control_, and _File_. The former two contain entries from the sample table's _Type_ column that specify which sample types are compared (by subtraction of PSI/PIR). _File_ species the path and file name of the `vast-tools diff` output table for each contrast. See examples.
### Optional
- **Specific event table**: Tab-deliminted file containing information about known regulation of events. If present, a bar graph will be produced showing the overlap of the changing events in each contrast with known differential events. E.g. neural-differential splicing would be encoded as the entries 'neural-LOW' or 'neural-HIGH' while non-differential events are NA. Columns are _EVENT_ (matching vast-tools events) and any number of columns with descritive names (e.g., 'neural', 'ASD') and events containing xxx-HIGH/xxx-LOW. The provided files should be considered examples and may not be accurate or comprehensive.
- **Lookup table for EVENT to GeneID**: If files for GO analysis should be produced, a tab-delimited file containing columns _EVENT_ and_EnsemblGeneID_ is required. Files for Mm2, Hsa and Hs2 are available in the 'EVENTtoGeneID' folder.

## Functionality
_ProcessVast.R_ will first perform filtering of raw PSI/PIR based on quality columns. Averages for each sample type are taken from `diff` outputs and averaged if a sample occurs in multiple ones. If more than half of the replicates in a type are NA for an event after filtering, the mean is set to NA.

For Alt5 and Alt3 events that represent alternative variants, a single score is calculated analogously to PSI, where -100 indicates a complete shift from the most distal to the most proximal splice set relative to the exon, and 100 the reverse.

Differential events are those that have have the specified absolute dPSI/dPIR (5, 10, or 15) and a significant change (>0) according to`vast-tools diff` 

## Output
- **Filtered INCLUSION... table** (INCLUSION_LEVELS_FULL-\*_clean<date>.tab.gz)
- **Changing events tables** (AS.DiffEvents_dPSI.\*.tab.gz): Contain differential events from all contrasts.
- **Sample type means** (ASmeans.tab.gz)
- **R archive of relevant tables** (AStables.Robj), including event info, means, dpsi, and diff tables at 5/10/15 cutoffs. Load into R using *load()*
- **Tarball with GO analysis input files** (GO.tar.gz; optional) and event type-specific backgrounds, unorderdered, at dPSI/dPIR cutoff of 10.
  
- **Bar graphs of changing events per type** (ChangingEvents.dPSI.\*_bars.pdf) 
- **Multi-dimensional scaling of single samples** (MDS_samplePSI.pdf) to assess replicate similarity
- **Clustered correlation heatmap of single samples** (PSI.cluster.corr.single.pdf) to assess replicate similarity
- **Clustered correlation heatmap of dPSI/dPIR** (dPSI.cluster.corr.pdf)
- **Scatterplots of dPSI** (dPSI.scatter.pdf) for up to nine contrasts against one another (or more, if --scatterForce is specified)
- **Clustered heatmap of dPSI/dPIR** (AS.cluster.dPSI.pdf) for all events. The legend for event types is provided in the file AS.cluster.legend.dPSI.pdf.
- **Clustered heatmap of dPSI/dPIR per event type** (AS.cluster_types_dPSI10.pdf)
- **Bar graphs of overlap of changes with known regulation** (SpecificEvents.dPSI.10.pdf; optional)

- **Log file** (vastResultsProcessing.log)

## Author
Please let me know if you have questions or encounter errors by emailing or raising an issue.

Ulrich Braunschweig, University of Toronto

[email](mailto:u.braunschweig@utoronto.ca)

## Reference
If using in published work, please cite the DOI of the release. The latest one is:
[![DOI](https://zenodo.org/badge/267636603.svg)](https://zenodo.org/badge/latestdoi/267636603)
  

