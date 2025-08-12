# Visualise CNAs from EPI2ME wf-human-variation

Rshiny app to visualise CNAs detected with the [EPI2ME wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) workflow. 
Two plots are included:
- `Multi-chromosome` plot with a small overview of all chromosomes - coverage and CNAs
- `Single-chromosome` plot with coverage and CNAs for the selected chromosome  

These plots show the mean coverage per genomic region (blue line) and the CNA predictions from wf-human-variation in green (gain) or red (loss). The small triangles show the start position of the CNA (sometimes the CNAs are so small that the corresponding line is not clear). In addition, the chromosome cytoband is included below the single chromosome plot. The plots are interactive.

![App overview](images/app.png)

## Installation

This app was tested with R v4.5.0
If setting up with Rstudio: File --> New project --> Version Control --> Repository URL: `https://github.com/CarolinaPB/visualise_cnv_wf-human-variation.git` and choose where you want to save it

The R Renv package is used to take care of package dependencies. Use `renv::restore()` to install the necessary packages.

### Prepare data

The app expects a directory that contains the output files from the EPI2ME wf-human-variation: `<SAMPLE>.wf_cnv.vcf.gz` and `<SAMPLE>.regions.bed.gz`.   
This directory can contain the data for as many samples as you want. 
By default, the app expects the data to be in the included `data` directory. The data directory can be configured in the config file or it can be specified in the app

#### Set up config

The config file `cnv_nanopore/configs/config.yaml` is used to define the path for the data directory and other files.   
By default, the app will use the `data` directory in the current project directory. You can leave it as is or add your own path. Inside the app you'll also be able to select a different directory.

```yaml
default:
    data_path: "data"
    cytoband_file: "resources/cytoBandIdeo.txt"
    fai_file: "resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai"
```

## Launch app

You can start the app by running

```R
shiny::runApp('cnv_nanopore')
```
