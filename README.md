# CNV Nanopore

Rshiny app to visualise CNVs detected with the [EPI2ME wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) workflow.



## Installation

This app was tested with R v4.5.0

```bash
git clone https://github.com/CarolinaPB/visualise_cnv_wf-human-variation.git
cd visualise_cnv_wf-human-variation
```

The R Renv package is used to take care of package dependencies. Use `renv::restore()` to set it up and install the necessary packages.

### Prepare data

The app expects a directory that contains the output files from the EPI2ME wf-human-variation: `<SAMPLE>.wf_cnv.vcf.gz` and `<SAMPLE>.regions.bed.gz`. This directory contains the data for as many samples as you want, all under the same directory. By default, the app expects the data to be in the included `data` directory.

### Set up config

The config file `cnv_nanopore/configs/config.yaml` is used to define the path where the data is. By default, the app will use the `data` directory in the current project directory. You can leave it as is or add your own path. Inside the app you'll also be able to select a different directory.

```yaml
default:
    data_path: "data"
    
testing:
    data_path: "/path/to/data"
```

## Launch app

You can start the app by running

```R
runApp('cnv_nanopore')
```