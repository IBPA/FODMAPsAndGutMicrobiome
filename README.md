# Dietary FODMAPs, Gut Microbiome and Health in Irritable Bowel Syndrome
Investigating the interplay between dietary FODMAPs, gut microbiome and health in irritable bowel syndrome (IBS) using data from multiple studies.

## Reproducability
* **Step 0:** Run ```install_packages.R``` to install required packages.
* **Step 1:** Place the microbiome data for each study under the corresponding directory under (```preprocessing/study_id/microbiome_data/```). For some studies (s1, s3, s4), the raw microbiome data can be fetched using a corresponding ```fetch_data.sh``` script while for other studies (s2, s5, s6) they were directly requested from the article authors. See table below for the list of studies analyzed. See ```preprocessing/study_id/microbiome_data/format.txt``` for the expected format.
* **Step 2:** Place the metadata for each study in ```preprocessing/study_id/metadata.csv``` with format specified at ```metadata_format.csv```.
* **Step 3:** Run ```run_all.sh``` and see the final figures under the ```results``` directory. Note that the preprocessing step was performed on cluster node with 64GB RAM and two Intel E5-2630 v3 2.4GHz CPU’s each having 8 cores and may take up to a day for a given study.

## Requirements
* [QIIME 2](https://docs.qiime2.org/2020.2/install/) (version ≥ 2018.11)
* [PICRUSt](https://picrust.github.io/picrust/install.html) (version 1.1.4)

## Code Architecture
```
├─ fig/
   ├─ Fig1.png
   ├─ Fig5A.png
├─ integrated_analysis/
   ├─ generate_Fig2.R
   ├─ generate_Fig3.R
   ├─ ...
├─ lib/
   ├─ abundance_management.R
   ├─ classification.R
   ├─ ...
├─ preprocessing/
   ├─ s1/
   ├─ s2/
   ├─ ...
├─ install_packages.sh
├─ run_all.sh
├─ settings.txt
```
