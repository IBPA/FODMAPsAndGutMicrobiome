# Dietary FODMAPs, Gut Microbiome and Health in Irritable Bowel Syndrome
Investigating the interplay between dietary FODMAPs, gut microbiome and health in irritable bowel syndrome (IBS) using data from multiple studies.

## Code Architecture
```
├─ installation/
   ├─ install.sh
   ├─ install_packages.R
├─ preprocessing/
   ├─ s1/
   ├─ s2/
   ├─ ...
├─ integrated_analysis/
   ├─ Fig2.R
   ├─ Fig3.R
   ├─ Fig4.R
   ├─ Fig5.R
├─ run_all.sh
```

### Reproducability
* **Step 0:** Run ```install.sh``` to install required software packages.
* **Step 1:** Place the microbiome data for each study under the corresponding directory under (```preprocessing/study_id/microbiome_data/```). For some studies (s1, 23, s4), the raw microbiome data can be fetched using a corresponding ```fetch_data.sh``` script while for other studies (s2, s5, s6) they were directly requested from the article authors. See table below for the list of studies analyzed. See ```preprocessing/study_id/microbiome_data/format.txt``` for the expected format.
* **Step 2:** Place the metadata for each study in ```preprocessing/study_id/metadata.csv``` with format specified at ```metadata_format.csv```.
* **Step 3:** Run ```run_all.sh```. Note that the preprocessing step takes the longest and were performed on cluster node with 64GB RAM and two Intel E5-2630 v3 2.4GHz CPU’s each having 8 cores. The final figures would be stored under a ```results``` directory.
