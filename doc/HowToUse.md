## Steps to Reproduce Findings

1. Download all files from the following links to the data/data_raw directory:
- https://www.dropbox.com/scl/fo/3filfi60qhxgokmjb6weq/AGLqwdegrXF5NUYqAl7EI6Q/QTL_results/eQTL_results/eQTL_finemapping_results?dl=0&rlkey=ihnv8hpbj1ilb96cxsm8fc42p&subfolder_nav_tracking=1
- https://data.mendeley.com/datasets/y7hyt83vxr/1

2. You should have two .zip folders in your data_raw directory.

3. In terminal, navigate to analysis and run the following command:
- `./processArchives.sh`

4. In terminal, run the following command:
- `./processSprimes.sh`

5. In terminal, run the following command:
- `./alignSNPs.py`
- alignSNPs takes the following parameters: desired PIP search value, desired match type, and desired search file, **in that order.**
    - PIP search value defaults to 0.5, and accepts any number
    - Match type defaults to Neanderthals and accepts any of the following in a case-insensitive manner:
        - Neanderthal match type: n, nmatch, neanderthal
        - Denisovan match type: d, dmatch, denisovan
        - Both match type: b, bmatch, both
        - Either match type: e, ematch, either
    - Search file defaults to eQTL_finemapping.significantAssociations.MAGE.v1.0.txt and accepts either of the following in a case-insensitive manner:
        - SNPs with Significant Associations: sig
        - SNPs with Any Association: all, any
- alignSNPs prints logs of each run to log.txt, stored in the analysis folder.
- alignSNPs outputs alignments stamped with search start date, start time, and parameters in a tab separated values format in the data/data_clean/Alignments directory.