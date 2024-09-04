## SyRI_IS 

Code created to detect SVs with multiple IS-related events.
SV detection by SyRI did not run well without masking nor with softmasking so hardmasking of IS loci was performed before running SyRI.

Run SyRI with hard masked genomes: `runSyRI_IS.py`.
Annotate IS insertion sites: `classify_IS_events.py`.
See run examples in `ipynb/*_multiple_runs.ipynb`.