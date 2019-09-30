# Snakemake Pipeline for TagSeq and RnaSeq

1. Alter the `templates/tagseq.json`, `templates/pe.json`, `templates/se.json` file.
    - be sure to specify the variable `running_locally` as True or False.
    - specify any other paths as necessary. 
    
2. Activate the proper version of snakemake:
    - `module load snakemake`
    - `source activate snakemake`    

3. Now you would likely want to run a set of samples (A) or an individual sample for a certain type of dataset (PE RNAseq,
    TAGseq, or SE RNAseq) so also be sure to specify the `type` as `PE`, `tagseq`, or `SE` in the respective json files. 
    - **A)** Run a set of samples:
        + Via list: `snakemake -s snakefile.py master_rule --configfile templates/pe.json --config samples=SampleAC1,SampleAC2,SampleAC3`
        + Via file (specified in .json): `snakemake -s snakefile.py master_rule --configfile templates/pe.json`
    - **B)** Run an individual sample:
        + `snakemake -s snakefile.py master_rule --configfile templates/pe.json --config samples=SampleAC1`
        
4. Finally for the samples that you have finished you will likely want to run stats on the hts processing.
    - **A)** Run a set of samples:
        + Via list: `snakemake hts_proc_stats --config sample=SampleAC1,SampleAC2,SampleAC3 type=PE`
        + Via file (master=True): `snakemake hts_proc_stats --config sample=samples.txt master=True type=PE`
    - **B)** Run an individual sample:
        + `snakemake hts_proc_stats --config sample=SampleAC1 type=PE`


A few notes about the pipeline: 
  - The `master_rule` looks at the variable `running_locally` and submits the rule `all` as either an sbatch or
   a local call for the variable `False` and `True` respectively.
  - The rule `all`, which is simply called by the rule `master_rule` looks at if hts_preproc or star have been completed
  (in that order).
  - If `running_locally` is `False` be sure to edit the file `templates/tagseq.json`, `templates/pe.json`, 
  `templates/se.json` for the specific sbatch parameters for the job. 

A few notes about snakemake:
  - If a job fails then the files marked as `output` for a rule will be removed. 
 
 
 #TODO 
 fix the glob function
 test the hts stats and star counts command
 get rid of the ID parameter 
 