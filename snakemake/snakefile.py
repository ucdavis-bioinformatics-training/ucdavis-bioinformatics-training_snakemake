"""
UC Davis Bioinformatics Core
Snakemake pipeline for RNA-seq read quantification using HTS PREPROCESSING and STAR
Author: Keith Mitchell (kgmitchell@ucdavis.edu)
"""

#TODO
#  fix path for the slurm out and the samples file
#  clean up the overalapping files (in STAR for example all output are the same files...) get from function
#  -F option for htspreproc, more parameters for this?
#  stranded option
#  add overlapper for pe
import json
import sys
import glob

args = {}
sbatch_args = {}

###########################################################################
# PROCESS PARAMETERS
###########################################################################

# register the path of the --configfile
args['configfile'] = sys.argv[sys.argv.index('--configfile') + 1]

# check if file is passed or if just string is passed to get samples.
if 'samples' not in config.keys():
    sys.stderr.write("Getting samples from file %s" % config["__default__"]['samples_file'] + '\n')
    SAMPLES = []
    try:
        with open(config["__default__"]['samples_file'], 'r') as samples_file:
            for line in samples_file:
                if line.strip('\n') != '':
                    SAMPLES.append(line.strip('\n'))
    except FileNotFoundError:
        sys.exit("Could not find or properly parse the 'sample_file' that was passed.")

else:
    try:
        SAMPLES = str(config['samples']).split(',')
    except ValueError:
        sys.exit("SYSTEM EXIT: Could not parse the 'samples' arguments that was passed.")

# check what type is being run (PE, tagseq, SE)
if config["__default__"]['type'] in ['tagseq', 'PE', 'SE']:
    args['type'] = config["__default__"]["type"]
else:
    sys.exit("Configuration error: type not found. Please select from ['tagseq', 'PE', 'SE']. Default = PE")

sys.stderr.write("SAMPLES: %s" % ', '.join(SAMPLES) + '\n')

###########################################################################
# CHECK IF RUNNING LOCALLY OR ON SLURM CLUSTER
###########################################################################
# once running as an sbatch can check "if SLURM_TASK_PID in os.environ.keys()"
if config["__default__"]["running_locally"]=="True":
    args["running_locally"] = True
    sys.stderr.write("Running Locally" + '\n')

else:
    args["running_locally"] = False
    sys.stderr.write("Running on Cluster" + '\n')
    import socket
    sys.stderr.write("SOCKET: %s" % socket.gethostname() + '\n')

###########################################################################
# CORE SETUP
###########################################################################
args['basename'] = config['project']['basename']
args['id'] = config['project']['id']
args['fastqs'] = args['basename'] + '/' + config['project']['fastqs']
args['htsout'] = args['basename'] + '/' + '01-HTS_Preproc'
args['starout'] = args['basename'] + '/' + '02-STAR_alignment'
args['multiqc_out'] = args['basename'] + '/' + '02-HTS_multiqc_report'
args['countsout'] = args['basename'] + '/' + '03-Counts'
args['human_rrna_ref'] = config['project']['human_rrna_ref']
args['star_ref'] = config['project']['star_ref']

###########################################################################
# MODULES
###########################################################################

shell.prefix("set -o pipefail; ")
shell.prefix("module load star; module load htstream/1.3.2; module load multiqc/htstream.dev0;")
shell("module list")

sys.stderr.write(json.dumps(args, indent=1) + '\n')
sys.stderr.flush()

###########################################################################
# MASTER RULE
###########################################################################
rule master_rule:
    run:
        for sample in SAMPLES:
            mystdout = config['hts_star']['output'].replace('.out', args['type'] + '_' + sample + '.out')
            mystderr = config['hts_star']['error'].replace('.err', args['type'] + '_' + sample + '.err')
            if args['running_locally'] == False:
                command = f"sbatch " \
                       f"--job-name={config['hts_star']['job-name'] + sample} " \
                       f"--ntasks={config['hts_star']['ntasks']} " \
                       f"--nodes={config['hts_star']['n']} " \
                       f"--partition={config['hts_star']['partition']} " \
                       f"--time={config['hts_star']['time']} " \
                       f"--mem={config['hts_star']['mem']} " \
                       f"--output={mystdout} " \ 
                       f"--error={mystderr} " \
                       f"--mail-type={config['hts_star']['mail-type']} " \
                       f"--mail-user={config['hts_star']['mail-user']} " \
                       f"--wrap='snakemake -s snakefile.py all --configfile {args['configfile']} --config samples={sample}'"
            else:
                command = f"snakemake -s snakefile.py all --configfile {args['configfile']} --config samples={sample}"
            sys.stderr.write(command + '\n')
            shell(command)


###########################################################################
# HTSPREPROCESSING
###########################################################################


rule htspreproc:
    input:
        # TODO make functions for error handling if files found are not whats expected?
        read1 = lambda wildcards: glob.glob('/{path}/{samp}/{samp}*R1*'.format(samp=wildcards.sample, path=args['fastqs']))
    output:
        hts_log = '%s/{sample}/{sample}_htsStats.log' % args['htsout'],
        fastq1 = '%s/{sample}/{sample}_R1.fastq.gz' % args['htsout']
    params:
        ref = args['human_rrna_ref'],
        read2 = lambda wildcards: glob.glob('/{path}/{samp}/{samp}*R2*'.format(samp=wildcards.sample, path=args['fastqs']))
    run:
        make_hts_out = "[[ -d %s ]] || mkdir %s" % (args['htsout'], args['htsout'])
        shell(make_hts_out)
        make_sample_out = "[[ -d %s/{wildcards.sample} ]] || mkdir %s/{wildcards.sample}" % (args['htsout'], args['htsout'])
        shell(make_sample_out)
        stats_se = "hts_Stats -L {output.hts_log} -U {input.read1} -N 'initial stats'"
        stats_pe = "hts_Stats -L {output.hts_log} -1 {input.read1} -2 {params.read2} -N 'initial stats'"
        seq_screen1 = "hts_SeqScreener -A {output.hts_log} -N 'screen phix'"
        seq_screen2 = "hts_SeqScreener -s {params.ref} -r -A {output.hts_log} -N 'count the number of rRNA reads'"
        adapter_trim = "hts_AdapterTrimmer -A {output.hts_log} -N 'trim adapters'"
        # TODO what is this doing
        cut_trim1 = "hts_CutTrim -n -a 20 -A {output.hts_log}" # TODO this is not present for SE/PE
        qwindow_trim = "hts_QWindowTrim -A {output.hts_log} -N 'quality trim the ends of reads'"
        ntrimmer = "hts_NTrimmer -A {output.hts_log} -N 'remove any remaining N characters'"
        length_filter = "hts_LengthFilter -n -m 50 -A {output.hts_log} -N 'remove reads < 50bp'" #
        length_filter_se = "hts_LengthFilter -n -m 45 -A {output.hts_log} -N 'remove reads < 45p'"
        stats2 = "hts_Stats -A {output.hts_log} -f %s/{wildcards.sample}/{wildcards.sample} -N 'final stats'" % (args['htsout'])
        super_deduper = "hts_SuperDeduper -A {output.hts_log} -N 'remove PCR duplicates'"
        if args['type'] == 'tagseq':
            sys.stderr.write("HTSPREPROC: TagSeq" + '\n')
            master_list = [stats_se, seq_screen1, seq_screen2, adapter_trim, cut_trim1, qwindow_trim, ntrimmer,
                           length_filter, stats2]
            command = ' | '.join(master_list)
            sys.stderr.write(command + '\n')
            shell(command)
        elif args['type'] == 'SE':
            sys.stderr.write("HTSPREPROC: SE" + '\n')
            master_list = [stats_se, seq_screen1, seq_screen2, adapter_trim, qwindow_trim, ntrimmer, length_filter_se, stats2]
            command = ' | '.join(master_list)
            sys.stderr.write(command + '\n')
            shell(command)
        elif args['type'] == 'PE':
            sys.stderr.write("HTSPREPROC: PE" + '\n')
            master_list = [stats_pe, seq_screen1, seq_screen2, super_deduper, adapter_trim, qwindow_trim, ntrimmer,
                           length_filter, stats2]
            command = ' | '.join(master_list)
            sys.stderr.write(command + '\n')
            shell(command)
        else:
            sys.exit("CONFIGURATION ERROR (STAR): type not found. Please select from ['tagseq', 'PE', 'SE']. Default = PE")



###########################################################################
# STAR ALIGNMENT
###########################################################################

# TAGSEQ AND SE STAR

rule star:
    input:
        log = '%s/{sample}/{sample}_htsStats.log' % args['htsout'],
        fastq1 = '%s/{sample}/{sample}_R1.fastq.gz' % args['htsout']
    output:
        counts = '%s/{sample}/{sample}_ReadsPerGene.out.tab' % args['starout'],
        bam = '%s/{sample}/{sample}_Aligned.sortedByCoord.out.bam' % args['starout'],
        log = '%s/{sample}/{sample}_Log.out' % args['starout'],
        prog = '%s/{sample}/{sample}_Log.progress.out' % args['starout']
    params:
        temp = '%s/{sample}/{sample}__STARtmp' % args['starout'],
        ref = args['star_ref'],
        stderr = '%s/{sample}/{sample}-STAR.stderr' % args['starout'],
        stdout = '%s/{sample}/{sample}-STAR.stdout' % args['starout'],
        fastq2 = '%s/{sample}/{sample}_R2.fastq.gz' % args['htsout']
    run:
        if args['type'] == 'tagseq' or args['type'] == 'SE':
            sys.stderr.write("STAR: SE/TagSeq" + '\n')
            command = f"STAR --runThreadN 8 --genomeDir {params.ref} " \
                             f"--outSAMtype BAM SortedByCoordinate " \
                             f"--readFilesCommand zcat " \
                             f"--readFilesIn {input.fastq1} " \
                             f"--quantMode GeneCounts " \
                             f"--outFileNamePrefix {args['starout']}/{wildcards.sample}/{wildcards.sample}_ " \
                             f"> {params.stdout} 2> {params.stderr} "
            sys.stderr.write(command + '\n')
            shell(command)
        elif args['type'] == 'PE':
            sys.stderr.write("STAR: PE" + '\n')
            command = f"STAR --runThreadN 8 --genomeDir {params.ref} " \
                      f"--outSAMtype BAM SortedByCoordinate " \
                      f"--readFilesCommand zcat " \
                      f"--readFilesIn {input.fastq1} {params.fastq2} " \
                      f"--quantMode GeneCounts " \
                      f"--outFileNamePrefix  {args['starout']}/{wildcards.sample}/{wildcards.sample}_ " \
                      f"> {params.stdout} 2> {params.stderr}"
            sys.stderr.write(command + '\n')
            shell(command)
        else:
            sys.exit("CONFIGURATION ERROR (STAR): type not found. Please select from ['tagseq', 'PE', 'SE']. Default = PE")

###########################################################################
# GENERAL RULES
###########################################################################
rule hts_proc_stats:
    input:
        expand('%s/{sample}/{sample}_htsStats.log' % args['htsout'], sample=SAMPLES)
    output:
        summary_file = '%s/summary_hts_%s.txt' % (args['basename'], args['type'])
    run:
        sample_string = ','.join(SAMPLES)
        command = f"python summarize_stats.py --samples {sample_string} --samples_is_file False " \
                  f"--type {args['type']} --input_dir {args['basename']} --output_dir {args['basename']}"
        sys.stderr.write(command + '\n')
        shell(command)

rule counts_table:
    input:
        expand('%s/{sample}/{sample}_ReadsPerGene.out.tab' % args['starout'], sample=SAMPLES)
    output:
        final_counts = '%s/rnaseq_workshop_counts.txt' % args['countsout']
    params:
        geneid_file = "%s/tmp/geneids.txt" % args['countsout'],
        temp_out = "%s/tmp/tmp.out" % args['countsout'],
        samp_file = '%s/samples.txt' % args['countsout']
    run:
        # TODO clean this up a bit more? add stranded option
        shell("mkdir %s/tmp" % args['countsout'])
        for sample in SAMPLES:
            shell(f"echo {sample} ")
            cmd = f"cat {args['starout']}/{sample}/{sample}_ReadsPerGene.out.tab | tail -n +5 | cut -f4 > " \
                  f"{args['countsout']}/tmp/{sample}.count "
            sys.stderr.write(cmd + '\n')
            shell(cmd)
        # can use any sample for this command (here we will use the one at the end of the for loop)
        gt_geneids = f"tail -n +5 {args['starout']}/{sample}/{sample}_ReadsPerGene.out.tab | cut -f1 > {params.geneid_file} "
        sys.stderr.write(gt_geneids + '\n')
        shell(gt_geneids)

        # finally combine all of the files using the `paste` command
        paste_cmd = f"paste {params.geneid_file} {args['countsout']}/tmp/*.count > {params.temp_out} "
        sys.stderr.write(paste_cmd + '\n')
        shell(paste_cmd)

        # create a row in the file of sample names
        with open(params.samp_file, 'w') as temp_file:
            for sample in SAMPLES:
                temp_file.write(sample + '\n')
        row_cmd = f"cat <(cat {params.samp_file} | sort | paste -s) {params.temp_out} > {output.final_counts} "
        sys.stderr.write(row_cmd + '\n')
        shell(row_cmd)
        shell(f"rm -rf {args['countsout']}/tmp ")

rule multiqc:
    input:
        log = '%s/{sample}/{sample}_htsStats.log' % args['htsout'],
        fastq1 = '%s/{sample}/{sample}_R1.fastq.gz' % args['htsout']
    run:
        shell("mkdir %s" % args['multiqc_out'])
        # todo is this the name of the report??
        shell("multiqc -i HTSMultiQC-cleaning-report -o %s %s" %(args['multiqc_out'], args['htsout']))

rule all:
    input:
        expand('%s/{sample}/{sample}_Aligned.sortedByCoord.out.bam'  % args['starout'], sample=SAMPLES),
        expand('%s/{sample}/{sample}_Log.out'  % args['starout'], sample=SAMPLES)

rule full_analysis:
    input:
        hts_summary = '%s/summary_hts_%s.txt' % (args['basename'], args['type']),
        final_counts = '%s/rnaseq_workshop_counts.txt' % args['countsout']
        #TODO what is this output?
        # multiqc_out = ""
