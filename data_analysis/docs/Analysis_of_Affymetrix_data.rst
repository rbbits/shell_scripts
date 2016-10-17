Introduction
~~~~~~~~~~~~

This page describes the steps to analyse a data set of Affymetrix microarray data.

Acute Myeloid Leukemia (AML) gene expression data from the Cancer Genome Atlas (TCGA) performed using the Affymetrix Human Genome U133 Plus 2.0 Array is used to describe this method, but the script can be used with other data sources as long as they are Affymetrix CEL files.

Data Input
~~~~~~~~~~

Only samples molecularly classified as “Normal Karyotype” were used for this analysis comprising 200 .CEL files. These files were provided by the Haematological Cancer Genetics Team (Vassiliou Faculty - 163) and are held in their team's NFS space.

For this analysis, the files are saved in the ``input`` directory.

Experimental Design
~~~~~~~~~~~~~~~~~~~

The groups used in the experimental design are described within tab-separated-values (TSV) targets file.

Targets file specifications
^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Notice how the numbers in ``arrayID`` don't have to be continuous

-  In the example, the name of the grouping column is ``Sample`` but it can have any name and then passed to the script using the ``--group_column`` option.

-  The group name must be a single word:

   -  **Don't leave spaces in the name** or other special characters (e.g. ?, /, +, #, !, <TAB>, <SPACE>, etc.).
   -  Use underscores ( _ ) to join separate words if necessary.
   -  Use only to separate columns, although other separators can be used.

 | ``arrayID  filename                          sampleID  bioRep  Sample``
 | ``9        TCGA-AB-2807-03A-01R-0757-21.CEL  s9        1       setbphi_2_Runx1mut``
 | ``30       TCGA-AB-2833-03A-01R-0757-21.CEL  s30       1       setbp1hi_2_RunxWT``
 | ``37       TCGA-AB-2843-03A-01R-0757-21.CEL  s37       1       setbp1hi_2_RunxWT``
 | ``56       TCGA-AB-2865-03A-01R-0757-21.CEL  s56       1       setbphi_2_Runx1mut``
 | ``76       TCGA-AB-2885-03A-01R-0757-21.CEL  s76       1       setbp1hi_2_RunxWT``
 | ``81       TCGA-AB-2890-03A-01R-0757-21.CEL  s81       1       setbphi_2_Runx1mut``
 | ``90       TCGA-AB-2899-03A-01R-0757-21.CEL  s90       1       setbphi_2_Runx1mut``
 | ``96       TCGA-AB-2908-03A-01R-0757-21.CEL  s96       1       setbp1hi_2_RunxWT``
 | ``105      TCGA-AB-2917-03A-01R-0757-21.CEL  s105      1       setbp1hi_2_RunxWT``
 | ``110      TCGA-AB-2922-03A-01R-0757-21.CEL  s110      1       setbp1hi_2_RunxWT``
 | ``115      TCGA-AB-2927-03A-01R-0757-21.CEL  s115      1       setbphi_2_Runx1mut``
 | ``120      TCGA-AB-2933-03A-01R-0757-21.CEL  s120      1       setbphi_2_Runx1mut``
 | ``122      TCGA-AB-2935-03A-01R-0757-21.CEL  s122      1       setbp1hi_2_RunxWT``
 | ``123      TCGA-AB-2936-03A-01R-0757-21.CEL  s123      1       setbphi_2_Runx1mut``
 | ``130      TCGA-AB-2943-03A-01R-0757-21.CEL  s130      1       setbp1hi_2_RunxWT``
 | ``136      TCGA-AB-2949-03B-01R-0757-21.CEL  s136      1       setbphi_2_Runx1mut``
 | ``143      TCGA-AB-2959-03A-01R-0757-21.CEL  s143      1       setbphi_2_Runx1mut``
 | ``173      TCGA-AB-2995-03A-01R-0757-21.CEL  s173      1       setbp1hi_2_RunxWT``
 | ``180      TCGA-AB-3002-03A-01R-0757-21.CEL  s180      1       setbp1hi_2_RunxWT``

Analysis Materials
~~~~~~~~~~~~~~~~~~

The script requires the following elements:

-  A **targets file**, usually stired in a directory called ``resources`` but this is not a requirement as long as the correct path (full or relative is passed to the ``--targets=`` argument.

-  A **contrast**. The contrast is a representation of the differential expression analysis that the script is going to perform in the form of a simple subtraction ``A-B`` where A is usually the group of mutant samples and B the wild type.

-  Directory structure of the working directory: It is recommended that the analysis is performed inside a working directory with a number of predefined directories, but others can be used as long as these are appropriately specified in the arguments.

   -  **input**: contains the CEL files
   -  **resources**: contains the targets file
   -  **output**: where the results are saved
   -  **bin**: contains the R script

The R script options
^^^^^^^^^^^^^^^^^^^^

The script makes an assumption of the required files are stored in those directories, **BUT** different parameter values can be passed to it via the following command line arguments:

-  ``--input_dir=/path/to/input/files/``
-  ``--output_dir=/path/where/results/are/saved/``

The paths can be absolute or relative e.g.: ``--input_dir=./a/dir/within/``

R Libraries
^^^^^^^^^^^

Currently the script uses R version 3.2.2 and the path to the necessary libraries are hard-coded. If you wish to use a different version or a different location for your libraries simply change the code to reflect this. Remember the version of R used to run Rscript should match the one that was used to compile and install the libraries.

Analysis
~~~~~~~~

To perform a differential expression analysis using the targets file and a specific contrast based on the group of interest contained in it, in this case column ``Sample`` use the R script provided with the correct arguments:

 ``/software/R-3.2.2/bin/Rscript run_affymetrix.R --targets_file=./resources/Targets_1.2.txt --contrast=setbphi_2_Runx1mut-setbp1hi_2_RunxWT --group_column=Sample --output_dir=./output --input_dir=./input``

Analysis using LSF
^^^^^^^^^^^^^^^^^^

For bigger analyses - e.g. involving more than a dozen of samples - memory and other resources may become an issue and prevent the analysis to be completed. In those cases it's advisable to use the Platform Load Sharing Facility (LSF) which will schedule the job and manage the resources accordingly based on the user's requests.

To run the ``Rscript`` command using LSF wrap it around with a ``bsub`` command with the appropriate request for resources, enclosing it within single or double quotes:

**Run in a single line**: delete backslashes (\\):

 | ``bsub \``
 | ``-J ``\ “``affyrun``”\ `` \``
 | ``-oo contrast_1.2.out \``
 | ``-eo contrast_1.2.err \``
 | ``-R ``\ \ ``'select[mem>8000] rusage[mem=8000]'``\ \ `` \``
 | ``-M 8000 \``
 | \ ``'/software/R-3.2.2/bin/Rscript ./bin/olly_dovey_200AMLsAffy.R --targets_file=./resources/Targets_1.2.txt --contrast=setbphi_2_Runx1mut-setbp1hi_2_RunxWT --group_column=Sample  --output_dir=./output --input_dir=./input'``\ 

Notice the values for ``bsub``'s options ``-oo`` and ``-eo`` (``contrast_1.2.out`` and ``contrast_1.2.err``) are the names for the output and error files respectively which ``bsub`` uses to log its output and errors that might have happened when trying to schedule/run the job. After setting the job off ``bsub`` will display a message similar to this:

 ``Job <XXXXXX> is submitted to default queue <normal>.``

Where *XXXXXX* is the **job ID** assigned by LSF. Use this number to monitor the job's status - possible values are: ``RUN``, ``PEND``, ``EXIT``, ``DONE``:

 ``bjobs XXXXXX``

Each of them mean:

``PEND``: the job is pending scheduling. If too long in this status (more than 24hr) seek advice.

``RUN``: the job has been scheduled successfully and it's running.

``EXIT``: the job has failed, seek help.

``DONE``: the job has been completed successfully, check your results.

Any other status: seek advice.

Results
~~~~~~~

Data has been sorted by significance (column adj.P.Val). The lower the adj.P.Val, the more significant the change in gene expression. Adj.P.Val and logFC (log Fold-Change: the amount by which gene expression changes) are strongly correlated, but not absolutely correlated. Two adj.P.Val cut-offs can be used in microarrays, (i) signiffcant genes are those where adj.P.Val is less than or equal to 0.01; (ii) signiffcant genes are those where adj.P.Val is less than or equal to 0.05. A positive logFC implies the gene is upregulated in the first of the conditions (usually ascribed to the mutant). Only one output file in plain text format is produced and saved in the ``output`` directory.

Quality Control
~~~~~~~~~~~~~~~

No QC performed as the files are of external origin but a snipet of the code to do this can be found commented in the script.
