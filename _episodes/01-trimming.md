---
title: "Trimming and Filtering"
teaching: 30
exercises: 25
questions:
- "How can I get rid of sequence data that doesn't meet my quality standards?"
objectives:
- "Clean FASTQ reads using cutadapt."
- "Select and set multiple options for command-line bioinformatic tools."
- "Write `for` loops with two variables."
keypoints:
- "The options you set for the command-line tools you use are important!"
- "Data cleaning is an essential step in a genomics workflow."
---

# Cleaning Reads

In the previous episode, we took a high-level look at the quality
of each of our samples using FastQC. We visualized per-base quality
graphs showing the distribution of read quality at each base across
all reads in a sample and extracted information about which samples
fail which quality checks. Some of our samples failed quite a few quality metrics used by FastQC. This doesn't mean,
though, that our samples should be thrown out! It's very common to have some quality metrics fail, and this may or may not be a problem for your downstream application. For our variant calling workflow, we will be removing some of the low quality sequences to reduce our false positive rate due to sequencing error.

We will use a program called
[cutadapt](https://github.com/marcelm/cutadapt)
to filter poor quality reads and trim poor quality bases from our samples.

## Cutadapt Options

Cutadapt has a variety of options to trim your reads. If we run the following command, we can see some of our options.

~~~
$ cutadapt
~~~
{: .bash}

Which will give you the following output:

~~~
This is cutadapt 2.10 with Python 3.7.3                   
Command line parameters:                                  
Run "cutadapt --help" to see command-line options.        
See https://cutadapt.readthedocs.io/ for full documentation.

cutadapt: error: You did not provide any input file names. Please give me something to do!
~~~
{: .output}

In order to see the parameter options we need to pass -h as follows:

~~~
$ cutadapt -h
~~~
{: .bash}

There are many parameters here so we will just show the top of the output which explains the usage, but you should look through and familiarise yourself with the options available. Importantly it shows you what are the required parameters and also the version of the software you have used which is important to keep note of. You should always record what versions of software you have used and what parameters you ran the software with in order to make your analysis reproducible

~~~
cutadapt version 2.10

Copyright (C) 2010-2020 Marcel Martin <marcel.martin@scilifelab.se>

cutadapt removes adapter sequences from high-throughput sequencing reads.

Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
~~~
{: .output}

This output shows us that there is a different usage for single end or paired end reads. For single ended data we need to use the (`-a`) flag to specify the sequence of the adaptor to be removed, (`-o`) to specify the name of the trimmed output, and then the name of the raw sequence we want to be trimmed after. For paired end data, cutadapt expects two adaptor sequences for trimming, specified by (`-a`) and (`-A`). The output again is specified by (`-o`) for one trimmed read file, and (`-p`) for the other trimmed read file. The two files to be trimmed are then given as the last arguements.

The first input file typically contains a `_1` or `_R1` in the name
The second input file typically contains a `_2` or `_R2` in the name

For more information about the arguments used see the [cutadapt documentation](https://cutadapt.readthedocs.io/en/stable/guide.html)

An example command for cutadapt will look like the command below. This command is an example and will not work, as we do not have the files it refers to:

~~~
$ cutadapt -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
           -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
           -o trimmed.R1.fastq.gz -p trimmed.R2.fastq.gz \
           reads.R1.fastq.gz reads.R2.fastq.gz \
           > sample_fastq.gz.log &
~~~
{: .bash}


In this example, we've told cutadapt:

| code   | meaning |
| ------- | ---------- |
|`-q 20`| Bases need to have a quality phred score of at least 20 to be maintained |
| `-a` | The sequence to trim from the forward reads, this is based on illuminas truseq adaptors|
| `-A` | The sequence to trim from the reverse reads, this is also based on illumina truseq adaptors |
| `-o trimmed.R1.fastq.gz` | the name of the output of the trimmed forward reads |
| `-p trimmed.R2.fastq.gz` | the name of the output of the trimmed reverse reads |
| `reads.R1.fastq.gz` | the input forward reads file |
| `reads.R2.fastq.gz` | the input reverse reads file |
| ` > sample_fastq.gz.log` | this writes the output log of cutadapt to a file |


> ## Multi-line commands
> Some of the commands we ran in this lesson are long! When typing a long
> command into your terminal, you can use the `\` character
> to separate code chunks onto separate lines. This can make your code more readable.
{: .callout}



## Running cutadapt

Now we will run cutadapt on our data. To begin, navigate to your `untrimmed_fastq` data directory:

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq
~~~
{: .bash}

We are going to run cutadapt on one of our paired-end samples.
While using FastQC we saw that Nextera adapters were present in our samples.
The adaptor sequences for both reads according to the illumina website is CTGTCTCTTATACACATCT

We will  remove bases if their phred score is below 20 (like in our example above). The ampersand (&) tells linux to run the process in the background and the greater than symbol (>) pipes the output of the command to a file

~~~
$ cutadapt -q 20 -a CTGTCTCTTATACACATCT \
           -A CTGTCTCTTATACACATCT \
           -o SRR2589044_1.trim.fastq.gz -p SRR2589044_2.trim.fastq.gz \
            SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz \
            > SRR2589044_fastq.gz.log &

~~~
{: .bash}

~~~
cutadapt output log to be put here:

~~~
{: .output}

> ## Exercise
>
> Use the output from your cutadapt command to answer the
> following questions.
>
> 1) What percent of reads did we discard from our sample?
> 2) What percent of reads did we keep both pairs?
>
>> ## Solution
>> 1) 0.23%
>> 2) 79.96%
> {: .solution}
{: .challenge}


We can confirm that we have our output files:

~~~
$ ls SRR2589044*
~~~
{: .bash}

~~~
SRR2589044_1.fastq.gz       SRR2589044_2.trim.fastq.gz
SRR2589044_1.trim.fastq.gz  SRR2589044_2.fastq.gz
~~~
{: .output}

The output files are also FASTQ files. It should be smaller than our
input file, because we've removed reads. We can confirm this:

~~~
$ ls SRR2589044* -l -h
~~~
{: .bash}

~~~
Need to rerun with cutadapt to update
-rw-rw-r-- 1 dcuser dcuser 124M Jul  6 20:22 SRR2589044_1.fastq.gz
-rw-rw-r-- 1 dcuser dcuser  94M Jul  6 22:33 SRR2589044_1.trim.fastq.gz
-rw-rw-r-- 1 dcuser dcuser 128M Jul  6 20:24 SRR2589044_2.fastq.gz
-rw-rw-r-- 1 dcuser dcuser  91M Jul  6 22:33 SRR2589044_2.trim.fastq.gz
~~~
{: .output}


We've just successfully run cutadapt on one of our FASTQ files!
However, there is some bad news. cutadapt can only operate on
one sample at a time and we have more than one sample. The good news
is that we can use a `for` loop to iterate through our sample files
quickly!

We unzipped one of our files before to work with it, let's compress it again before we run our for loop.

~~~
gzip SRR2584863_1.fastq
~~~
{: .bash}

~~~
$ for infile in *_1.fastq.gz
> do
>   base=$(basename ${infile} _1.fastq.gz)
>   cutadapt -q 20 -a CTGTCTCTTATACACATCT \
>           -A CTGTCTCTTATACACATCT \
>           -o ${base}_1.trim.fastq.gz -p ${base}_2.trim.fastq.gz \
>            ${base}_1.fastq.gz ${base}_2.fastq.gz \
>            > v${base}_fastq.gz.log &
> done

~~~
{: .bash}


Go ahead and run the for loop. It should take a few minutes for
cutadapt to run for each of our six input files. Once it's done
running, take a look at your directory contents. You'll notice that even though we ran cutadapt on file `SRR2589044` before running the for loop, there is only one set of files for it. Because we matched the ending `_1.fastq.gz`, we re-ran cutadapt on this file, overwriting our first results. That's ok, but it's good to be aware that it happened.

~~~
$ ls
~~~
{: .bash}

~~~
SRR2584866_1.fastq.gz         SRR2589044_1.trim.fastq.gz
SRR2584863_1.fastq.gz         SRR2584866_1.trim.fastq.gz
SRR2584863_1.trim.fastq.gz    SRR2589044_2.fastq.gz
SRR2584866_2.fastq.gz         SRR2589044_2.trim.fastq.gz
SRR2584863_2.fastq.gz         SRR2584866_2.trim.fastq.gz
SRR2584863_2.trim.fastq.gz    SRR2589044_1.fastq.gz
~~~
{: .output}



We've now completed the trimming and filtering steps of our quality
control process! Before we move on, let's move our trimmed FASTQ files
to a new subdirectory within our `data/` directory.

~~~
$ cd ~/dc_workshop/data/untrimmed_fastq
$ mkdir ../trimmed_fastq
$ mv *.trim* ../trimmed_fastq
$ cd ../trimmed_fastq
$ ls
~~~
{: .bash}

~~~
SRR2584863_1.trim.fastq.gz    SRR2584866_1.trim.fastq.gz    SRR2589044_1.trim.fastq.gz
SRR2584863_2.trim.fastq.gz    SRR2584866_2.trim.fastq.gz    SRR2589044_2.trim.fastq.gz
~~~
{: .output}

> ## Bonus Exercise (Advanced)
>
> Now that our samples have gone through quality control, they should perform
> better on the quality tests run by FastQC. Go ahead and re-run
> FastQC on your trimmed FASTQ files and visualize the HTML files
> to see whether your per base sequence quality is higher after
> trimming.
>
>> ## Solution
>>
>> In your AWS terminal window do:
>>
>> ~~~
>> $ fastqc ~/dc_workshop/data/trimmed_fastq/*.fastq*
>> ~~~
>> {: .bash}
>>
>> In a new tab in your terminal do:
>>
>> ~~~
>> $ mkdir ~/Desktop/fastqc_html/trimmed
>> $ scp dcuser@ec2-34-203-203-131.compute-1.amazonaws.com:~/dc_workshop/data/trimmed_fastq/*.html ~/Desktop/fastqc_html/trimmed
>> ~~~
>> {: .bash}
>>
>> Then take a look at the html files in your browser.
>>
>> Remember to replace everything between the `@` and `:` in your scp
>> command with your AWS instance number.
>>
>> After trimming and filtering, our overall quality is much higher,
>> we have a distribution of sequence lengths, and more samples pass
>> adapter content. However, quality trimming is not perfect, and some
>> programs are better at removing some sequences than others. Because our
>> sequences still contain 3' adapters, it could be important to explore
>> other trimming tools like [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to remove
>> these, depending on your downstream application. Cutadapt did pretty well though, and its performance
>> is good enough for our workflow.
> {: .solution}
{: .challenge}
