---
title: "Automating a Variant Calling Workflow"
teaching: 30
exercises: 15
questions:
- "How can I make my workflow more efficient and less error-prone?"
objectives:
- "Write a shell script with multiple variables."
- "Incorporate a `for` loop into a shell script."
keypoints:
- "We can combine multiple commands into a shell script to automate a workflow."
- "Use `echo` statements within your scripts to get an automated progress update."
---
# What is a shell script?

You wrote a simple shell script in a [previous lesson](http://www.datacarpentry.org/shell-genomics/05-writing-scripts/) that we used to extract bad reads from our
FASTQ files and put them into a new file.

Here's the script you wrote:

~~~
grep -B1 -A2 NNNNNNNNNN *.fastq > scripted_bad_reads.txt

echo "Script finished!"
~~~
{: .bash}

That script was only two lines long, but shell scripts can be much more complicated
than that and can be used to perform a large number of operations on one or many
files. This saves you the effort of having to type each of those commands over for
each of your data files and makes your work less error-prone and more reproducible.
For example, the variant calling workflow we just carried out had about eight steps
where we had to type a command into our terminal. Most of these commands were pretty
long. If we wanted to do this for all six of our data files, that would be forty-eight
steps. If we had 50 samples (a more realistic number), it would be 400 steps! You can
see why we want to automate this.

We've also used `for` loops in previous lessons to iterate one or two commands over multiple input files.
In these `for` loops, the filename was defined as a variable in the `for` statement, which enabled you to run the loop on multiple files. We will be using variable assignments like this in our new shell scripts.

Here's the `for` loop you wrote for unzipping `.zip` files:

~~~
$ for filename in *.zip
> do
> unzip $filename
> done
~~~
{: .bash}

And here's the one you wrote for running cutadapt on all of our `.fastq` sample files:


~~~
$ for infile in SRR2584863 SRR2584866 SRR2589044
> do
> cutadapt -q 20 -a CTGTCTCTTATACACATCT \
>           -A CTGTCTCTTATACACATCT \
>           -o $infile\_1.trim.fastq.gz -p $infile\_2.trim.fastq.gz \
>            $infile\_1.fastq.gz $infile\_2.fastq.gz \
>            >> $infile\_fastq.gz.log
> done
~~~

{: .bash}

Notice that in this `for` loop, we used one variable, `infile`, which was defined in the `for` statement.

> ## Creating Variables
> Within the Bash shell you can create variables at any time (as we did
> above, and during the 'for' loop lesson). Assign any name and the
> value using the assignment operator: '='. You can check the current
> definition of your variable by typing into your script: echo $variable_name.
{: .callout}

In this lesson, we'll use two shell scripts to automate the variant calling analysis: one for FastQC analysis (including creating our summary file), and a second for the remaining variant calling. To write a script to run our FastQC analysis, we'll take each of the commands we entered to run FastQC and process the output files and put them into a single file with a `.sh` extension. The `.sh` is not essential, but serves as a reminder to ourselves and to the computer that this is a shell script.

# analysing Quality with FastQC

We will use the command `touch` to create a new file where we will write our shell script. We will create this script in a new
directory called `scripts/`. Previously, we used
`nano` to create and open a new file. The command `touch` allows us to create a new file without opening that file.

~~~
$ mkdir -p ~/cs_course/scripts
$ cd ~/cs_course/scripts
$ touch read_qc.sh
$ ls
~~~
{: .bash}

~~~
read_qc.sh
~~~
{: .output}

We now have an empty file called `read_qc.sh` in our `scripts/` directory. We will now open this file in `nano` and start
building our script.

~~~
$ nano read_qc.sh
~~~
{: .bash}

**Enter the following pieces of code into your shell script (not into your terminal prompt).**

Our first line will ensure that our script will exit if an error occurs, and is a good idea to include at the beginning of your scripts. The second line will move us into the `untrimmed_fastq/` directory when we run our script.

~~~
set -e
cd ~/cs_course/data/untrimmed_fastq/
~~~
{: .output}

These next two lines will give us a status message to tell us that we are currently running FastQC, then will run FastQC
on all of the files in our current directory with a `.fastq` extension.

~~~
echo "Running FastQC ..."
fastqc *.fastq*
~~~
{: .output}

Our next line will create a new directory to hold our FastQC output files. Here we are using the `-p` option for `mkdir` again. It is a good idea to use this option in your shell scripts to avoid running into errors if you don't have the directory structure you think you do.

~~~
mkdir -p ~/cs_course/results/fastqc_untrimmed_reads
~~~
{: .output}

Our next three lines first give us a status message to tell us we are saving the results from FastQC, then moves all of the files
with a `.zip` or a `.html` extension to the directory we just created for storing our FastQC results.

~~~
echo "Saving FastQC results..."
mv *.zip ~/cs_course/results/fastqc_untrimmed_reads/
mv *.html ~/cs_course/results/fastqc_untrimmed_reads/
~~~
{: .output}

The next line moves us to the results directory where we've stored our output.

~~~
cd ~/cs_course/results/fastqc_untrimmed_reads/
~~~
{: .output}

The next five lines should look very familiar. First we give ourselves a status message to tell us that we're unzipping our ZIP
files. Then we run our for loop to unzip all of the `.zip` files in this directory.

~~~
echo "Unzipping..."
for filename in *.zip
    do
    unzip $filename
    done
~~~
{: .output}

Next we concatenate all of our summary files into a single output file, with a status message to remind ourselves that this is
what we're doing.

~~~
echo "Saving summary..."
cat */summary.txt > ~/cs_course/docs/fastqc_summaries.txt
~~~
{: .output}

> ## Using `echo` statements
>
> We've used `echo` statements to add progress statements to our script. Our script will print these statements
> as it is running and therefore we will be able to see how far our script has progressed.
>
{: .callout}

Your full shell script should now look like this:

~~~
cd ~/cs_course/data/untrimmed_fastq/

echo "Running FastQC ..."
fastqc *.fastq*

mkdir -p ~/cs_course/results/fastqc_untrimmed_reads

echo "Saving FastQC results..."
mv *.zip ~/cs_course/results/fastqc_untrimmed_reads/
mv *.html ~/cs_course/results/fastqc_untrimmed_reads/

cd ~/cs_course/results/fastqc_untrimmed_reads/

echo "Unzipping..."
for filename in *.zip
    do
    unzip $filename
    done

echo "Saving summary..."
cat */summary.txt > ~/cs_course/docs/fastqc_summaries.txt
~~~
{: .output}

Save your file and exit `nano`. We can now run our script. Use pwd to check which directory you are in, if you are not in the scripts directory cd into the correct place as shown below, before running the script:

~~~
$ cd ~/cs_course/scripts/
$ bash read_qc.sh
~~~
{: .bash}

~~~
Running FastQC ...
Started analysis of SRR2584866.fastq
Approx 5% complete for SRR2584866.fastq
Approx 10% complete for SRR2584866.fastq
Approx 15% complete for SRR2584866.fastq
Approx 20% complete for SRR2584866.fastq
Approx 25% complete for SRR2584866.fastq
.
.
.
~~~
{: .output}

For each of your sample files, FastQC will ask if you want to replace the existing version with a new version. This is
because we have already run FastQC on this samples files and generated all of the outputs. We are now doing this again using
our scripts. Go ahead and select `A` each time this message appears. It will appear once per sample file (six times total).

~~~
replace SRR2584866_fastqc/Icons/fastqc_icon.png? [y]es, [n]o, [A]ll, [N]one, [r]ename:
~~~
{: .output}


# Automating the Rest of our Variant Calling Workflow

We can extend these principles to the entire variant calling workflow. To do this, we will take all of the individual commands that we wrote before, put them into a single file, add variables so that the script knows to iterate through our input files and write to the appropriate output files. This is very similar to what we did with our `read_qc.sh` script, but will be a bit more complex.

Our variant calling workflow has the following steps:

1. Index the reference genome for use by bwa and samtools.
2. Align reads to reference genome.
3. Convert the format of the alignment to sorted BAM, with some intermediate steps.
4. Calculate the read coverage of positions in the genome.
5. Detect the single nucleotide polymorphisms (SNPs).
6. Filter and report the SNP variants in VCF (variant calling format).

Let's go through this script together. Lets first cd into the scripts folder and open nano with the following command:



~~~
$ cd ~/cs_course/scripts
$ nano run_variant_calling.sh
~~~
{: .bash}

The script should look like this. You should copy and paste the contents of the output box below and save the script. We can then go through the script line by line:

~~~
cd ~/cs_course/results

bwa index ../data/ecoli_rel606.fasta

for file in SRR2584863 SRR2584866 SRR2589044
do
	echo "working with file $file"
	bwa mem ../data/ecoli_rel606.fasta ../data/trimmed_fastq_small/$file\_1.trim.sub.fastq ../data/trimmed_fastq_small/$file\_2.trim.sub.fastq > $file.aligned.sam
	samtools view -S -b $file.sam > $file.aligned.bam
	samtools sort -o $file.aligned.sorted.bam $file.aligned.bam
	samtools index $file.aligned.sorted.bam
	bcftools mpileup -O b -o $file\_raw.bcf -f ../data/ecoli_rel606.fasta $file.aligned.sorted.bam
	bcftools call --ploidy 1 -m -v -o $file\_variants.vcf $file\_raw.bcf
	vcfutils.pl varFilter $file\_variants.vcf > $file\_final_variants.vcf

done

~~~
{: .output}

Now, we'll go through each line in the script before running it.

First, notice that we change our working directory so that we can create new results subdirectories
in the right location.

~~~
cd ~/cs_course/results
~~~
{: .output}

Next we index our reference genome for BWA:

~~~
bwa index ../data/ecoli_rel606.fasta
~~~
{: .output}

Then, we use a loop to run the variant calling workflow on each of our FASTQ files. The full list of commands
within the loop will be executed once for each of the FASTQ files in the
`data/trimmed_fastq_small/` directory.
We will include a few `echo` statements to give us status updates on our progress.

The first thing we do is assign the name of the FASTQ file we're currently working with to a variable called `$file` and
tell the script to `echo` the filename back to us so we can check which file we're on.
In this variable $file we are just giving a list of the sample name without the suffix

~~~
for file in SRR2584863 SRR2584866 SRR2589044
do
	echo "working with file $file"
~~~
{: .bash}

We are using the base of this name and in order to access both the _1 and _2 input files we use
$file and the rest of the file name _2.trim.sub.fastq for the _2 input file. The first time through this loop the computer will interpret '$file\_2.trim.sub.fastq' as SRR2584863_2.trim.sub.fastq. The '\' is a special character which allows us to add the variable name $file to the string to _2.trim.sub.fastq. These files are one folder up in the data folder, and then within a folder within this, called trimmed_fastq_small.

~~~
bwa mem ../data/ecoli_rel606.fasta ../data/trimmed_fastq_small/$file\_1.trim.sub.fastq ../data/trimmed_fastq_small/$file\_2.trim.sub.fastq > $file.aligned.sam
~~~
{: .bash}

And finally, the actual workflow steps:

1) align the reads to the reference genome and output a `.sam` file:

~~~
bwa mem ../data/ecoli_rel606.fasta ../data/trimmed_fastq_small/$file\_1.trim.sub.fastq ../data/trimmed_fastq_small/$file\_2.trim.sub.fastq > $file.aligned.sam
~~~
{: .output}

2) convert the SAM file to BAM format:

~~~
    samtools view -S -b $file.sam > $file.aligned.bam
~~~
{: .output}

3) sort the BAM file:

~~~
    samtools sort -o $file.aligned.sorted.bam $file.aligned.bam
~~~
{: .output}

4) index the BAM file for display purposes:

~~~
    samtools index $file.aligned.sorted.bam
~~~
{: .output}

5) calculate the read coverage of positions in the genome:

~~~
    bcftools mpileup -O b -o $file\_raw.bcf -f ../data/ecoli_rel606.fasta $file.aligned.sorted.bam
~~~
{: .output}

6) call SNPs with bcftools:

~~~
    bcftools call --ploidy 1 -m -v -o $file\_variants.vcf $file\_raw.bcf
~~~
{: .output}

7) filter and report the SNP variants in variant calling format (VCF):

~~~
    vcfutils.pl varFilter $file\_variants.vcf > $file\_final_variants.vcf
~~~
{: .output}



> ## Commenting your code
> It's a good idea to add comments to your code so that you (or a collaborator) can make sense of what you did later.
> Look through your existing script and add comments to describe what it does. Anything following
> a `#` character will be interpreted as a comment - the terminal will not try to run these comments as code.
{: .callout}


Now we can run our script:

~~~
$ bash run_variant_calling.sh
~~~
{: .bash}


> ## Exercise
>
> The samples we just performed variant calling on are part of the long-term evolution experiment introduced at the
> beginning of our variant calling workflow. From the metadata table, we know that SRR2589044 was from generation 5000,
> SRR2584863 was from generation 15000, and SRR2584866 was from generation 50000. How did the number of mutations per sample change
> over time? Examine the metadata table. What is one reason the number of mutations may have changed the way they did?
>
> If you are attending an instructor-led workshop, discuss these questions in your breakout rooms. Nominate someone to put a summary of your thoughts onto the Padlet.
>
> Hint: You can find a copy of the output files for the subsampled trimmed FASTQ file variant calling in the
> `~/.solutions/wrangling-solutions/variant_calling_auto/` directory.
>
>> ## Solution
>>
>> ~~~
>> $ for infile in ~/cs_course/results/vcf/*_final_variants.vcf
>> > do
>> >     echo ${infile}
>> >     grep -v "#" ${infile} | wc -l
>> > done
>> ~~~
>> {: .bash}
>>
>> For SRR2589044 from generation 5000 there were 10 mutations, for SRR2584863 from generation 15000 there were 25 mutations,
>> and SRR2584866 from generation 766 mutations. In the last generation, a hypermutable phenotype had evolved, causing this
>> strain to have more mutations.
> {: .solution}
{: .challenge}


> ## Bonus Exercise
>
> If you have time after completing the previous exercise, use `run_variant_calling.sh` to run the variant calling pipeline
> on the full-sized trimmed FASTQ files. You should have a copy of these already in `~/cs_course/data/trimmed_fastq`, but if
> you don't, there is a copy in `~/.solutions/wrangling-solutions/trimmed_fastq`. Does the number of variants change per sample?
{: .challenge}
