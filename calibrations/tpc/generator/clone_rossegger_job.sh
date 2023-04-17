#!/bin/bash

# This shell script follows the pattern of the valid sample job, with some hints
# about what parts are what, and generates a job by bXing number for all other
# source files in the provided filelist.
# If you are changing the format of the filename, you may need to alter the code
# so that it continues to extract and replace the correct items.
#
# in: (hardcoded) example file, hints about what to replace, and filelist
# out: set of job_NNNNN.job files, each of which can be submitted to condor.

#specify what parts of the sample to replace:
sample_job=condor_rossegger_test.job
sample_job_signature=26 36 40 101 2

#specify what to replace them with (NNNNN means 'number of sets', IIIII means 'this set id')
this_job_signature=26 36 40 NNNNN IIIII

njobs=100

#loop over the elements in that list

for ((i=0; i<njobs; i++))
do
    
    #find the root number of the crossing
    crossing=`echo $sourcefile | grep -Eo '[0-9]{6,9}' `
    #make the name of the specific job
    jobname=jobs/job_ross_${i}.job
    cp $sample_job $jobname
    echo assigning job: $jobname with id: $i 
    echo $sample_job_signature becoming $this_job_signature
    #filter the instances of the original filename with the new filename in the specific job
    sed -i "s|${sample_job_signature}|${this_job_signature}|g" $jobname
    sed -i "s|NNNNN|${njobs}|g" $jobname
    sed -i "s|IIIII|${i}|g" $jobname

    #uncomment to run just one, for testing:
    #exit
done
