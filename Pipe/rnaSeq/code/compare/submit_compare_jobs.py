#!/usr/bin/python3
## file: submit_compare_jobs.py
import os
import csv, subprocess

parameter_file_full_path = "rnaSeq/code/compare/compare_job_params.csv"

with open(parameter_file_full_path, mode='r', newline='', encoding='utf-8-sig') as csvfile:
    reader = csv.reader(csvfile)
    for job in reader:


        submit_command = ("sbatch "
            "--job-name=compare_{1}_{0} "
            "--output=rnaSeq/output/compare/{1}/{0}/{1}-{0}.out "
            "--export=DATASET={0},COMPAREMETHOD={1} rnaSeq/code/compare/compare_job.sh").format(*job)

        # print(submit_command) # Uncomment this line when done testing to use the submit command created
        # uncomment the following 3 lines when done testing to submit the jobs
        exit_status = subprocess.call(submit_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job compare_{0} failed to submit".format(submit_command))


print("Done submitting jobs!")

