#!/usr/bin/env python

# Copyright (C) 2022- The University of Notre Dame
# This software is distributed under the GNU General Public License.
# See the file COPYING for details.

# This example shows some of the data handling features of taskvine.
# It performs a BLAST search of the "Landmark" model organism database.
# It works by constructing tasks that download the blast executable
# and landmark database from NCBI, and then performs a short query.

# Each task in the workflow performs a query of the database using
# 16 (random) query strings generated at the manager.
# Both the downloads are automatically unpacked, cached, and shared
# with all the same tasks on the worker.

import ndcctools.taskvine as vine
import random
import argparse
import getpass
import os
# Permitted letters in an amino acid sequence


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="vine_example_bbduk.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--B', action="store_true")
    parser.add_argument('--C', action="store_true")
    args = parser.parse_args()


    m = vine.Manager(port=9123)
    m.set_name("qcreport")
    # m.declare_poncho('env.tar.gz')
    
    print("Declaring files...")

    run_content = None
    run_path = "/afs/crc.nd.edu/user/j/jzhou24/ND_BIO60132_Sp24/Analysis/01_Quality_Control/single_read/test.sh"
    run_path = "./test.sh"
    with open(run_path, "r") as f:
        run_content = f.read()
    run = m.declare_buffer(run_content)
    
    fastqc = m.declare_file("../../Software/FastQC/")   
    bbmap = m.declare_file("../../Software/bbmap/")
    sickle = m.declare_file("../../Software/sickle/")
    trimmomatic = m.declare_file("../../Software/Trimmomatic-0.39/")
    references = m.declare_file("../References")
    # phix_adapters = m.declare_file("../References/phix_adapters.fa.gz")

    print("Declaring tasks...")

    directory = "../../../Input_data/"
    for filename in os.listdir(directory):
        fullpath = os.path.join(directory, filename)
        if os.path.isfile(fullpath) and fullpath.endswith("fastq"):
            
            ############ B
            if args.B:
                input_filename = os.path.basename(fullpath)
                input_data = m.declare_file(fullpath)
                
                t = vine.Task(
                    command = f"""
                        ./fastqc/fastqc -o ./fastqc {filename}
                    """,

                    inputs = {
                        fastqc: {"remote_name": "fastqc"},

                        input_data: {"remote_name": input_filename},
                    },
                )
                output_filename = input_filename.replace('.fastq', '_fastqc.zip')
                output_data = m.declare_file(output_filename)
                t.add_output(output_data, f"./fastqc/{output_filename}", watch=True)
        
                task_id = m.submit(t)
                print(f"submitted task {t.id}: {t.command}")

            
            ############ C
            if args.C:
                input_filename = os.path.basename(fullpath)
                input_data = m.declare_file(fullpath)
                input_filename = os.path.basename(fullpath)
                if "R1" in input_filename:
                    base, _ = os.path.splitext(input_filename)
                    base = base.split("_R1")[0]

                    trim_output_filename = f"{base}_1.trimclean.fq"
                    trim_log_filename = f"{base}.trimlog.txt"

                    sickle_output_filename = f"{base}_1.trimclean.sickleclean.fq"

                    bbduk_output_filename = sickle_output_filename.split("fq")[0] + "spikeclean.fq"

                    bbwrap_output_filename = bbduk_output_filename.split("fq")[0] + "hostclean.fq"

                    final_output_filename = bbduk_output_filename.split("_1.trim")[0] + "_1.final.clean.fq"
                    final_outputm_filename = bbduk_output_filename.split("_1.trim")[0] + "_1.reads_that_match_rRNA.fq"

                    t = vine.Task(
                        inputs = {
                            trimmomatic: {"remote_name": "trimmomatic"},
                            sickle: {"remote_name": "sickle"},
                            bbmap: {"remote_name": "bbmap"},
                            references: {"remote_name": "references"},
                            # phix_adapters: {"remote_name": "phix_adapters.fa.gz"},

                            input_data: {"remote_name": input_filename},
                        },

                        command = f"""
                            java -jar ./trimmomatic/trimmomatic-0.39.jar SE {input_filename} \
                            -threads 16 \
                            -trimlog {trim_log_filename} \
                            {trim_output_filename} \
                            ILLUMINACLIP:../References/adaptors.fa:1:50:30 \
                            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60
                            
                            ./sickle/sickle se \
                            -n \
                            -f {trim_output_filename} \
                            -o {sickle_output_filename} \
                            -t sanger \
                            -q 20 \
                            -l 60
                            
                            ./bbmap/bbduk.sh \
                            threads=8 \
                            in={sickle_output_filename} \
                            k=31 \
                            ref=./references/phix_adapters.fa.gz \
                            out1={bbduk_output_filename} \
                            minlength=60

                            ./bbmap/bbwrap.sh \
                            threads=16 \
                            minid=0.95 \
                            maxindel=3 \
                            bwr=0.16 \
                            bw=12 \
                            quickmatch \
                            fast \
                            minhits=2 \
                            qtrim=rl \
                            trimq=20 \
                            minlength=60 \
                            in={bbduk_output_filename} \
                            outu1={bbwrap_output_filename} \
                            path=./references/human

                            ./bbmap/bbduk.sh \
                            in={bbwrap_output_filename} \
                            ref=./references/smr_v4.3_default_db.fasta \
                            out={final_output_filename} \
                            outm={final_outputm_filename}

                        """,
                    )
                    

                trim_output = m.declare_file(trim_output_filename)
                trim_log = m.declare_file(trim_log_filename)
                sickle_output = m.declare_file(sickle_output_filename)
                bbduk_output = m.declare_file(bbduk_output_filename)
                bbwrap_output = m.declare_file(bbwrap_output_filename)
                final_output = m.declare_file(final_output_filename)
                final_outputm = m.declare_file(final_outputm_filename)

                t.add_output(trim_output, trim_output_filename, watch=True)
                t.add_output(trim_log, trim_log_filename, watch=True)
                t.add_output(sickle_output, sickle_output_filename, watch=True)
                t.add_output(bbduk_output, bbduk_output_filename, watch=True)
                t.add_output(bbwrap_output, bbwrap_output_filename, watch=True)
                t.add_output(final_output, final_output_filename, watch=True)
                t.add_output(final_outputm, final_outputm_filename, watch=True)

                task_id = m.submit(t)
                print(f"submitted task {t.id}: {t.command}")

            ############
            
                
    

    print(f"TaskVine listening for workers on {m.port}")

    print("Waiting for tasks to complete...")
    while not m.empty():
        t = m.wait(5)
        if t:
            if t.successful():
                print(f"task {t.id} result: {t.std_output}")
            elif t.completed():
                print(
                    f"task {t.id} completed with an executin error, exit code {t.exit_code}"
                )
            else:
                print(f"task {t.id} failed with status {t.result}")

    print("all tasks complete!")
# vim: set sts=4 sw=4 ts=4 expandtab ft=python:

