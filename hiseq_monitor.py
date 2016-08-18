"""
Data conversion script for Illumina HiSeq short read runs

Set up to be routinely executed by a cron job at qhatever frequency required
by your organization's level of throughput

The script will check the output RTA directory of any unconverted Illumina 
HiSeq run. A file in the conversion output directory has the history of all 
previously converted runs.
If the run is found to be new then two files that signify RTA transfer is
 complete are checked. Once all data is deemed transfered, a shell script is
generated in the user's tmp directory. This shell script is submit to a cluster
via qsub (open Grid Scheduler) and monitored for successful return status.

The shell script adds to the path the required bcl2fastq script, change dir to
the raw data dir, and then executes the bcl2fastq program. Parameters can be 
changed by modifying the cmd in this script.

Upon successful return status, the name of the directory is added to the history
file so that upon next invocation of this monitor script there is no redundant
computing.

A log file with any HiSeq runs sent to the grid during each time the script is
executed also gets updated.

Crontab:
# run hiseq monitor every 30 minutes to convert bcl2fastq
2 * * * * /usr/bin/python /path/to/hiseq_monitor.py

"""

import os
import sys
import subprocess
import datetime

HISEQ_DIR = "/path/out/raw/HiSeq/"
HISEQ_DATA_DIR = "/path/to/hiseq_fastq/"
analyzed_file = HISEQ_DATA_DIR + "hiseq_complete.txt"

analyzed = []

with open(analyzed_file, 'r') as h:
	analyzed = h.readlines()
	analyzed = [x.rstrip() for x in analyzed]

hiseq_files = os.listdir(HISEQ_DIR)
fastq_files = [x for x in analyzed]

for hiseq_file in hiseq_files: 
	if hiseq_file.startswith("."):
		hiseq_files.remove(hiseq_file)

for hiseq_file in fastq_files:
	if hiseq_file.startswith("."):
        	fastq_files.remove(hiseq_file)

status = 1
for hiseq_file in hiseq_files:
	if hiseq_file not in fastq_files:
		check_1, check_2 = False, False
		for f in os.listdir(HISEQ_DIR + hiseq_file):
			if f == 'Basecalling_Netcopy_complete.txt':
				check_1 = True
			if f == 'ImageAnalysis_Netcopy_complete.txt':
				check_2 = True
		if check_1 and check_2:
			print(hiseq_file)
			conversion_file = hiseq_file
			tmpdir = "/hptmp/acristo/"
			cmd = "/add_to_path .bcl2fastq2-2.17.1.14\n"
			cmd += "cd {hs}".format(hs = HISEQ_DIR) + hiseq_file + "/\n"
			cmd += "bcl2fastq --with-failed-reads -o {hs}{fastq_dir} --sample-sheet {raw}{in_dir}/Samplesheet.csv --input-dir {raw}{in_dir} -p 8 -d 8".format(hs = HISEQ_DIR, raw = HISEQ_DATA_DIR, fastq_dir = hiseq_file, in_dir = hiseq_file)
			with open(tmpdir + "{in_dir}.sh".format(in_dir = hiseq_file), 'w') as op:
				op.write(cmd)
			subcmd = "qsub -sync y -q short -o " + tmpdir + " -e " + tmpdir + " " + tmpdir + hiseq_file +".sh"
			status = subprocess.call(subcmd.split())
			subprocess.call("mkdir " + HISEQ_DATA_DIR + "{fastq_dir}/logs/".format(fastq_dir = hiseq_file))
			print(status)
		if status == 0:
			with open(analyzed_file, 'a') as af:
				af.write(conversion_file + "\n")
	else:
		continue

monitor_log = HISEQ_DATA_DIR + "hiseq_monitor.log"
with open(monitor_log, 'a') as m: 
	time = str(datetime.datetime.now())
	if status == 0:
		m.write("Did conversion on {file} at ".format(file = conversion_file) + time + "\n")
	else:
		m.write("Found no new HiSeq run at " + time + "\n")
