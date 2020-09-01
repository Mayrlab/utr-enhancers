#!/home/fanslerm/software/miniconda3/bin/python
import os
import sys

from snakemake.utils import read_job_properties

_DEFAULT = {
    "queue": "cpuqueue",
    "threads": 1,
    "mem": 8,
    "walltime": 4
    }

# Extract job properties
job_script = sys.argv[1]
job_properties = read_job_properties(job_script)
print(job_properties)

rule_name = job_properties["rule"]
threads = job_properties.get('threads', _DEFAULT['threads'])
job_resources = job_properties.get('resources', _DEFAULT)
queue = job_resources.get('queue', _DEFAULT['queue'])
mem = job_resources.get('mem', _DEFAULT['mem'])
walltime = job_resources.get('walltime', _DEFAULT['walltime'])

# Logs
os.makedirs("logs/lsf/%s" % rule_name, exist_ok=True)
stdlog = "logs/lsf/{}/%J.stdout".format(rule_name)
errlog = "logs/lsf/{}/%J.err".format(rule_name)

cmd_str = 'bsub -q {queue} -n {threads} -R"span[hosts=1] rusage[mem={mem}]" -W {walltime}:00 -o {stdlog} -e {errlog} -J {job_name} {script}'.format(
    queue=queue, threads=threads, mem=mem, walltime=walltime, stdlog=stdlog, errlog=errlog, script=job_script, job_name=rule_name)

print("System call:\n>%s" % cmd_str)

os.system(cmd_str)
