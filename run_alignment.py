from subprocess import run
import os

for f in os.listdir():
    if not f.endswith('.fastq.gz'):
        continue
    base = f.split('_R1')[0]
    target = f'alignment/NCBI_genome/{base}_NCBI.sam'
    if os.path.isfile(target):
        continue
    args = ['bowtie2', '-x', 'alignment/genomeBuild/SPE_NCBI', '-p', '8', '-U', f, '-S', target] 
    print(' '.join(args))
    result = run(args=args, capture_output=True)
    if result.returncode != 0:
        print(result.stderr.decode())
        print(result.stdout.decode())
        exit()
#for f in os.listdir():
#    if not f.endswith('.fastq.gz'):
#        continue
#    base = f.split('_R1')[0]
#    args = ['bowtie2', '-x', 'alignment/genomeBuild/SPE_selfAssembly', '-p', '8', '-U', f, '-S', f'alignment/selfAssembly_genome/{base}_selfAssembly.sam']
#    print(' '.join(args))
#    result = run(args=args, capture_output=True)
#    if result.returncode != 0:
#        print(result.stderr.decode())
#        print(result.stdout.decode())
#        exit()
#
