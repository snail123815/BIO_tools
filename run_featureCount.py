from subprocess import run
import os
for f in os.listdir('alignment/NCBI_genome'):
    if not f.endswith('.sam'):
        continue
    b = os.path.splitext(os.path.split(f)[-1])[0]
    args = ['featureCounts', '-T', '8', '-a', 'alignment/genomeBuild/NZ_CP022438_S.peucetius_subsp.caesius_ATCC_27952_chromosome.gtf', '-t', 'gene', '-g', 'locus_tag', '-o', f'genecounts/{b}.txt', f'alignment/NCBI_genome/{f}']
    print(' '.join(args))
    res = run(args, capture_output=True)
    if res.returncode != 0:
        print(res.stdout.decode())
        print(res.stderr.decode())
        break
