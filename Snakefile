def load_datasets():
    import os

    dataset = config['datasets']
    dat_dir = os.path.dirname(dataset)

    with open(dataset) as f:
        datasets = [line.strip().split() for line in f if not line.startswith('#')]
        names = [item[0] for item in datasets]
        alignment = [os.path.join(dat_dir, item[1]) for item in datasets]
        cdrs = [os.path.join(dat_dir, item[2]) for item in datasets]

    return names, alignment, cdrs

NAMES, ALIGNMENTS, CDRS = load_datasets()
DIR = '.' #os.path.dirname(__file__)
OUTP = config['output']


def fail():
    print('Error')
    exit(1)


rule all:
    input:
        expand('{outp}/data/{name}/candidate.csv', outp=OUTP, name=NAMES)


rule create_inputs:
    input:
        alignments=ALIGNMENTS,
        cdrs=CDRS
    output:
        alignments=expand('{outp}/data/{name}/alignment.csv', outp=OUTP, name=NAMES),
        cdrs=expand('{outp}/data/{name}/cdrs.csv', outp=OUTP, name=NAMES)
    run:
        import subprocess
        import os

        for i in range(len(NAMES)):
            command = ['ln', '-s', os.path.abspath(input.alignments[i]), output.alignments[i]]
            if subprocess.run(command).returncode:
                fail()

            command2 = ['ln', '-s', os.path.abspath(input.cdrs[i]), output.cdrs[i]]
            if subprocess.run(command2).returncode:
                fail()
    

rule errors:
    input:
        '%s/data/{name}/alignment.csv' % OUTP
    output:
        '%s/data/{name}/errors.csv' % OUTP
    shell:
        '%s/alignment_errors.py -i {input} -o {output}' % DIR


rule coverage:
    input:
        '%s/data/{name}/alignment.csv' % OUTP
    output:
        '%s/data/{name}/segment_coverage.csv' % OUTP
    shell:
        r"cat {input} | sed -n '3~4p' | sed 's/.*GENE:\([^|]*\)|.*/\1/' | sort | uniq -c "
        r"| awk '{{t = $1; $1 = $2; $2 = t; print}}' > {output}"

rule filter:
    input:
        errors='%s/data/{name}/errors.csv' % OUTP,
        coverage='%s/data/{name}/segment_coverage.csv' % OUTP
    output:
        '%s/data/{name}/filtered.csv' % OUTP
    shell:
        '%s/filter_errors.py -e {input.errors} -c {input.coverage} -o {output}' % DIR

rule candidate:
    input:
        '%s/data/{name}/filtered.csv' % OUTP
    output:
        '%s/data/{name}/candidate.csv' % OUTP
    shell:
        '%s/candidate_polymorphisms_gap.py -i {input} -o {output} -g 0.4' % DIR
