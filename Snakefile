def load_datasets():
    import os

    dataset = config['datasets']
    dat_dir = os.path.dirname(dataset)

    with open(dataset) as f:
        datasets = [line.strip().split() for line in f if not line.startswith('#')]
        names = [item[0] for item in datasets]
        alignment = [os.path.abspath(os.path.join(dat_dir, item[1])) for item in datasets]
        cdrs = [os.path.abspath(os.path.join(dat_dir, item[2])) for item in datasets]

    return names, alignment, cdrs

NAMES, ALIGNMENTS, CDRS = load_datasets()
DIR = '.' #os.path.dirname(__file__)
OUTP = config['output']
import os


def fail():
    print('Error')
    exit(1)


rule all:
    input:
#        expand('{outp}/data/{name}/expanded.{iterations}.csv', outp=OUTP, name=NAMES, iterations=config['iterations']),
#        expand('{outp}/data/{name}/labels.csv', outp=OUTP, name=NAMES),
        expand('{outp}/data/{name}/combined.csv', outp=OUTP, name=NAMES)
    threads: 16


def element_by_name(arr):
    import os
    def arr_el(wildcards):
        return arr[NAMES.index(wildcards['name'])]
    return arr_el


rule create_alignment_input:
    input:
        element_by_name(ALIGNMENTS)
    output:
        '%s/data/{name}/alignment.fa' % OUTP
    shell:
        'ln -s {input} {output}'


rule create_cdrs_input:
    input:
        element_by_name(CDRS)
    output:
        '%s/data/{name}/cdrs.fa' % OUTP
    shell:
        'ln -s {input} {output}'


rule errors:
    input:
        '%s/data/{name}/alignment.fa' % OUTP
    output:
        '%s/data/{name}/errors.csv' % OUTP
    shell:
        '%s/alignment_errors.py -i {input} -o {output}' % DIR


rule coverage:
    input:
        '%s/data/{name}/alignment.fa' % OUTP
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
        '%s/data/{name}/expanded.0.csv' % OUTP
    shell:
        '%s/candidate_polymorphisms_gap.py -i {input} -o {output} -g 0.4' % DIR
        

def expansion_name(wildcards):
    return '%s/data/%s/expanded.%d.csv' % (OUTP, wildcards['name'], int(wildcards['num']) - 1)


rule expand:
    input:
        filtered='%s/data/{name}/filtered.csv' % OUTP,
        prev=expansion_name
    output:
        '%s/data/{name}/expanded.{num}.csv' % OUTP
    wildcard_constraints:
        num='[1-9]\d*'
    shell:
        '%s/expand_candidate.py -f {input.filtered} -c {input.prev} -o {output}' % DIR

rule compile:
    input:
        os.path.abspath('%s/src/{name}.cpp' % DIR)
    output:
        os.path.abspath('%s/bin/{name}' % DIR)
    shell:
        'g++ {input} -std=c++1y -O3 -o {output}'


rule cdr_hamming:
    input:
        bin=os.path.abspath('%s/bin/hamming' % DIR),
        f='%s/data/{name}/cdrs.fa' % OUTP
    output:
        '%s/data/{name}/labels.csv' % OUTP
    shell:
        r"cat {input.f} | sed 's/^>\(.*\)$/>\1\t/' | tr -d '\n' | tr '>' '\n' | %s/hamming_labels.py -i - -o {output}" % DIR


#rule create_datasets_file:
#    output:
#        '%s/data/datasets.txt' % OUTP
#    run:
#        with open(output[0], 'w') as outp:
#            outp.write('name\talignment\tlabels\n')
#            for name in NAMES:
#                outp.write('{name}\t{name}/filtered.csv\t{name}/labels.csv\n'.format(name=name))

rule combine_candidate:
    input:
        filtered='%s/data/{name}/filtered.csv' % OUTP,
        candidate='%s/data/{name}/expanded.%d.csv' % (OUTP, config['iterations']),
        labels='%s/data/{name}/labels.csv' % OUTP,
        segments=config['segments']
    output:
        '%s/data/{name}/combined.csv' % OUTP
    shell:
        '%s/combine_errors.py -f {input.filtered} -c {input.candidate} -l {input.labels} '
        '-s {input.segments} -o {output}' % DIR

