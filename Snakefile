def load_datasets():
    import os

    dataset = config['input']
    dat_dir = os.path.dirname(dataset)

    with open(dataset) as f:
        datasets = [line.strip().split() for line in f if not line.startswith('#')]
        names = [item[0] for item in datasets]
        alignment = [os.path.abspath(os.path.join(dat_dir, item[1])) for item in datasets]
        cdrs = [os.path.abspath(os.path.join(dat_dir, item[2])) for item in datasets]

    return names, alignment, cdrs

NAMES, ALIGNMENTS, CDRS = load_datasets()
DIR = config['script_dir']
OUTP = config['output']
import os


rule all:
    input:
        expand('{outp}/data/merged.csv', outp=OUTP)
    threads: config['threads']


def element_by_name(arr):
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


rule alignment_errors:
    input:
        '%s/data/{name}/alignment.fa' % OUTP
    output:
        '%s/data/{name}/errors.csv' % OUTP
    message:
        'Alignment fasta to csv: {wildcards.name}'
    shell:
        '%s/alignment_errors.py -i {input} -o {output}' % DIR


rule segment_coverage:
    input:
        '%s/data/{name}/alignment.fa' % OUTP
    output:
        '%s/data/{name}/segment_coverage.csv' % OUTP
    message:
        'Calculating segment coverage: {wildcards.name}'
    shell:
        r"cat {input} | sed -n '3~4p' | sed 's/.*GENE:\([^|]*\)|.*/\1/' | sort | uniq -c "
        r"| awk '{{t = $1; $1 = $2; $2 = t; print}}' > {output}"


rule filter_errors:
    input:
        errors='%s/data/{name}/errors.csv' % OUTP,
        coverage='%s/data/{name}/segment_coverage.csv' % OUTP
    output:
        '%s/data/{name}/filtered.csv' % OUTP
    message:
        'Filtering errors: {wildcards.name}'
    shell:
        '{dir}/filter_errors.py -e {{input.errors}} -c {{input.coverage}} -o {{output}} ' \
        '--range {range[0]} {range[1]} --threshold {threshold}{keep_n}'.format(dir=DIR,
                                                                       range=config['range'],
                                                                       threshold=config['segment_coverage'],
                                                                       keep_n=' --keep-n' if config['keep_n'] else '')

rule compile:
    input:
        os.path.abspath('%s/src/{name}.cpp' % DIR)
    output:
        os.path.abspath('%s/bin/{name}' % DIR)
    message:
        'Compiling {wildcards.name}.cpp'
    shell:
        'g++ {input} -std=c++1y -O3 -o {output}'


rule hamming_labels:
    input:
        bin=os.path.abspath('%s/bin/hamming' % DIR),
        f='%s/data/{name}/cdrs.fa' % OUTP
    output:
        '%s/data/{name}/labels.csv' % OUTP
    message:
        'Merging read labels using Hamming graph: {wildcards.name}'
    shell:
        r"cat {{input.f}} | sed 's/^>\(.*\)$/>\1\t/' | tr -d '\n' | tr '>' '\n' | " \
        '{dir}/hamming_labels.py -i - -o {{output}} --tau {tau};\n' \
        r"sed -i '3~1s/\t/\t{{wildcards.name}}_/' {{output}}".format(dir=DIR, tau=config['labels_tau'])
    

rule combine_errors:
    input:
        filtered='%s/data/{name}/filtered.csv' % OUTP,
        labels='%s/data/{name}/labels.csv' % OUTP,
        segments=config['segments']
    output:
        '%s/data/{name}/combined.csv' % OUTP
    message:
        'Combine errors: {wildcards.name}'
    shell:
        '{dir}/combine_errors.py -f {{input.filtered}} -l {{input.labels}} '
        '-s {{input.segments}} -o {{output}} --peaks {peaks} '
        '--coverage {coverage} --spl-coverage {spl_coverage} --similarity {similarity}' \
            .format(dir=DIR, coverage=config['subset_coverage'], peaks=config['peaks'],
                    spl_coverage=config['split_coverage'], similarity=config['similarity'])


rule clip_and_filter:
    input:
        bin=os.path.abspath('%s/bin/clip_and_filter' % DIR),
        combined='%s/data/{name}/combined.csv' % OUTP,
        segments=config['segments'] if config['clipped_segments'] is None else config['clipped_segments']
    output:
        '%s/data/{name}/clipped_filtered.csv' % OUTP
    message:
        'Clip and filter sequences: {wildcards.name}'
    shell:
        'tail -n+3 {{input.combined}} | cut -f4 | ' \
        '{dir}/clip_and_filter.py -s - -S {{input.segments}} -o {{output}} ' \
        '--range {range[0]} {range[1]} --min-dist {min_dist}{no_clipping}' \
                .format(dir=DIR,
                        range=config['range'],
                        min_dist=config['min_dist'],
                        no_clipping='' if config['clipped_segments'] is None else ' --no-clipping')


rule concat_sequences:
    input:
        expand('{outp}/data/{name}/clipped_filtered.csv', outp=OUTP, name=NAMES)
    output:
        '%s/data/concat.csv' % OUTP
    shell:
        'echo "seq\tlabel" > {output}\n' \
        'for f in {input}; do tail -n+3 $f >> {output}; done'


rule hamming_sequences:
    input:
        bin=os.path.abspath('%s/bin/hamming' % DIR),
        f='%s/data/concat.csv' % OUTP
    output:
        '%s/data/components.csv' % OUTP
    message:
        'Merging sequences using Hamming graph'
    shell:
        '{dir}/hamming_labels.py -i {{input.f}} -o {{output}} --tau {tau}{diff_lengths}' \
            .format(dir=DIR, tau=config['seqs_tau'],
                    diff_lengths='' if config['seqs_same_lengths'] else ' --accept-diff-lengths')


rule merge_and_sort:
    input:
        combined=expand('{outp}/data/{name}/combined.csv', outp=OUTP, name=NAMES),
        components='%s/data/components.csv' % OUTP
    output:
        '%s/data/merged.csv' % OUTP,
        '%s/data/merged_full.csv' % OUTP
    message:
        'Merging novel sequences and sorting them'
    shell:
        '{dir}/merge_and_sort.py -i {{input.combined}} -c {{input.components}} -o {{output}} --min-coverage {min_coverage}' \
            .format(dir=DIR, min_coverage=config['label_coverage'])

