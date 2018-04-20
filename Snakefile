def load_datasets():
    import os

    dataset = config['input']
    dat_dir = os.path.dirname(dataset)

    with open(dataset) as f:
        datasets = [line.strip().split() for line in f if not line.startswith('#') if not line.startswith('!')]
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
        '%s/found.csv' % OUTP
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


rule merge_reads:
    input:
        errors='%s/data/{name}/errors.csv' % OUTP,
        labels='%s/data/{name}/labels.csv' % OUTP
    output:
        merged='%s/data/{name}/merged.csv' % OUTP,
        coverage='%s/data/{name}/segment_coverage.csv' % OUTP
    message:
        'Merging reads: {wildcards.name}'
    shell:
        '%s/merge_reads.py -e {input.errors} -l {input.labels} -o {output.merged} -c {output.coverage}' % DIR


rule filter_errors:
    input:
        errors='%s/data/{name}/merged.csv' % OUTP,
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


rule divide_reads:
    input:
        filtered='%s/data/{name}/filtered.csv' % OUTP,
        bin='$s/bin/divide_reads' % DIR 
    output:
        reads='%s/data/{name}/subset_reads.csv' % OUTP,
        summary='%s/data/{name}/subset_summary.csv' % OUTP
    message:
        'Dividing read sets: {wildcards.name}'
    shell:
        '{dir}/divide_reads.py -i {{input.filtered}} -r {{output.reads}} -s {{output.summary}} ' \
        '-S {significance} -p {pairs} -c {coverage}' \
            .format(dir=DIR, significance=config['significance'], pairs=config['pairs'],
                    coverage=config['subset_coverage'])


rule add_sequences:
    input:
        summary='%s/data/{name}/subset_summary.csv' % OUTP,
        germline=config['segments']
    output:
        full='%s/data/{name}/full.fa' % OUTP,
        cropped='%s/data/{name}/cropped.fa' % OUTP
    message:
        'Converting subsets into fasta: {wildcards.name}'
    shell:
        '{dir}/add_sequences.py -i {{input.summary}} -s {{input.germline}} -o {{output.full}} {{output.cropped}} ' \
        '-r {range[0]} {range[1]} --keep-empty -c 0 -l 0' \
            .format(dir=DIR, range=config['range'])


rule fitting_alignment:
    input:
        cropped='%s/data/{name}/cropped.fa' % OUTP,
        germline=config['segments']
    output:
        '%s/data/{name}/fit.csv' % OUTP
    message:
        'Performing fitting alignment on read subsets: {wildcards.name}'
    shell:
        '{dir}/fitting_alignment.py -i {{input.cropped}} -s {{input.germline}} -o {{output}} ' \
        '-M {mismatch} -O {open_gap} -E {extend_gap}' \
            .format(dir=DIR, mismatch=config['mismatch'], open_gap=config['open_gap'], extend_gap=config['extend_gap'])


rule summarize:
    input:
        subsets='%s/data/{name}/subset_summary.csv' % OUTP,
        fit='%s/data/{name}/fit.csv' % OUTP,
        full='%s/data/{name}/full.fa' % OUTP,
        cropped='%s/data/{name}/cropped.fa' % OUTP
    output:
        '%s/data/{name}/summary.csv' % OUTP
    message:
        'Summarizing {wildcards.name}'
    shell:
        '{dir}/summarize.py -s {{input.subsets}} -f {{input.fit}} -S {{input.full}} {{input.cropped}} ' \
        '-o {{output}} -D 0 -I 0' \
            .format(dir=DIR)

rule combine_and_classify:
    input:
        summaries=expand('{outp}/data/{name}/summary.csv', outp=OUTP, name=NAMES),
        datasets=config['input']
    output:
        success='%s/found.csv' % OUTP,
        fail='%s/discarded.csv' % OUTP
    message:
        'Combining summaries and classifying putative segments'
    shell:
        '{dir}/combine_and_classify.py -i {{input.summaries}} -d {{input.datasets}} ' \
        '-o {{output.success}} {{output.fail}} {plots} --prob {prob} {isoclines}' \
            .format(dir=DIR, plots='' if config['skip_plots'] else '-p %s/plots' % OUTP,
                    prob=config['prob'],
                    isoclines='' if config['skip_plots'] else '--isoclines %s' % ' '.join(map(str, config['isoclines'])))

