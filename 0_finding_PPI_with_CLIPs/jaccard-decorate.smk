from pathlib import Path
# copied from notebook, just to find input
skipper_dir = Path('/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/20220728_encode3/k562/output/reproducible_enriched_windows')
rbps_to_include = list(skipper_dir.glob('SF3B*')
                )+list(skipper_dir.glob('FXR*')
                )+list(skipper_dir.glob('FMR*')
                )+list(skipper_dir.glob('IGF*'))
# here we are to include some HepG2 CLIPs (processed with different annotations from above)
skipper_dir = Path('/projects/ps-yeolab3/eboyle/encode/pipeline/05_20220720/20220728_encode3/hepg2/output/reproducible_enriched_windows')
rbps_to_include += list(skipper_dir.glob('FXR*'))

print(rbps_to_include)



# extract sample ID, list comprehension again with zip() paring ID to path and dict() converting a list of items pairs to a dictionary
name_to_file = dict(
                    zip(
                        [f.name.split('.')[0] for f in rbps_to_include],
                        [str(f) for f in rbps_to_include]
                        )
                    ) # dictionary pointing ID 'SF3B4_K562' -> path '/.../SF3B4_K562.reproducible_windows.bed'
print(name_to_file)

# get combinations of two RBPs
from itertools import combinations
rbp_combinations = list(combinations(name_to_file.keys(), 2))
print(rbp_combinations)

import os
try:
    os.mkdir('error_files')
    os.mkdir('stdout')
except:
    pass

# rule "all" is the special rule that defines all the intended outputs
rule all:
    input:
        expand("output/{sample_combinations}.jaccard.txt",
                sample_combinations = [f'{i[0]}.{i[1]}' for i in rbp_combinations]
        )

# now we are to implement stuffs that gets us the jaccard index
# we want to use bedtools jaccard
# https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html
# which takes two bed files
# but skipper's reproducible enriched window is a tsv with HEADER, which is NOT a bed file
# so we gonna make it bed-ish

rule make_skipper_become_bed:
    input:
        lambda wildcards: name_to_file[wildcards.sample_label]
    output:
        "output/{sample_label}.bed"
    params:
        error_out_file = "error_files/make_bed.{sample_label}",
        out_file = "stdout/make_bed.{sample_label}",
        run_time = "00:05:00",
        cores = 1,
    shell:
        """
        module load bedtools
        zcat {input} |  tail -n +2 | sort -k1,1 -k2,2n -k3,3n > {output}
        """

rule jaccard:
    input:
        sample_one_bed = "output/{sample_label_1}.bed",
        sample_two_bed = "output/{sample_label_2}.bed"
    output:
        "output/{sample_label_1}.{sample_label_2}.jaccard.txt"
    params:
        error_out_file = "error_files/jaccard.{sample_label_1}.{sample_label_2}",
        out_file = "stdout/jaccard.{sample_label_1}.{sample_label_2}",
        run_time = "00:05:00",
        cores = 1,
    shell:
        """
        module load bedtools
        bedtools jaccard -a {input.sample_one_bed} -b {input.sample_two_bed} -F 0.5 > {output}
        """
