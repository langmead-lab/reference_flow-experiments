import os
import altair as alt
import pandas as pd
import argparse

def is_snp(fields):
    return int(fields[2]) - int(fields[1]) == 1


def get_variant_set(filename):
    f = open(filename, 'r')
    
    # Each element in `variant_set` is represented as "<chr>_<start_pos>"
    variant_set = set()
    for line in f:
        fields = line.split()
        if is_snp(fields):
            variant_set.add(f'{fields[0]}_{fields[1]}')
        else:
            raise ValueError('Variant is not a SNP', line)
    return variant_set


def get_number_of_overlapped_variants(set_a, set_b):
    return len(set_a.intersection(set_b))


def get_variant_sets_and_fractions(method, repeats, prefix=''):
    set_all = get_variant_set(f'{prefix}{method}-NA12878-biased_sites.bed')
    repeat_sets = {}
    fractions = {}
    for _, repeat in enumerate(repeats):
        current_set = get_variant_set(f'{prefix}{method}-rmsk-{repeat}.bed')

        # Repeat classes should not overlap.
        for j_rc, other_set in enumerate(repeat_sets):
            num_overlap = get_number_of_overlapped_variants(current_set, other_set)
            if num_overlap != 0:
                raise ValueError('There are overlapping variants in two sets',
                                 repeat, list(repeat_sets.keys())[j_rc])
        repeat_sets[repeat] = current_set
        fractions[repeat] = len(repeat_sets[repeat]) / len(set_all)
    return repeat_sets, fractions


def summarize_repeat_classes(prefix, alignment_methods, output_class, map_methods):
    repeat_classes = ['SINE', 'LINE', 'LTR', 'Simple_repeat', 'Satellite']
    
    methods = []
    fractions = []
    numbers = []
    rclasses = []
    for method in alignment_methods:
        repeat_sets, fractions_of_each_class = get_variant_sets_and_fractions(
            method=method, repeats=repeat_classes, prefix=prefix)
        fractions.extend(list(fractions_of_each_class.values()))
        numbers.extend([len(s) for s in repeat_sets.values()])
        # Remove the space in "Simple repeat".
        rclasses.extend(['SINE', 'LINE', 'LTR', 'Simple repeat', 'Satellite'])
        methods.extend([map_methods[method]] * len(list(fractions_of_each_class.values())))

    df_repeat_classes = pd.DataFrame()
    df_repeat_classes['Number'] = numbers
    df_repeat_classes['Fraction'] = fractions
    df_repeat_classes['Repeat Class'] = rclasses
    df_repeat_classes['Alignment Method'] = methods
    df_repeat_classes.to_csv(output_class, sep='\t', index=None)

    # Y: number of biased HETs.
    alt.Chart(df_repeat_classes).mark_bar().encode(
        x=alt.X('Alignment Method:N', sort=alignment_methods, title=None),
        y='Number',
        color='Alignment Method:N',
        column=alt.Column('Repeat Class:N',
                        sort=['SINE', 'LINE', 'LTR', 'Satellite', 'Simple repeat'])
    )


def summarize_repeat_families(prefix, alignment_methods, output_family, map_methods):
    repeat_families = ['L1', 'Alu', 'ERVL', 'L2', 'ERV1', 'MIR']
    methods = []
    fractions = []
    numbers = []
    rfamilies = []
    for method in alignment_methods:
        repeat_sets, fractions_of_each_family = get_variant_sets_and_fractions(
            method=method, repeats=repeat_families, prefix=prefix)
        fractions.extend(list(fractions_of_each_family.values()))
        numbers.extend([len(s) for s in repeat_sets.values()])
        rfamilies.extend(repeat_families)
        methods.extend([map_methods[method]] * len(list(fractions_of_each_family.values())))

    df_repeat_families = pd.DataFrame()
    df_repeat_families['Number'] = numbers
    df_repeat_families['Fraction'] = fractions
    df_repeat_families['Repeat Family'] = rfamilies
    df_repeat_families['Alignment Method'] = methods
    df_repeat_families.to_csv(output_family, sep='\t', index=None)

    # Y: fraction of all biased HETs using the same method.
    alt.Chart(df_repeat_families).mark_bar().encode(
        x=alt.X('Alignment Method:N', sort=alignment_methods, title=None),
        y='Fraction',
        color='Alignment Method:N',
        column=alt.Column('Repeat Family:N', sort=repeat_families)
    )

    # Y: number of biased HETs.
    alt.Chart(df_repeat_families).mark_bar().encode(
        x=alt.X('Alignment Method:N', sort=alignment_methods, title=None),
        y='Number',
        color='Alignment Method:N',
        column=alt.Column('Repeat Family:N', sort=repeat_families)
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-p', '--prefix',
        help='Prefix of extracted BED files'
    )
    parser.add_argument(
        '-of', '--output_family',
        help='Output TSV file for repeat families'
    )
    parser.add_argument(
        '-oc', '--output_class',
        help='Output TSV file for repeat classes'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    alignment_methods = ['GRC', 'major', 'randflow_ld', 'vg', 'per']
    map_methods = {
        'GRC': 'GRC', 'major': 'Major',
        'randflow_ld': 'RandFlow-LD', 'vg': 'vg', 'per': 'Personalized'}
    summarize_repeat_classes(args.prefix, alignment_methods, args.output_class, map_methods)
    summarize_repeat_families(args.prefix, alignment_methods, args.output_family, map_methods)