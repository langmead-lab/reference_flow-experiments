rule build_individual_het:
    input:
        vcf = PREFIX_VCF_F + '.vcf'
    output:
        os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}.vcf')
    shell:
        '{BCFTOOLS} view -s {wildcards.INDIV} {input.vcf} | \
         {BCFTOOLS} view -i "AC>0" -v snps -g het -m2 -M2 > {output}'

rule get_het_with_indel:
    input:
        vcf = PREFIX_VCF_F + '.vcf'
    output:
        het = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_withindel.vcf')
    shell:
        '{BCFTOOLS} view -s {wildcards.INDIV} {input.vcf} | {BCFTOOLS} view -i "AC>0" -g het -m2 -M2 > {output.het}'

rule get_het_with_indel_processed:
    input:
        het = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_withindel.vcf')
    output:
        het = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf')
    shell:
        'cat {input.het} | python {DIR_SCRIPTS}/remove_het_overlapping_indel.py > {output.het}'

rule calc_major_bias:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_FIRST_PASS, '{}-h37maj-liftover-sorted.sam'.format(CHROM))
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'major-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'major-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'major-refbias.txt')
    shell:
        'echo "major" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
            -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_grc_bias:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        sam = os.path.join(DIR_FIRST_PASS, '{}-grch37-sorted.sam'.format(CHROM))
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'grch37-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'grch37-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'grch37-refbias.txt')
    shell:
        'echo "grch37" > {output.list_id};'
        'ls {input.sam} > {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
           -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_per_bias:
    input:
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        samA = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapA-liftover-sorted.sam'.format(CHROM)),
        samB = os.path.join(DIR_FIRST_PASS, '{}-per-merged-hapB-liftover-sorted.sam'.format(CHROM))
        # sam = os.path.join(DIR_FIRST_PASS, '{}-per-liftover-sorted.sam'.format(CHROM))
    output:
        list_path = os.path.join(DIR_FIRST_PASS, 'per-refbias.paths'),
        list_id = os.path.join(DIR_FIRST_PASS, 'per-refbias.ids'),
        bias = os.path.join(DIR_FIRST_PASS, 'per-refbias.txt')
    shell:
        'echo "perA" > {output.list_id};'
        'echo "perB" >> {output.list_id};'
        'ls {input.samA} > {output.list_path};'
        'ls {input.samB} >> {output.list_path};'
        '{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
            -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}'

rule calc_refflow_bias:
    input:
        # vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}.vcf'),
        vcf = os.path.join(DIR_FIRST_PASS, CHROM + '_{INDIV}_het_no_overlaps.vcf'),
        maj = os.path.join(DIR_FIRST_PASS,
            '{}-h37maj-mapqgeq{}-liftover-sorted.sam'.format(CHROM, ALN_MAPQ_THRSD)),
        second_maj = os.path.join(DIR_SECOND_PASS, '2ndpass-h37maj-liftover-sorted.sam'),
        second_pop = [os.path.join(DIR_SECOND_PASS,'2ndpass-') + 
            g + '-liftover-sorted.sam' for g in GROUP]
    output:
        list_path = os.path.join(DIR_SECOND_PASS, 'refflow-{}.paths'.format(POP_DIRNAME)),
        list_id = os.path.join(DIR_SECOND_PASS, 'refflow-{}.ids'.format(POP_DIRNAME)),
        bias = os.path.join(DIR_SECOND_PASS, 'refflow-{}.txt'.format(POP_DIRNAME))
    run:
        #: prepare list_path and list_id
        second_pass_group = ['h37maj']
        for g in GROUP:
            second_pass_group.append(g)
        with open(output.list_id, 'w') as f:
            for g in second_pass_group:
                f.write(g + '\n')
            f.write('h37maj-mapqgeq{}\n'.format(ALN_MAPQ_THRSD))
        list_second_pass_lifted_sam = [
            # os.path.join(DIR_SECOND_PASS, '2ndpass-') + g +
            os.path.join(DIR, 'experiments/' + wildcards.INDIV +
            '/' + POP_DIRNAME + '/2ndpass-') + g +
            '-liftover-sorted.sam' for g in second_pass_group]
        with open(output.list_path, 'w') as f:
            for s in list_second_pass_lifted_sam:
                # sys.stderr.write(s + '\n')
                f.write(s + '\n')
            f.write(input.maj + '\n')
        shell('cat {output.list_path};')
        shell('{PYTHON} {DIR_SCRIPTS}/refbias/lift_ref_flow.py -v {input.vcf} \
               -s {output.list_path} -n {output.list_id} -f {GENOME} -o {output.bias}')

rule check_refbias:
    input:
        refflow = expand(os.path.join(DIR_SECOND_PASS,
            'refflow-{}.txt'.format(POP_DIRNAME)), INDIV = INDIV),
        grch37 = expand(os.path.join(DIR_FIRST_PASS,
            'grch37-refbias.txt'), INDIV = INDIV),
        major = expand(os.path.join(DIR_FIRST_PASS,
            'major-refbias.txt'), INDIV = INDIV),
        per = expand(os.path.join(DIR_FIRST_PASS,
            'per-refbias.txt'), INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'refbias.done')))
        
