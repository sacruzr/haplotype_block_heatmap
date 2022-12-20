import pandas as pd

# load configuration file
configfile: "./config.yaml"

# read roi specification data
rois = pd.read_csv(config['rois_table_path'])
rois.set_index('roi_name', drop=False, inplace=True)

print(rois)

# get target output files
def get_final_output(wildcards):
    roi_names = rois['roi_name'].tolist()
    return ['{out_dir}/results/plots/{roi_name}.pdf'.format(
            out_dir = config['out_dir'],
            roi_name = roi_name    
            ) for roi_name in roi_names] 

# get the corresponding gwas out file according to roi
def get_gwas_output_by_roi(wildcards):
    gwas_file = rois.loc[wildcards.roi_name, 'gwas_out_filename']
    return '{gwas_dir}/{gwas_file}'.format(gwas_dir=config['gwas_dir'], gwas_file=gwas_file)
    


rule all:
    input: get_final_output

rule plot_GWAS_LD_heatmap_by_roi:
    output: '{out_dir}/results/plots/{{roi_name}}.pdf'.format(out_dir=config['out_dir'])
    input:
        haplotype_ld = '{out_dir}/haploview/target_ld_blocks/{{roi_name}}.ld/'.format(out_dir=config['out_dir']),
        gwas_out = get_gwas_output_by_roi

# split vcf arround roi to detect ld blocks
rule delimitate_vcf_by_roi:
    output: '{out_dir}/vcfs/window_roi/{{roi_name}}_splitted.vcf.gz'.format(out_dir=config['out_dir'])
    input:
        vcf = config['vcf_path']
    params:
        chrom = lambda wildcards: rois.loc[wildcards.roi_name, 'Chrom'],
        i_pos = lambda wildcards: rois.loc[wildcards.roi_name, 'Pos']  - int(rois.loc[wildcards.roi_name, 'initial_window']/2),
        e_pos = lambda wildcards: rois.loc[wildcards.roi_name, 'Pos']  + int(rois.loc[wildcards.roi_name, 'initial_window']/2)
    shell:
        """
        if [ ! -f {config[vcf_path]}.tbi ]; then
            tabix -p vcf {input.vcf}
        fi        

        bcftools view -r {params.chrom}:{params.i_pos}-{params.e_pos} -Oz -o {output} {input.vcf}
        """



# convert splitted vcf into Haploview format
rule convert_splitted_to_haploview_by_roi:
    output: 
        ped = '{out_dir}/haploview/data/{{roi_name}}.ped'.format(out_dir=config['out_dir']),
        info = '{out_dir}/haploview/data/{{roi_name}}.info'.format(out_dir=config['out_dir'])
    input:
        vcf = '{out_dir}/vcfs/window_roi/{{roi_name}}_splitted.vcf.gz'.format(out_dir=config['out_dir'])
    shell:
        """
        plink --vcf {input.vcf} \
        --double-id \
        --allow-extra-chr \
        --snps-only just-acgt \
        --recode HV \
        --out {config[out_dir]}/haploview/data/{wildcards.roi_name} && \
        mv {config[out_dir]}/haploview/data/{wildcards.roi_name}.chr*.ped {output.ped}
        mv {config[out_dir]}/haploview/data/{wildcards.roi_name}.chr*.info {output.info}
        rm {config[out_dir]}/haploview/data/{wildcards.roi_name}.nosex 
        rm {config[out_dir]}/haploview/data/{wildcards.roi_name}.log
        """
# detect LD blocks usng Haploview
rule get_ld_blocks_by_roi:
    output:'{out_dir}/haploview/ld_blocks/{{roi_name}}.{hv_method}/'.format(out_dir=config['out_dir'], hv_method = config['haploview_method'])
    input:
        ped = '{out_dir}/haploview/data/{{roi_name}}.ped/'.format(out_dir=config['out_dir']),
        info = '{out_dir}/haploview/data/{{roi_name}}.info/'.format(out_dir=config['out_dir'])
    shell:
        

# select the block where roi is located
rule get_target_ld_block_by_roi:
    output: '{out_dir}/haploview/target_ld_blocks/{{roi_name}}.ld/'.format(out_dir=config['out_dir'], hv_method = config['haploview_method'])
    input:
        vcf = '{out_dir}/vcfs/window_roi/{{roi_name}}_splitted.vcf.gz'.format(out_dir=config['out_dir']),
        ld_blocks = '{out_dir}/haploview/ld_blocks/{{roi_name}}.{hv_method}/'.format(out_dir=config['out_dir'],
                                                                    hv_method = config['haploview_method'])

