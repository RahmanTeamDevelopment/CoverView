
def region_size(region_name, region_coords):
    x = region_coords[region_name]
    return int(x[2])-int(x[1])

def create_flag_statistics(regions_data, region_coords):

    by_gene = {}
    counter_flagged_regions = 0
    for r in regions_data:
        region_name = r['region']

        if '_' in region_name:
            gene_name = region_name.split('_')[0]
        else:
            gene_name = region_name

        flag = (r['pass_or_flag'] == 'FLAG')
        if flag:
            counter_flagged_regions += 1

        '''
        if region_size(region_name, region_coords) == 1:
            continue
        '''

        if gene_name not in by_gene:
            by_gene[gene_name] = {
                'exons': 0,
                'flaggedexons': 0,
                'chrom': region_coords[region_name][0]
            }

        by_gene[gene_name]['exons'] += 1
        if flag:
            by_gene[gene_name]['flaggedexons'] += 1

    all_genes = []
    by_chrom = {}
    counter_flagged_genes = 0
    for key, value in by_gene.iteritems():
        chrom = value['chrom']
        if chrom not in by_chrom:
            by_chrom[chrom] = []
        g = {
            'gene': key,
            'exons': value['exons'],
            'flaggedexons': value['flaggedexons']
        }
        all_genes.append(g)
        by_chrom[chrom].append(g)
        if value['flaggedexons'] > 0:
            counter_flagged_genes += 1

    return all_genes, by_chrom, '{}  [of {} regions]'.format(counter_flagged_regions, len(regions_data)), '{}  [of {} genes]'.format(counter_flagged_genes, len(by_gene))
