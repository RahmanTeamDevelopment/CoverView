
def region_size(region_name, region_coords):
    x = region_coords[region_name]
    return int(x[2])-int(x[1])

def create_fail_statistics(regions_data, region_coords):

    by_gene = {}
    counter_failed_regions = 0
    for r in regions_data:
        region_name = r['region']

        if '_' in region_name:
            gene_name = region_name.split('_')[0]
        else:
            gene_name = region_name

        fail = (r['pass_or_fail'] == 'FAIL')
        if fail:
            counter_failed_regions += 1

        if region_size(region_name, region_coords) == 1:
            continue

        if gene_name not in by_gene:
            by_gene[gene_name] = {
                'exons': 0,
                'failedexons': 0,
                'chrom': region_coords[region_name][0]
            }

        by_gene[gene_name]['exons'] += 1
        if fail:
            by_gene[gene_name]['failedexons'] += 1

    all_genes = []
    by_chrom = {}
    counter_failed_genes = 0
    for key, value in by_gene.iteritems():
        chrom = value['chrom']
        if chrom not in by_chrom:
            by_chrom[chrom] = []
        g = {
            'gene': key,
            'exons': value['exons'],
            'failedexons': value['failedexons']
        }
        all_genes.append(g)
        by_chrom[chrom].append(g)
        if value['failedexons'] > 0:
            counter_failed_genes += 1

    return all_genes, by_chrom, '{}  [of {} regions]'.format(counter_failed_regions, len(regions_data)), '{}  [of {} genes]'.format(counter_failed_genes, len(by_gene))
