#this is a python script loosely based on Kumar and Konrad's effort here: https://github.com/mkveerapen/covid19_sequencing

import hail as hl
import hail.expr.aggregators as agg
from bokeh.io import export_png #this requires selenium (conda install -c bokeh selenium) and phantomjs (conda install phantomjs)
import time

hl.init(spark_conf=None, log='02.ion_hail_JJ.log', tmp_dir='tmp/')

pathNorm = 'ion.normID.noChrM.vqsr.flt.vcf.gz'
pathMT   = 'ion.normID.noChrM.vqsr.flt.vcf.mt'
pathQC   = 'ion.normID.GTflt.AB.noChrM.vqsr.flt.vcf.bgz'

hl.import_vcf(pathNorm, min_partitions=4, reference_genome='GRCh38', force_bgz=True).write(pathMT, overwrite=True)

metaData = hl.get_vcf_metadata(pathNorm)
mt = hl.read_matrix_table(pathMT)
#mt.count()

print("#####################################################################################")
print(" 2.0 Sample QC (Sex Imputation filtering / additional filtering / outlier filtering)")
print("#####################################################################################")
### 2.0.1 call_rate filter
mt = hl.sample_qc(mt)
mt = mt.filter_cols(mt.sample_qc.call_rate >= 0.96)
n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])

print("#####################################################################################")
print(" 2.0.2 Sex Imputation")
print("#####################################################################################")
### First filter for high quality calls for sex QC
mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                            (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
                            (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.97))
imputed_sex = hl.impute_sex(mt.GT,aaf_threshold=0.05, female_threshold=0.5, male_threshold=0.75)
mt = mt.annotate_cols(impute_sex = imputed_sex[mt.s])

print("#####################################################################################")
print(" 2.0.3 Additional filters")
print("#####################################################################################")
mt = mt.annotate_cols(aneuploid= ((mt.impute_sex.f_stat >= 0.5) ) | (hl.is_missing(mt.impute_sex.f_stat)) | 
                      ((mt.impute_sex.f_stat >= 0.4) & (mt.impute_sex.f_stat <= 0.6)) ,
        sex_aneuploidy=(mt.impute_sex.f_stat < 0.4))

mt = mt.filter_cols( (mt.sample_qc.call_rate >= 0.96) &
                    (mt.sample_qc.dp_stats.mean > 20) & (hl.is_defined(mt.aneuploid))  )
n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])


print("#####################################################################################")
print(" 2.0.4 Relatedness filter")
print("#####################################################################################")
pca_eigenvalues, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=False)
relatedness_ht = hl.pc_relate(mt.GT, min_individual_maf=0.01, scores_expr=pca_scores[mt.col_key].scores,
                                      block_size=4096, min_kinship=0.1, statistics='all')
pairs = relatedness_ht.filter(relatedness_ht['kin'] > 0.088)
related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)
n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])

print("#####################################################################################")
print("2.0.6 Outlier Detection")
print("#####################################################################################")
metric_values = hl.agg.collect(mt.sample_qc.n_snp)
metric_median = hl.median(metric_values)
metric_mad = 1.4826 * hl.median(hl.abs(metric_values - metric_median))
outlier_metric=hl.struct( median=metric_median,
            mad=metric_mad,
            upper=metric_median + 4 * metric_mad,
            lower=metric_median - 4 * metric_mad)

mt = mt.annotate_globals(metrics_stats=mt.aggregate_cols(outlier_metric))
mt = mt.filter_cols((mt.sample_qc.n_snp <= mt.metrics_stats.upper) | 
                    (mt.sample_qc.n_snp >=  mt.metrics_stats.lower))
n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])

# 2.1 Variant Quality Control
#mt = mt.filter_rows(hl.len(mt.filters) == 0)
#mt.count()


print("#####################################################################################")
print(" 2.2 Genotype Quality Control")
mt = mt.annotate_entries(AB = (mt.AD[1] / hl.sum(mt.AD) ))
print("#####################################################################################")
mt = mt.filter_entries((mt.GQ>=20) & (mt.GQ >= 10) & 
                       ((mt.GT.is_hom_ref() & (mt.AB <= 0.1)) | 
                        (mt.GT.is_het() & (mt.AB >= 0.25) & (mt.AB <= 0.75)) | 
                        (mt.GT.is_hom_var() & (mt.AB >= 0.9))))
#mt.count()
mt = hl.variant_qc(mt)
mt = mt.filter_rows((mt.variant_qc.gq_stats.mean > 10) & (mt.variant_qc.dp_stats.mean > 5) & (mt.variant_qc.call_rate > 0.8)  & (mt.variant_qc.p_value_hwe>5e-8))
n = mt.count()

pprint('n samples:')
print(n[1])
pprint('n variants:')
print(n[0])


#plot different QC metrics
p_sample_call=hl.plot.histogram(mt.sample_qc.call_rate,legend='Sample Call Rate')
p_sample_dp=hl.plot.histogram(mt.sample_qc.dp_stats.mean,legend='Sample Mean DP')
p_variant_call=hl.plot.histogram(mt.variant_qc.call_rate,legend='Variant Call Rate')
p_mean_variant_dp=hl.plot.histogram(mt.variant_qc.dp_stats.mean,legend='Variant Mean DP')
p_variant_hwe_p=hl.plot.histogram(mt.variant_qc.p_value_hwe,legend='Variant HWE P-values')
export_png(p_sample_call, filename="fig_ion_hail_newcovid19/p_sample_call.png")
export_png(p_sample_dp, filename="fig_ion_hail_newcovid19/p_sample_dp.png")
export_png(p_variant_call, filename="fig_ion_hail_newcovid19/p_variant_call.png")
export_png(p_mean_variant_dp, filename="fig_ion_hail_newcovid19/p_mean_variant_dp.png")
export_png(p_variant_hwe_p, filename="fig_ion_hail_newcovid19/p_variant_hwe_p.png")

#export back to vcf          
hl.export_vcf(mt, pathQC, metadata=metaData, tabix=True)
