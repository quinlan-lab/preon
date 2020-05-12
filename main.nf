params.project = false
params.fasta = false
params.manifest = false
params.width = 1000000
params.chunksize = 30000000
params.caddsnv = false
params.caddindel = false
params.outdir = './results'
params.intervals = false
params.sbgproject = false
params.bed = false
params.gff = false
params.upload = false
params.authtoken = "$AUTH_TOKEN"

// somatic VCF
params.minquality = 20
params.mindepth = 15
params.maxao = 2
params.somaticmask = false

if (params.authtoken == "") {
    // exit 1, "--authtoken for SevenBridges must be defined for data upload. See: https://docs.sevenbridges.com/docs/get-your-authentication-token"
    log.info("Uploading to SevenBridges will be skipped without --authtoken or empty \$AUTH_TOKEN")
}

if (params.authtoken && !params.sbgproject) {
    exit 1, "--sbgproject needs to specify the SevenBridges project when specifying --authtoken or \$AUTH_TOKEN"
}

intervals = false
if (params.intervals) {
    intervals = file(params.intervals)
}

fasta = file(params.fasta)
faidx = file("${params.fasta}.fai")
exclude = params.exclude.tokenize(',')
outdir = file(params.outdir)
bed = file(params.bed)
gff = file(params.gff)
sexchroms = params.sexchroms ?: 'X,Y'
sexchroms = sexchroms.replaceAll(" ", "")

// CADD references for VEP
cadd_snv = file(params.caddsnv) ?: false
cadd_snv_idx = cadd_snv ? file("${params.caddsnv}.tbi") : false
cadd_indel = file(params.caddindel) ?: false
cadd_indel_idx = cadd_indel ? file("${params.caddindel}.tbi") : false


// ingests the sample/file manifest
// format of the TSV is: patient, sample, status, lane, r1, r2
Channel.fromPath(params.manifest)
    .splitCsv(sep: '\t')
    .map { row ->
        def sample_id = row[1]
        def status = row[2].toInteger()
        def run_id = row[3]
        def r1 = file(row[4])
        def r2 = file(row[5])
        [[sample_id, status, run_id], r1, r2]
    }
    // expects files as last 2 elements
    .splitFastq(by: params.chunksize, pe:true, file:true)
    .map { row ->
        def identifiers = row[0]
        def sample_id = identifiers[0]
        def status = identifiers[1]
        def run_id = identifiers[2]
        def idx = row[1].baseName.split(/\\.fastq|\\.fq/)[0].split("\\.")[-1]
        [sample_id, status, run_id, idx, row[1], row[2]]
    }
    .set { fastq_ch }

// grab the indexes for bwa
Channel
    .value(file("${params.fasta}.{amb,ann,bwt,pac,sa}"))
    .set { bwaidx_ch }

// creates intervals across the genome for freebayes and vep
intervals_ch = Channel
    .from(params.intervals ? intervals : faidx)
    .splitCsv(sep: '\t')
    .map { row ->
        // rows of interval lists
        if (row[0][0] != "@") {
            def interval_start = row[1].toLong()
            def interval_length = row[2].toLong()
            long start
            long end
            int width = params.width

            if (!params.intervals) {
                // update interval start and length for .fai
                interval_start = 0
                interval_length = row[1].toLong()
            }

            while(interval_start < interval_length) {
                start = interval_start
                // add a slight overlap
                end = interval_start + width + 500
                interval_start = end - 500
                if (end > interval_length) {
                    end = interval_length
                    interval_start = end
                }
                // add the interval to the channel
                intervals_ch.bind( "${row[0]}:${start}-${end}" )
            }
        }
    }

// not currently in use
// grabs chromosome IDs from the given reference
// Channel
//     .fromPath("${params.fasta}.fai")
//     .splitCsv(sep: "\t", strip: true)
//     .map { row -> "${row[0]}" }
//     .filter( ~/(?!${exclude.collect {".*$it.*"}.join("|")})([a-zA-Z0-9_]+)/ )
//     .into { chr_ch; gather_ch }


process map_reads {
    tag { "${sample_id}_${run_id}" }

    input:
    set sample_id, status, run_id, idx, file(r1), file(r2) from fastq_ch
    file(fasta)
    file(bwaidx) from bwaidx_ch

    output:
    set sample_id, file("${sample_id}_${run_id}_${idx}.bam") into bwa_ch

    script:
    rg = "@RG\\tID:${run_id}\\tPU:${run_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina"
    """
    bwa mem -K 100000000 -R \"${rg}\" -t ${task.cpus} -M $fasta $r1 $r2 \
        | samtools sort -n --threads ${task.cpus} -m 2G --output-fmt BAM -o ${sample_id}_${run_id}_${idx}.bam
    """
}


process mark_duplicates {
    tag "${sample_id}"
    publishDir path: "$outdir/alignments"

    input:
    set sample_id, file(bam) from bwa_ch.groupTuple()

    output:
    file("${sample_id}.md.bam") into (alignments_ch, alignments_upload_ch, create_call_bams_ch)
    file("${sample_id}.md.bam.bai") into (alignment_indexes_ch, alignment_indexes_upload_ch, indexcov_input_ch)

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g MarkDuplicatesSpark \
        ${bam.collect { "--input $it" }.join(" ")} \
        --output ${sample_id}.md.bam \
        --tmp-dir . \
        --spark-master \'local[*]\'
    gatk --java-options -Xmx${task.memory.toGiga()}g BuildBamIndex \
        --INPUT ${sample_id}.md.bam \
        --OUTPUT ${sample_id}.md.bam.bai \
        --TMP_DIR .
    """
}


process upload_alignments {
    label 'upload'

    input:
    val token from params.authtoken
    val sbgproject from params.sbgproject
    file(bam) from alignments_upload_ch
    file(bai) from alignment_indexes_upload_ch

    when:
    token != ""

    script:
    """
    sbg-uploader.sh -t $token -p $sbgproject -f WGS/bam/ $bam $bai
    """
}


process run_freebayes {
    tag "$interval"

    input:
    file(aln) from alignments_ch.collect()
    file(idx) from alignment_indexes_ch.collect()
    each interval from intervals_ch
    file(fasta)
    file(faidx)

    output:
    // need to rename due to colon in intervals
    file("${params.project}_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz") into (unannotated_ch, vcf_ch)

    script:
    """
    freebayes \
        --fasta-reference ${fasta} \
        ${params.freebayesoptions} \
        --region ${interval} \
        ${aln.collect { "--bam $it" }.join(" ")} \
        | bgzip -c > ${params.project}_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz
    """
}


process merge_vcfs {
    publishDir path: "$outdir/freebayes"

    input:
    file(vcf) from unannotated_ch.collect()
    file(fasta)
    file(faidx)

    output:
    file("${params.project}.vcf.gz") into freebayesvcf_ch
    file("${params.project}.vcf.gz.tbi") into freebayesvcfidx_ch

    script:
    """
    gunzip -cd $vcf | vcffirstheader | bgzip -c > ${params.project}_dirty.vcf.gz
    gsort ${params.project}_dirty.vcf.gz $faidx | vcfuniq \
        | bgzip -c > ${params.project}_dirty_sorted.vcf.gz
    bcftools norm -c all -f $fasta --multiallelics - --threads ${task.cpus} \
        --output ${params.project}.vcf.gz --output-type z \
        ${params.project}_dirty_sorted.vcf.gz
    tabix -p vcf ${params.project}.vcf.gz
    """
}


process upload_vcfs {
    label 'upload'

    input:
    val token from params.authtoken
    val sbgproject from params.sbgproject
    file(vcf) from freebayesvcf_ch
    file(tbi) from freebayesvcfidx_ch

    when:
    token != ""

    script:
    """
    sbg-uploader.sh -t $token -p $sbgproject -f WGS/freebayes/ $vcf $tbi
    """
}


process annotate_vcf {
    input:
    file(vcf) from vcf_ch
    file(cadd_indel)
    file(cadd_indel_idx)
    file(cadd_snv)
    file(cadd_snv_idx)

    output:
    file("${vcf.baseName}.vep.vcf.gz") into annotated_ch

    script:
    cadd = (params.caddsnv && params.caddindel) ? "--plugin CADD,${cadd_snv},${cadd_indel}" : ""
    """
    vep --input_file $vcf --output_file ${vcf.baseName}.vep.vcf --format vcf --assembly GRCh38 \
        --cache --cache_version 95 --dir_cache /.vep --fork 2 $cadd --offline \
        --everything --filter_common --per_gene --total_length --vcf --force_overwrite \
        --buffer_size 15000
    gzip ${vcf.baseName}.vep.vcf
    """
}


process merge_annotated_vcfs {
    publishDir path: "$outdir/vep"

    input:
    file(vcf) from annotated_ch.collect()
    file(fasta)
    file(faidx)

    output:
    file("${params.project}.vcf.gz") into {vepvcf_ch; somaticvcf_ch}
    file("${params.project}.vcf.gz.tbi") into vepvcfidx_ch

    script:
    """
    gunzip -cd $vcf | vcffirstheader | bgzip -c > ${params.project}_dirty.vcf.gz
    gsort ${params.project}_dirty.vcf.gz $faidx | vcfuniq \
        | bgzip -c > ${params.project}_dirty_sorted.vcf.gz
    bcftools norm -c all -f $fasta --multiallelics - --threads ${task.cpus} \
        --output ${params.project}.vcf.gz --output-type z \
        ${params.project}_dirty_sorted.vcf.gz
    tabix -p vcf ${params.project}.vcf.gz
    """
}


process upload_annotated_vcfs {
    label 'upload'

    input:
    val token from params.authtoken
    val sbgproject from params.sbgproject
    file(vcf) from vepvcf_ch
    file(tbi) from vepvcfidx_ch

    when:
    token != ""

    script:
    """
    sbg-uploader.sh -t $token -p $sbgproject -f WGS/vep/ $vcf $tbi
    """
}


// create channels for SV calling
create_call_bams_ch
    .map { file -> tuple(file.getSimpleName(), file, file + ("${file}".endsWith('.cram') ? '.crai' : '.bai')) }
    .into { call_bams; genotype_bams }


process smoove_call {
    publishDir path: "$outdir/sv/smoove/called", pattern: "*.vcf.gz*"
    publishDir path: "$outdir/sv/logs", pattern: "*-stats.txt"
    publishDir path: "$outdir/sv/logs", pattern: "*-smoove-call.log"

    input:
    env SMOOVE_KEEP_ALL from params.sensitive
    set sample, file(bam), file(bai) from call_bams
    file fasta
    file faidx
    file bed

    output:
    file("${sample}-smoove.genotyped.vcf.gz") into call_vcfs
    file("${sample}-smoove.genotyped.vcf.gz.csi") into call_idxs
    file("${sample}-stats.txt") into variant_counts
    file("${sample}-smoove-call.log") into sequence_counts

    script:
    excludechroms = params.exclude ? "--excludechroms \"${params.exclude}\"" : ""
    filters = params.sensitive ? "--noextrafilters" : ""
    """
    smoove call --genotype --name $sample --processes ${task.cpus} \
        --fasta $fasta --exclude $bed $excludechroms $filters \
        $bam 2> >(tee -a ${sample}-smoove-call.log >&2)
    bcftools stats ${sample}-smoove.genotyped.vcf.gz > ${sample}-stats.txt
    """
}


process smoove_merge {
    input:
    file vcf from call_vcfs.collect()
    file idx from call_idxs.collect()
    file fasta
    file faidx

    output:
    file("${params.project}.sites.vcf.gz") into sites

    script:
    """
    smoove merge --name ${params.project} --fasta $fasta $vcf
    """
}


process smoove_genotype {
    publishDir path: "$outdir/sv/smoove/genotyped"

    input:
    env SMOOVE_KEEP_ALL from params.sensitive
    set sample, file(bam), file(bai) from genotype_bams
    file sites
    file fasta
    file faidx

    output:
    file("${sample}-smoove.genotyped.vcf.gz.csi") into genotyped_idxs
    file("${sample}-smoove.genotyped.vcf.gz") into genotyped_vcfs

    script:
    """
    wget -q https://raw.githubusercontent.com/samtools/samtools/develop/misc/seq_cache_populate.pl
    perl seq_cache_populate.pl -root \$(pwd)/cache $fasta 1> /dev/null 2> err || (cat err; exit 2)
    export REF_PATH=\$(pwd)/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
    export REF_CACHE=xx

    samtools quickcheck -v $bam
    smoove genotype --duphold --processes ${task.cpus} --removepr --outdir ./ --name ${sample} --fasta $fasta --vcf $sites $bam
    """
}


process smoove_square {
    publishDir path: "$outdir/sv/smoove/annotated", pattern: "*.vcf.gz*"
    publishDir path: "$outdir/reports/bpbio", pattern: "*.html"

    input:
    file vcf from genotyped_vcfs.collect()
    file idx from genotyped_idxs.collect()
    file gff

    output:
    file("${params.project}.smoove.square.anno.vcf.gz") into square_vcf
    file("${params.project}.smoove.square.anno.vcf.gz.csi") into square_idx
    file("svvcf.html") into svvcf

    script:
    smoovepaste = "smoove paste --outdir ./ --name ${params.project} $vcf"
    if( vcf.collect().size() < 2 ) {
        paste = "cp $vcf ${params.project}.smoove.square.vcf.gz && cp $idx ${params.project}.smoove.square.vcf.gz.csi"
    }
    """
    $smoovepaste

    smoove annotate --gff $gff ${params.project}.smoove.square.vcf.gz | bgzip --threads ${task.cpus} -c > ${params.project}.smoove.square.anno.vcf.gz
    bcftools index ${params.project}.smoove.square.anno.vcf.gz
    bpbio plot-sv-vcf ${params.project}.smoove.square.anno.vcf.gz
    """
}


process run_indexcov {
    publishDir path: "$outdir/reports/indexcov"

    input:
    file idx from indexcov_input_ch.collect()
    file faidx

    output:
    file("${params.project}*.png")
    file("*.html")
    file("${params.project}*.bed.gz") into bed_ch
    file("${params.project}*.ped") into indexcov_ped_ch
    file("${params.project}*.roc") into roc_ch

    script:
    excludepatt = params.exclude ? "--excludepatt \"${params.exclude}\"" : ""
    """
    goleft indexcov --sex $sexchroms $excludepatt --directory ${params.project} --fai $faidx $idx
    mv ${params.project}/* .
    """
}


process build_covviz_report {
    publishDir path: "$outdir/reports", mode: "copy", pattern: "*.html"
    label 'covviz'
    cache 'lenient'

    input:
    file ped from indexcov_ped_ch
    file bed from bed_ch
    file gff

    output:
    file("covviz_report.html")

    script:
    """
    covviz --min-samples ${params.minsamples} --sex-chroms $sexchroms --exclude '${params.exclude}' \
        --skip-norm --z-threshold ${params.zthreshold} --distance-threshold ${params.distancethreshold} \
        --slop ${params.slop} --ped $ped --gff $gff $bed
    """
}


process make_somatic {
    publishDir path:  "$outdir/vep", mode: "copy"

    input:
    file vcf from somaticvcf_ch
    file mask from params.somaticmask

    output:
    file("somatic.vcf")

    script:
    // entirely serial as each is dependent upon previous
    // need to test failure rate of individual steps
    """
    bedtools intersect -a $vcf -b $mask -wa -u -header > masked.vcf
    SnpSift filter "(QUAL >= ${params.minquality}) & (GEN[ALL].DP >= ${params.mindepth})" masked.vcf > filtered.vcf
    bcftools view --min-alleles 2 --max-alleles 2 --output-file biallelic.vcf filtered.vcf
    snpEff -Xmx16g -i vcf -o vcf -noStats GRCh38.86 biallelic.vcf > snpeff.vcf
    SnpSift filter "(GEN[0].AO<=${params.maxao})" snpeff.vcf | bgzip > somatic.vcf.gz
    tabix -p vcf somatic.vcf.gz
    """
}
