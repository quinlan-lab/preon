params.project = false
params.fasta = false
params.manifest = false
params.width = 1000000
params.chunksize = 30000000
params.caddsnv = false
params.caddindel = false
params.outdir = './results'
params.intervals = false

intervals = false
if (params.intervals) {
    intervals = file(params.intervals)
}

fasta = file(params.fasta)
faidx = file("${params.fasta}.fai")
exclude = params.exclude.tokenize(',')
outdir = file(params.outdir)

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
    tag { sample_id + "_" + run_id }

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
    tag "$sample_id"
    publishDir path: "$outdir/alignments"

    input:
    set sample_id, file(bam) from bwa_ch.groupTuple()

    output:
    file("${sample_id}.md.bam") into alignments_ch
    file("${sample_id}.md.bam.bai") into alignment_indexes_ch

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
    publishDir path: "${params.outdir}/freebayes"

    input:
    file(vcf) from unannotated_ch.collect()
    file(fasta)
    file(faidx)

    output:
    file("${params.project}.vcf.gz")
    file("${params.project}.vcf.gz.tbi")

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
        --cache --cache_version 98 --dir_cache /.vep --fork 2 $cadd --offline \
        --everything --filter_common --per_gene --total_length --vcf --force_overwrite \
        --buffer_size 15000
    gzip ${vcf.baseName}.vep.vcf
    """
}


process merge_annotated_vcfs {
    publishDir path: "${params.outdir}/vep"

    input:
    file(vcf) from annotated_ch.collect()
    file(fasta)
    file(faidx)

    output:
    file("${params.project}.vcf.gz")
    file("${params.project}.vcf.gz.tbi")

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
