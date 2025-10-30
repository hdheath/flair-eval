#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// CLASS DEFINITIONS - Define dataset and test configuration objects
// =============================================================================

class Dataset {
    String name
    def reads  // Can be String (single file) or List (multiple files)
    String genome
    String gtf
    String junctions
    String cage
    String quantseq
    
    // Constructor
    Dataset(String name, Map config) {
        this.name = name
        this.reads = config.reads
        this.genome = config.genome
        this.gtf = config.gtf
        this.junctions = config.junctions
        this.cage = config.cage
        this.quantseq = config.quantseq
    }
    
    // Helper methods
    boolean hasJunctions() { return junctions != null }
    boolean hasCage() { return cage != null }
    boolean hasQuantseq() { return quantseq != null }
    
    // Get reads as a list (handles both single file and multiple files)
    List<String> getReadsList() {
        return reads instanceof List ? reads : [reads]
    }
    
    String description() {
        return "Dataset[${name}]"
    }
}

class TestSet {
    String name
    Dataset dataset
    List<List> alignOptions
    List<List> transcriptomeOptions
    List<List> partitionOptions
    List<List> correctOptions
    List<List> collapseOptions
    
    // Constructor
    TestSet(String name, Dataset dataset, Map options) {
        this.name = name
        this.dataset = dataset
        this.alignOptions = options.align
        this.transcriptomeOptions = options.transcriptome ?: []
        this.partitionOptions = options.partition
        this.correctOptions = options.correct
        this.collapseOptions = options.collapse
    }
    
    // Calculate total number of jobs this test set will produce
    int totalJobs() {
        return alignOptions.size() * 
               partitionOptions.size() * 
               correctOptions.size() * 
               collapseOptions.size()
    }
    
    String description() {
        return "TestSet[${name}, ${totalJobs()} jobs]"
    }
}

// =============================================================================
// PROCESS DEFINITIONS
// =============================================================================

process FlairAlign {
    publishDir "results/align", mode: 'copy'
    errorStrategy 'ignore'  // Allow other test sets to continue if this fails
    
    input:
    tuple val(test_name), val(dataset_name), path(reads), val(align_label), val(align_args)
    path genome
    
    output:
    tuple val(test_name), val(dataset_name), val(align_label), 
          path("${dataset_name}_${align_label}.bam"), 
          path("${dataset_name}_${align_label}.bam.bai"),
          path("${dataset_name}_${align_label}.bed"), emit: alignments
    
    script:
    """
    flair align ${align_args} \\
        -r ${reads} \\
        -g ${genome} \\
        -o ${dataset_name}_${align_label}
    """
}

process FlairTranscriptome {
    publishDir "results/transcriptome/${test_name}", mode: 'copy'
    errorStrategy 'ignore'  // Allow correct/collapse branch to continue if transcriptome fails
    
    input:
    tuple val(test_name), val(dataset_name), val(align_label), path(bam), path(bai), path(bed), val(transcriptome_args)
    path genome
    
    output:
    tuple val(test_name), val(dataset_name), val(align_label),
          path("${dataset_name}_${align_label}.transcriptome.isoforms.bed"),
          path("${dataset_name}_${align_label}.transcriptome.isoforms.gtf"),
          path("${dataset_name}_${align_label}.transcriptome.isoforms.fa"),
          path("${dataset_name}_${align_label}.transcriptome.isoform.counts.txt"),
          path("${dataset_name}_${align_label}.transcriptome.isoform.read.map.txt"), emit: transcriptome
    
    script:
    """
    flair transcriptome \\
        -b ${bam} \\
        --genome ${genome} \\
        ${transcriptome_args} \\
        -o ${dataset_name}_${align_label}.transcriptome
    """
}

process FlairPartition {
    publishDir "results/partition/${test_name}", mode: 'copy'
    
    input:
    tuple val(test_name), val(dataset_name), val(align_label), path(bam), path(bai), path(bed), val(partition_label), val(partition_args)
    path genome
    path gtf
    path cage_peaks
    path quantseq_peaks
    path junctions
    path target_regions
    
    output:
    tuple val(test_name), val(dataset_name), val(align_label), val(partition_label), 
          path("${dataset_name}_${align_label}_${partition_label}.bam"), 
          path("${dataset_name}_${align_label}_${partition_label}.bam.bai"),
          path("${dataset_name}_${align_label}_${partition_label}.bed"),
          path("${dataset_name}_${align_label}_${partition_label}_genome.fa"),
          path("${dataset_name}_${align_label}_${partition_label}_annotation.gtf"),
          path("${dataset_name}_${align_label}_${partition_label}_cage.bed"),
          path("${dataset_name}_${align_label}_${partition_label}_quantseq.bed"),
          path("${dataset_name}_${align_label}_${partition_label}_junctions.bed"),
          path("${dataset_name}_${align_label}_${partition_label}_targets.bed"), emit: partitioned
    
    script:
    def output_prefix = "${dataset_name}_${align_label}_${partition_label}"
    
    // Check if this is a passthrough (--all) or region slice
    if (partition_args.contains('--all')) {
        """
        # Passthrough mode - symlink original files
        ln -s ${bam} ${output_prefix}.bam
        ln -s ${bai} ${output_prefix}.bam.bai
        ln -s ${bed} ${output_prefix}.bed
        ln -s ${genome} ${output_prefix}_genome.fa
        ln -s ${gtf} ${output_prefix}_annotation.gtf
        
        # Handle optional files
        if [ "${cage_peaks.name}" != "NO_CAGE" ]; then
            ln -s ${cage_peaks} ${output_prefix}_cage.bed
        else
            touch ${output_prefix}_cage.bed
        fi
        
        if [ "${quantseq_peaks.name}" != "NO_QUANTSEQ" ]; then
            ln -s ${quantseq_peaks} ${output_prefix}_quantseq.bed
        else
            touch ${output_prefix}_quantseq.bed
        fi
        
        if [ "${junctions.name}" != "NO_JUNCTIONS" ]; then
            ln -s ${junctions} ${output_prefix}_junctions.bed
        else
            touch ${output_prefix}_junctions.bed
        fi
        
        if [ "${target_regions.name}" != "NO_TARGETS" ]; then
            ln -s ${target_regions} ${output_prefix}_targets.bed
        else
            touch ${output_prefix}_targets.bed
        fi
        """
    } else {
        // Extract region from partition_args (format: --region chr1:1-100000)
        """
        # Parse region from arguments
        REGION=\$(echo "${partition_args}" | sed -n 's/.*--region *\\([^ ]*\\).*/\\1/p')
        CHR=\$(echo "\$REGION" | cut -d: -f1)
        START=\$(echo "\$REGION" | cut -d: -f2 | cut -d- -f1)
        END=\$(echo "\$REGION" | cut -d: -f2 | cut -d- -f2)
        
        echo "Partitioning to region: \$CHR:\$START-\$END"
        
        # Slice BAM (input BAM is already indexed from flair align)
        samtools view -b ${bam} "\$REGION" | samtools sort -o ${output_prefix}.bam
        samtools index ${output_prefix}.bam
        
        # Slice BED (0-based coordinates)
        awk -v chr="\$CHR" -v start="\$START" -v end="\$END" \\
            '\$1==chr && \$2>=start && \$3<=end' ${bed} > ${output_prefix}.bed
        
        # Slice GTF (1-based coordinates, keep comments)
        awk -v chr="\$CHR" -v start="\$START" -v end="\$END" \\
            'BEGIN {OFS="\\t"} /^#/ {print; next} \$1==chr && \$4>=start && \$5<=end' \\
            ${gtf} > ${output_prefix}_annotation.gtf
        
        # Slice genome FASTA
        samtools faidx ${genome} "\$REGION" > ${output_prefix}_genome.fa || touch ${output_prefix}_genome.fa
        
        # Slice optional files if they exist
        if [ "${cage_peaks.name}" != "NO_CAGE" ]; then
            awk -v chr="\$CHR" -v start="\$START" -v end="\$END" \\
                '\$1==chr && \$2>=start && \$3<=end' ${cage_peaks} > ${output_prefix}_cage.bed
        else
            touch ${output_prefix}_cage.bed
        fi
        
        if [ "${quantseq_peaks.name}" != "NO_QUANTSEQ" ]; then
            awk -v chr="\$CHR" -v start="\$START" -v end="\$END" \\
                '\$1==chr && \$2>=start && \$3<=end' ${quantseq_peaks} > ${output_prefix}_quantseq.bed
        else
            touch ${output_prefix}_quantseq.bed
        fi
        
        if [ "${junctions.name}" != "NO_JUNCTIONS" ]; then
            awk -v chr="\$CHR" -v start="\$START" -v end="\$END" \\
                '\$1==chr && \$2>=start && \$3<=end' ${junctions} > ${output_prefix}_junctions.bed
        else
            touch ${output_prefix}_junctions.bed
        fi
        
        if [ "${target_regions.name}" != "NO_TARGETS" ]; then
            awk -v chr="\$CHR" -v start="\$START" -v end="\$END" \\
                '\$1==chr && \$2>=start && \$3<=end' ${target_regions} > ${output_prefix}_targets.bed
        else
            touch ${output_prefix}_targets.bed
        fi
        """
    }
}

process FlairCorrect {
    publishDir "results/correct/${test_name}", mode: 'copy', saveAs: { filename -> 
        task.exitStatus == 0 ? filename : null  // Only publish successful outputs
    }
    errorStrategy 'ignore'  // Allow other correct options and transcriptome to continue
    
    input:
    tuple val(test_name), val(dataset_name), val(align_label), val(partition_label), 
          path(bam), path(bed), path(genome), path(gtf),
          path(cage), path(quantseq), path(junctions_bed), path(targets),
          val(correct_label), val(correct_args)
    
    output:
    tuple val(test_name), val(dataset_name), val(align_label), val(partition_label), val(correct_label),
          path("${dataset_name}_${align_label}_${partition_label}_${correct_label}_all_corrected.bed", optional: true),
          path("${dataset_name}_${align_label}_${partition_label}_${correct_label}_all_inconsistent.bed", optional: true), emit: corrected
    
    script:
    """
    flair correct ${correct_args} \\
        -q ${bed} \\
        -f ${gtf} \\
        -o ${dataset_name}_${align_label}_${partition_label}_${correct_label}
    """
}

process FlairCollapse {
    publishDir "results/collapse/${test_name}", mode: 'copy'
    errorStrategy 'ignore'  // Allow other test sets to continue if this fails
    
    input:
    tuple val(test_name), val(dataset_name), val(align_label), val(partition_label), val(correct_label),
          path(corrected_bed), path(inconsistent_bed), path(reads),
          path(genome), path(gtf), val(collapse_label), val(collapse_args)
    
    output:
    tuple val(test_name), val(dataset_name), val(align_label), val(partition_label), val(correct_label), val(collapse_label),
          path("${dataset_name}_${align_label}_${partition_label}_${correct_label}_${collapse_label}.isoforms.bed"),
          path("${dataset_name}_${align_label}_${partition_label}_${correct_label}_${collapse_label}.isoforms.fa"),
          path("${dataset_name}_${align_label}_${partition_label}_${correct_label}_${collapse_label}.isoforms.gtf"),
          path("${dataset_name}_${align_label}_${partition_label}_${correct_label}_${collapse_label}.isoform.read.map.txt"),
          path("${dataset_name}_${align_label}_${partition_label}_${correct_label}_${collapse_label}.isoform.counts.txt"), emit: collapsed
    
    script:
    """
    flair collapse ${collapse_args} \\
        -r ${reads} \\
        -q ${corrected_bed} \\
        -g ${genome} \\
        -f ${gtf} \\
        --generate_map \\
        -o ${dataset_name}_${align_label}_${partition_label}_${correct_label}_${collapse_label}
    """
}

process FlairEval {
    publishDir "results/eval/${test_name}", mode: 'copy'
    errorStrategy 'ignore'  // Allow other test sets to continue if this fails
    
    input:
    tuple val(test_name), val(dataset_name), val(align_label), val(partition_label), val(process_label),
          path(isoforms_bed), path(reads_bed), path(gtf)
    
    output:
    tuple val(test_name), val(dataset_name), val(align_label), val(partition_label), val(process_label),
          path("${dataset_name}_${align_label}_${partition_label}_${process_label}.eval_summary.txt"), emit: eval_results
    
    script:
    """
    #!/usr/bin/env python3
    import sys
    
    SINGLE_EXON_END_WINDOW = 100
    
    BED_READS_FILE = '${reads_bed}'
    BED_READS_INTERVALS_FILE = '${dataset_name}_${align_label}_${partition_label}_${process_label}.reads.intervals.bed'
    ANNOT_FILE = '${gtf}'
    ISOFORMS_FILE = '${isoforms_bed}'
    
    def get_chromtoint(file):
        totreads = 0
        chromtoint = {}
        for line in open(file): 
            line = line.rstrip().split('\\t')
            chrom, strand = line[0], line[5]
            iso, start, end = line[3], int(line[1]), int(line[2])
            if chrom not in chromtoint: chromtoint[chrom] = []
            chromtoint[chrom].append((start, end))
            totreads += 1
        return chromtoint, totreads
    
    def get_regions(chromtoint, outfilename):
        totregions= 0
        out = open(outfilename, 'w')
        for chrom in chromtoint:
            intlist = sorted(chromtoint[chrom])
            newints = []
            laststart, lastend = 0, 0
            for s, e in intlist:
                if s > lastend:
                    if lastend != 0:
                        newints.append((laststart, lastend))
                    laststart = s
                if e > lastend:
                    lastend = e
            newints.append((laststart, lastend))
            for s, e in newints:     
                out.write(f'{chrom}\\t{s}\\t{e}\\n')   
                totregions += 1 
        out.close()
        return totregions
    
    def get_intersect_count(filea, fileb):
        import subprocess
        c = 0
        result = subprocess.run(['bedtools', 'intersect', '-f', '0.5', '-u', '-a', filea, '-b', fileb],
                              capture_output=True, text=True)
        for line in result.stdout.rstrip('\\n').split('\\n'):
            if line:
                c += 1
        return c
    
    def extract_sj_info(file):
        read_sjc, read_se_ends = {}, {}
        for line in open(file):
            line = line.rstrip().split('\\t')
            chrom, strand = line[0], line[5]
            iso, start, end = line[3], int(line[1]), int(line[2])
            esizes, estarts = [int(x) for x in line[10].rstrip(',').split(',')], [int(x) for x in line[11].rstrip(',').split(',')]
            exons = [(start+estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
            introns = tuple([(exons[x][1], exons[x+1][0]) for x in range(len(exons)-1)])
            cs = chrom
            if cs not in read_sjc:
                read_sjc[cs] = {}
                read_se_ends[cs] = {}
            if len(introns) > 0:
                if introns not in read_sjc[cs]: read_sjc[cs][introns] = 0
                read_sjc[cs][introns] += 1
            else:
                roundedends = (SINGLE_EXON_END_WINDOW * round(start/SINGLE_EXON_END_WINDOW), 
                                SINGLE_EXON_END_WINDOW * round(end/SINGLE_EXON_END_WINDOW))
                if roundedends not in read_se_ends[cs]: read_se_ends[cs][roundedends] = 0
                read_se_ends[cs][roundedends] += 1
        return read_sjc, read_se_ends
    
    # Open output file
    outfile = open('${dataset_name}_${align_label}_${partition_label}_${process_label}.eval_summary.txt', 'w')
    
    # EVALUATE GENIC REGIONS
    outfile.write("=== GENIC REGION EVALUATION ===\\n")
    outfile.write("file\\ttotal_read_regions\\tfound_regions\\ttotal_reads\\tgenic_reads\\n")
    
    chromtoint, totreads = get_chromtoint(BED_READS_FILE)
    totregions = get_regions(chromtoint, BED_READS_INTERVALS_FILE)
    
    chromtoint_iso, _ = get_chromtoint(ISOFORMS_FILE)
    iso_intervals_file = ISOFORMS_FILE.replace('.bed', '.intervals.bed')
    get_regions(chromtoint_iso, iso_intervals_file)
    
    foundregions = get_intersect_count(BED_READS_INTERVALS_FILE, iso_intervals_file)
    genicreads = get_intersect_count(BED_READS_FILE, iso_intervals_file)
    
    outfile.write(f"{ISOFORMS_FILE}\\t{totregions}\\t{foundregions}\\t{totreads}\\t{genicreads}\\n")
    
    # EVALUATE SPLICE JUNCTIONS
    outfile.write("\\n=== SPLICE JUNCTION EVALUATION ===\\n")
    outfile.write("file\\ttotal_sjc\\tsupported_sjc\\tsubset_sjc\\ttotal_se\\tsupported_se\\n")
    
    read_sjc, read_se_ends = extract_sj_info(BED_READS_FILE)
    found_sjc, found_se_ends = extract_sj_info(ISOFORMS_FILE)
    
    found_subsets = {}
    for cs in found_sjc:
        found_subsets[cs] = set()
        for sjc in found_sjc[cs]:
            for slen in range(len(sjc)-1, 0, -1):
                for i in range(0, len(sjc)-slen+1):
                    found_subsets[cs].add(sjc[i:i+slen])
    
    tot_sjc, sup_sjc, tot_se, sup_se = 0, 0, 0, 0
    subset_sjc = 0
    for cs in read_sjc:
        if cs in found_sjc:
            for sjc in read_sjc[cs]:
                if sjc in found_sjc[cs]:
                    sup_sjc += read_sjc[cs][sjc]
                elif sjc in found_subsets[cs]:
                    subset_sjc += read_sjc[cs][sjc]
                tot_sjc += read_sjc[cs][sjc]
        else:
            for sjc in read_sjc[cs]:
                tot_sjc += read_sjc[cs][sjc]
    
    for cs in read_se_ends:
        if cs in found_se_ends:
            for se in read_se_ends[cs]:
                if se in found_se_ends[cs]:
                    sup_se += read_se_ends[cs][se]
                tot_se += read_se_ends[cs][se]
        else:
            for se in read_se_ends[cs]:
                tot_se += read_se_ends[cs][se]
    
    outfile.write(f"{ISOFORMS_FILE}\\t{tot_sjc}\\t{sup_sjc}\\t{subset_sjc}\\t{tot_se}\\t{sup_se}\\n")
    
    # EVALUATE TRANSCRIPT CLASSIFICATION
    outfile.write("\\n=== TRANSCRIPT CLASSIFICATION (vs ANNOTATION) ===\\n")
    outfile.write("file\\tFSM\\tISM\\tNIC\\tNNC\\tSEM\\tSEN\\n")
    
    transcripttoexons = {}
    for line in open(ANNOT_FILE):
        if line[0] != '#':
            line = line.split('\\t')
            if line[2] == 'exon':
                chrom, strand = line[0], line[6]
                start, end = int(line[3]), int(line[4])
                tname = line[-1].split('transcript_id "')[1].split('"')[0]
                if (chrom, strand) not in transcripttoexons: transcripttoexons[(chrom, strand)] = {}
                if tname not in transcripttoexons[(chrom, strand)]: transcripttoexons[(chrom, strand)][tname] = []
                transcripttoexons[(chrom, strand)][tname].append((start, end))
    
    refjuncs, refjuncchains, refseends = {}, {}, {}
    for chrom, strand in transcripttoexons:
        refjuncs[(chrom, strand)] = set()
        refjuncchains[(chrom, strand)] = set()
        refseends[(chrom, strand)] = set()
        for tname in transcripttoexons[(chrom, strand)]:
            exons = sorted(transcripttoexons[(chrom, strand)][tname])
            if len(exons) > 1:
                introns = [(exons[x][1], exons[x+1][0]-1) for x in range(len(exons)-1)]
                refjuncchains[(chrom, strand)].add(tuple(introns))
                refjuncs[(chrom, strand)].update(set(introns))
            else:
                refseends[(chrom, strand)].add(exons[0])
    
    fsm, ism, nic, nnc, sem, sen, tot = 0, 0, 0, 0, 0, 0, 0
    for line in open(ISOFORMS_FILE):
        line = line.rstrip().split('\\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])
        esizes, estarts = [int(x) for x in line[10].rstrip(',').split(',')], [int(x) for x in line[11].rstrip(',').split(',')]
        exons = [(start+estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
        introns = tuple([(exons[x][1], exons[x+1][0]) for x in range(len(exons)-1)])
        tot += 1
        if len(introns) > 0:
            if introns in refjuncchains[(chrom, strand)]:
                fsm += 1
            else:
                isISM = False
                myjuncstring = str(introns)[1:-1]
                for juncchain in refjuncchains[(chrom, strand)]:
                    if myjuncstring in str(juncchain):
                        isISM = True
                        ism += 1
                        break
                if not isISM:
                    allFound = True
                    for j in introns:
                        if j not in refjuncs[(chrom, strand)]:
                            allFound = False
                            break
                    if allFound: nic += 1
                    else: nnc += 1
        else:
            endsMatch = False
            for refstart, refend in refseends[(chrom, strand)]:
                if abs(start-refstart) < SINGLE_EXON_END_WINDOW and abs(end-refend) < SINGLE_EXON_END_WINDOW:
                    endsMatch = True
                    break
            if endsMatch: sem += 1
            else: sen += 1
    
    outfile.write(f"{ISOFORMS_FILE}\\t{fsm}\\t{ism}\\t{nic}\\t{nnc}\\t{sem}\\t{sen}\\n")
    outfile.close()
    
    print(f"Evaluation complete. Results written to ${dataset_name}_${align_label}_${partition_label}_${process_label}.eval_summary.txt")
    """
}

workflow {
    // =============================================================================
    // DEFINE DATASETS AS OBJECTS
    // =============================================================================
    
    def human_dataset = new Dataset('human', [
        reads: '/private/groups/brookslab/hdheath/projects/flair-eval/cache/datasets/human/WTC11.100reads.fasta',
        genome: '/private/groups/brookslab/reference_genomes/gencode_human/GRCh38.primary_assembly.genome.fa',
        gtf: '/private/groups/brookslab/reference_genomes/gencode_human/gencode.v47.annotation.gtf',
        junctions: '/private/groups/brookslab/hdheath/projects/flair-eval/cache/datasets/human/junctions.bed',
        cage: '/private/groups/brookslab/hdheath/projects/flair-eval/cache/datasets/human/hg38_fair+new_CAGE_peaks_phase1and2.bed',
        quantseq: '/private/groups/brookslab/hdheath/projects/flair-eval/cache/datasets/human/Atlas.clusters.2.0.GRCh38.96.bed'
    ])
    
    def a549_dataset = new Dataset('A549_cDNA', [
        reads: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11.100reads.fasta',
        genome: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/GRCh38.primary_assembly.genome.fa',
        gtf: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/gencode.v48.annotation.gtf',
        junctions: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11_all.SJ.out.tab',
        cage: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/CAGE_TSS_human.bed',
        quantseq: '/private/groups/brookslab/hdheath/projects/test_suite/flair-test-suite/flair-test-suite/tests/data/WTC11_all_polyApeaks_fixed.bed'
    ])
    
    def mouse_dataset = new Dataset('mouse', [
        reads: '/path/to/mouse/reads.fasta',
        genome: '/path/to/mouse/genome.fa',
        gtf: '/path/to/mouse/annotation.gtf',
        junctions: null,
        cage: null,
        quantseq: null
    ])
    
    // =============================================================================
    // DEFINE STANDARD/GLOBAL OPTIONS - default test matrix
    // =============================================================================
    
    def standard_options = [
        align: [
            ['default', ''],
            ['with_nvrna', '--nvrna'],
            ['stringent', '-m 200 --quality 20']
        ],
        partition: [
            ['all', '--all'],
            ['chr1', '--region chr1'],
            ['chr1_1_100k', '--region chr1:1-100000']
        ],
        correct: [
            ['with_gtf', ''],
            ['with_gtf_and_nvrna', '--nvrna'],
            ['stringent', '-w 20']
        ],
        collapse: [
            ['default', ''],
            ['stringent', '-s 4 -e 100'],
            ['check_splice', '--check_splice']
        ]
    ]
    
    // =============================================================================
    // DEFINE TEST SETS
    // =============================================================================
    // Use standard_options when you want the full test matrix on a dataset
    // Use custom options when you want targeted/specific tests
    
    // Test 1: Quick validation - custom minimal options
    def quick_test = new TestSet('quick_test', a549_dataset, [
        align: [['default', '']],
        transcriptome: [['with_gtf', "-f ${a549_dataset.gtf}"]],
        partition: [['chr1_1_100k', '--region chr1:1-100000']],
        correct: [['with_gtf', '']],
        collapse: [['default', '']]
    ])
    // 1 * 1 * 1 * 1 = 1 job
    
    // Test 2: Comprehensive human test - uses STANDARD options for full matrix
    def comprehensive_test = new TestSet('comprehensive_test', human_dataset, standard_options)
    // 3 * 3 * 3 * 3 = 81 jobs - runs all combinations!
    
    // Test 3: Partition comparison - custom focused on partitions
    def partition_test = new TestSet('partition_test', a549_dataset, [
        align: [['default', '']],
        partition: [
            ['all', '--all'],
            ['chr1', '--region chr1'],
            ['chr1_1_100k', '--region chr1:1-100000'],
            ['chr1_100k_200k', '--region chr1:100000-200000']
        ],
        correct: [['with_gtf', '']],
        collapse: [['default', '']]
    ])
    // 1 * 4 * 1 * 1 = 4 jobs
    
    // Test 4: Correct method comparison - custom focused on correct methods
    def correct_comparison_test = new TestSet('correct_comparison', human_dataset, [
        align: [['default', '']],
        partition: [['chr1_1_100k', '--region chr1:1-100000']],
        correct: [
            ['with_gtf', ''],
            ['with_gtf_and_nvrna', '--nvrna'],
            ['stringent', '-w 20']
        ],
        collapse: [['default', '']]
    ])
    // 1 * 1 * 3 * 1 = 3 jobs
    
    // Test 5: Mouse standard test - uses STANDARD options
    def mouse_standard_test = new TestSet('mouse_standard', mouse_dataset, standard_options)
    // 3 * 3 * 3 * 3 = 81 jobs on mouse dataset
    
    // =============================================================================
    // COLLECT ALL TEST SETS TO RUN
    // =============================================================================
    
    def test_sets = [
        quick_test,
        // comprehensive_test,  // Commented out - too many jobs!
        //partition_test,
        //correct_comparison_test,
        // mouse_test  // Commented out - files don't exist yet
    ]
    
    // Print summary
    println("\n=== FLAIR Test Suite ===")
    test_sets.each { test ->
        println("${test.name}: ${test.totalJobs()} jobs using dataset '${test.dataset.name}'")
    }
    println("Total jobs: ${test_sets.sum { it.totalJobs() }}")
    println("========================\n")
    
    // =============================================================================
    // BUILD CHANNELS FROM TEST SETS
    // =============================================================================
    
    // Create align inputs from all test sets
    align_inputs = Channel.from(test_sets)
        .flatMap { test_set ->
            // Get reads as list (handles both single and multiple files)
            def reads_list = test_set.dataset.getReadsList()
            
            // Create one alignment job per read file × align option
            reads_list.collectMany { reads_file ->
                test_set.alignOptions.collect { align_option ->
                    tuple(
                        test_set.name,              // test_name
                        test_set.dataset.name,      // dataset_name
                        file(reads_file),           // reads (single file)
                        align_option[0],            // align_label
                        align_option[1]             // align_args
                    )
                }
            }
        }
    
    // Get unique genomes (in case multiple datasets use same genome)
    unique_genome = Channel.from(test_sets)
        .map { it.dataset.genome }
        .unique()
        .map { file(it) }
        .first()
    
    // Run alignment
    FlairAlign(align_inputs, unique_genome)
    
    // Create partition inputs
    partition_inputs = FlairAlign.out.alignments
        .flatMap { test_name, dataset_name, align_label, bam, bai, bed ->
            // Find the test set that produced this alignment
            def test_set = test_sets.find { it.name == test_name }
            
            test_set.partitionOptions.collect { partition_option ->
                tuple(
                    test_name,
                    dataset_name,
                    align_label,
                    bam,
                    bai,
                    bed,
                    partition_option[0],  // partition_label
                    partition_option[1]   // partition_args
                )
            }
        }
    
    // Prepare optional files using dataset objects
    cage_files = Channel.from(test_sets)
        .map { test_set ->
            def cage_path = test_set.dataset.cage
            tuple(test_set.dataset.name, cage_path ? file(cage_path) : file('NO_CAGE'))
        }
        .unique { it[0] }  // Unique by dataset name
    
    quantseq_files = Channel.from(test_sets)
        .map { test_set ->
            def qs_path = test_set.dataset.quantseq
            tuple(test_set.dataset.name, qs_path ? file(qs_path) : file('NO_QUANTSEQ'))
        }
        .unique { it[0] }
    
    junctions_files = Channel.from(test_sets)
        .map { test_set ->
            def junc_path = test_set.dataset.junctions
            tuple(test_set.dataset.name, junc_path ? file(junc_path) : file('NO_JUNCTIONS'))
        }
        .unique { it[0] }
    
    target_files = Channel.value(file('NO_TARGETS'))
    
    // Get unique GTF for partition process
    unique_gtf = Channel.from(test_sets)
        .map { it.dataset.gtf }
        .unique()
        .map { file(it) }
        .first()
    
    // Run partition
    FlairPartition(
        partition_inputs,
        unique_genome,
        unique_gtf,
        cage_files.map { it[1] }.first(),
        quantseq_files.map { it[1] }.first(),
        junctions_files.map { it[1] }.first(),
        target_files
    )
    
    // Create transcriptome inputs from partitioned outputs
    transcriptome_inputs = FlairPartition.out.partitioned
        .flatMap { test_name, dataset_name, align_label, partition_label, bam, bai, bed, genome, gtf,
                   cage, quantseq, junctions_bed, targets ->
            // Find the test set that produced this partition
            def test_set = test_sets.find { it.name == test_name }
            
            // If no transcriptome options, skip this partition
            if (!test_set.transcriptomeOptions || test_set.transcriptomeOptions.isEmpty()) {
                return []
            }
            
            test_set.transcriptomeOptions.collect { transcriptome_option ->
                tuple(
                    test_name,
                    dataset_name,
                    "${align_label}_${partition_label}",  // combined label
                    bam,
                    bai,
                    bed,
                    transcriptome_option[1]  // transcriptome_args
                )
            }
        }
    
    // Run transcriptome generation on partitioned outputs
    FlairTranscriptome(transcriptome_inputs, unique_genome)
    
    // Create correct inputs
    correct_inputs = FlairPartition.out.partitioned
        .flatMap { test_name, dataset_name, align_label, partition_label, bam, bai, bed, genome, gtf,
                   cage, quantseq, junctions_bed, targets ->
            def test_set = test_sets.find { it.name == test_name }
            
            test_set.correctOptions.collect { correct_option ->
                tuple(
                    test_name,
                    dataset_name,
                    align_label,
                    partition_label,
                    bam,
                    bed,
                    genome,
                    gtf,
                    cage,
                    quantseq,
                    junctions_bed,
                    targets,
                    correct_option[0],  // correct_label
                    correct_option[1]   // correct_args
                )
            }
        }
    
    // Run correct
    FlairCorrect(correct_inputs)
    
    // Create collapse inputs - filter out failed correct tasks
    collapse_inputs = FlairCorrect.out.corrected
        .filter { test_name, dataset_name, align_label, partition_label, correct_label, corrected_bed, inconsistent_bed ->
            // Only proceed if corrected_bed exists and has content
            corrected_bed != null && corrected_bed.size() > 0
        }
        .flatMap { test_name, dataset_name, align_label, partition_label, correct_label, corrected_bed, inconsistent_bed ->
            def test_set = test_sets.find { it.name == test_name }
            
            test_set.collapseOptions.collect { collapse_option ->
                tuple(
                    test_name,
                    dataset_name,
                    align_label,
                    partition_label,
                    correct_label,
                    corrected_bed,
                    inconsistent_bed,
                    file(test_set.dataset.reads),
                    file(test_set.dataset.genome),
                    file(test_set.dataset.gtf),
                    collapse_option[0],  // collapse_label
                    collapse_option[1]   // collapse_args
                )
            }
        }
    
    // Run collapse
    FlairCollapse(collapse_inputs)
    
    // =============================================================================
    // EVALUATION: Run FlairEval on both transcriptome and collapse outputs
    // =============================================================================
    
    // Prepare align outputs for joining - extract test_name, dataset_name, and reads BED
    align_for_eval = FlairAlign.out.alignments
        .map { test_name, dataset_name, align_label, bam, bai, bed ->
            tuple(test_name, dataset_name, bed)
        }
    
    // Prepare transcriptome outputs for evaluation
    eval_transcriptome_inputs = FlairTranscriptome.out.transcriptome
        .map { test_name, dataset_name, align_partition_label, isoforms_bed, isoforms_gtf, isoforms_fa, 
               isoform_counts, isoform_read_map ->
            // Split the combined label back into align_label and partition_label
            def labels = align_partition_label.split('_', 2)  // Split at first underscore
            def align_label = labels[0]
            def partition_label = labels.size() > 1 ? labels[1] : ''
            tuple(test_name, dataset_name, align_label, partition_label, 'transcriptome', isoforms_bed)
        }
        .combine(align_for_eval, by: [0, 1])  // Join by test_name and dataset_name
        .map { test_name, dataset_name, align_label, partition_label, process_label, isoforms_bed, reads_bed ->
            // Get GTF from test set
            def test_set = test_sets.find { it.name == test_name }
            tuple(
                test_name,
                dataset_name,
                align_label,
                partition_label,
                process_label,
                isoforms_bed,
                reads_bed,
                file(test_set.dataset.gtf)
            )
        }
    
    // Prepare collapse outputs for evaluation
    eval_collapse_inputs = FlairCollapse.out.collapsed
        .map { test_name, dataset_name, align_label, partition_label, correct_label, collapse_label,
               isoforms_bed, isoforms_fa, isoforms_gtf, isoform_read_map, isoform_counts ->
            def process_label = "${correct_label}_${collapse_label}"
            tuple(
                test_name,
                dataset_name,
                align_label,
                partition_label,
                process_label,
                isoforms_bed
            )
        }
        .combine(align_for_eval, by: [0, 1])  // Join by test_name and dataset_name
        .map { test_name, dataset_name, align_label, partition_label, process_label, isoforms_bed, reads_bed ->
            // Get GTF from test set
            def test_set = test_sets.find { it.name == test_name }
            tuple(
                test_name,
                dataset_name,
                align_label,
                partition_label,
                process_label,
                isoforms_bed,
                reads_bed,
                file(test_set.dataset.gtf)
            )
        }
    
    // Combine both evaluation inputs and run FlairEval once
    eval_all_inputs = eval_transcriptome_inputs.mix(eval_collapse_inputs)
    FlairEval(eval_all_inputs)
}
