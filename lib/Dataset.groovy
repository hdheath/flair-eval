/**
 * Dataset - Represents a single long-read RNA-seq dataset with associated
 * reference files and orthogonal validation data.
 *
 * Used by the FLAIR evaluation pipeline to track all input files for a sample.
 */
class Dataset {
    String name
    def reads  // Can be String (single file) or List (multiple files)
    String genome
    String gtf
    String cage
    String quantseq
    String library_type  // Library prep / sequencing type (e.g., pacbio_cDNA, ont_cDNA, ont_dRNA)
    String cage_signal_plus   // Optional: CAGE plus-strand bedGraph signal
    String cage_signal_minus  // Optional: CAGE minus-strand bedGraph signal
    String quantseq_signal_plus   // Optional: QuantSeq plus-strand bedGraph signal
    String quantseq_signal_minus  // Optional: QuantSeq minus-strand bedGraph signal
    String bam  // Optional: pre-aligned BAM file
    String bai  // Optional: BAM index file (auto-constructed from bam path)
    String junction_tab  // Optional: short-read junction file for flair transcriptome

    // Constructor
    Dataset(String name, Map config) {
        this.name = name
        this.reads = config.reads
        this.genome = config.genome
        this.gtf = config.gtf
        this.library_type = config.library_type ?: 'unknown'
        this.cage = config.cage
        this.quantseq = config.quantseq
        this.cage_signal_plus = config.cage_signal_plus
        this.cage_signal_minus = config.cage_signal_minus
        this.quantseq_signal_plus = config.quantseq_signal_plus
        this.quantseq_signal_minus = config.quantseq_signal_minus
        this.bam = config.bam
        // Resolve index path from BAM path. Prefer .bai, fall back to .csi if needed.
        if (config.bam) {
            def bamPath = config.bam.toString()
            def baiPath = "${bamPath}.bai"
            def csiPath = "${bamPath}.csi"
            if (new File(baiPath).exists()) {
                this.bai = baiPath
            } else if (new File(csiPath).exists()) {
                this.bai = csiPath
            } else {
                this.bai = baiPath
            }
        } else {
            this.bai = null
        }
        this.junction_tab = config.junction_tab
    }

    // Helper methods
    boolean hasBam() { return bam != null }

    // Get reads as a list (handles both single file and multiple files)
    List<String> getReadsList() {
        return reads instanceof List ? reads : [reads]
    }
}
