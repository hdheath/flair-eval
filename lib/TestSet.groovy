/**
 * TestSet - Represents a named collection of modes to evaluate for a given Dataset.
 *
 * Encapsulates the alignment, partition, transcriptome, Bambu, IsoQuant,
 * IsoSeq, FLAMES, and StringTie2 mode configurations that define one
 * "test run" of the evaluation pipeline.
 */
class TestSet {
    String name
    Dataset dataset
    Map alignModes
    Map partitionModes
    Map transcriptomeModes
    Map bambuModes        // Optional: Bambu transcriptome assembler modes
    Map isoquantModes     // Optional: IsoQuant transcriptome assembler modes
    Map isoseqModes       // Optional: PacBio IsoSeq assembler modes
    Map flamesModes       // Optional: FLAMES assembler modes
    Map stringtie2Modes   // Optional: StringTie2 assembler modes

    // Constructor with simplified parameter structure
    TestSet(String name, Dataset dataset, Map modes) {
        this.name = name
        this.dataset = dataset
        this.alignModes = modes.align ?: [:]
        this.partitionModes = modes.partition ?: [:]
        this.transcriptomeModes = modes.transcriptome ?: [:]
        this.bambuModes = modes.bambu ?: [:]
        this.isoquantModes = modes.isoquant ?: [:]
        this.isoseqModes = modes.isoseq ?: [:]
        this.flamesModes = modes.flames ?: [:]
        this.stringtie2Modes = modes.stringtie2 ?: [:]
    }

    // Calculate total number of jobs this test set will produce
    int totalJobs() {
        def align_count = alignModes.size() ?: 1
        def partition_count = partitionModes.size() ?: 1

        return align_count * partition_count
    }
}
