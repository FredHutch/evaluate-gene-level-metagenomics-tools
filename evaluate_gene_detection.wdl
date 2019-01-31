import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/77dcf3de35d4e2d8d11810762a37e88880001176/WDL/denovo-assembly-plass/denovo-assembly-plass.wdl" as plass
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/86ee57e3cb8234bee734652b0f93fb10f8fba7ba/tools/make_unique_fasta_headers/make_unique_fasta_headers.wdl?token=AE-VSF9RE3XNpRsIhI9h2sU4ocNuBZ8kks5cXKzVwA" as unique_fasta
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/782ca251f64312bd58bbbeae7e5b115c1508400f/tools/calculate_gene_accuracy/calculate_gene_accuracy.wdl?token=AE-VSDMFO_SAP840GfQFuTBXP4MZKTXJks5cXK5hwA%3D%3D" as calc_acc

workflow evaluateGeneDetection {

  File sample_sheet
  Array[Object] samples = read_objects(sample_sheet)

  scatter (sample in samples) {
    
    call plass.plass {
      input:
        input_fastq=sample.sim_reads
    }
    call unique_fasta.makeUniqueFastaHeaders as plass_unique {
      input:
        fasta_input=plass.contigs
    }
    call calc_acc.calculateGeneAccuracy_wf as plass_calc {
      input:
        ref_fasta=sample.ref_fasta,
        ref_abund=sample.ref_abund,
        detected_fasta=plass.contigs
    }


  }

  output {
    Array[File] plass=plass_calc.aln
  }

}