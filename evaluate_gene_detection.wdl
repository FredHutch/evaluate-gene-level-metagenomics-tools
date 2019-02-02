import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/77dcf3de35d4e2d8d11810762a37e88880001176/WDL/denovo-assembly-plass/denovo-assembly-plass.wdl" as plass
import "https://raw.githubusercontent.com/FredHutch/reproducible-workflows/3c4d7f3f5125b93981315e1c3202dce4952e5fcd/WDL/align-proteins-famli/align-proteins-famli.wdl" as famli
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/1a4791d42b50b4fb87b445a5999cb193651a6693/tools/calculate_gene_accuracy/calculate_gene_accuracy.wdl?token=AE-VSEiy1FYZRGNLtRVAC2o-hcYnRHK7ks5cXeX4wA" as calc_acc
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/1f26244639445b030282c63d60ed4f9d5365a179/tools/extract_detected_genes/extract_detected_genes_famli.wdl?token=AE-VSN7aw0nY2XpbGuT3ZEdMkFbVTqPaks5cXiLrwA" as extract_famli_genes
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/be533c85a18a7807313f84518f80a857bdfd6fd4/tools/make_unique_fasta_headers/make_unique_fasta_headers.wdl?token=AE-VSOKDLxrQh5a8n42Do3DhQNFY5tREks5cXh9LwA" as clean_headers
import "https://raw.githubusercontent.com/FredHutch/evaluate-gene-level-metagenomics-tools/5dfcafb254884d132480a0a470d42190695c6b3c/tools/cluster_proteins_by_identity/cluster_proteins_by_identity.wdl?token=AE-VSPdWLIiGQWDCprmvobg8NMGKhksdks5cXiI8wA" as clust

workflow evaluateGeneDetection {

  File sample_sheet
  String identity="0.9"
  Array[Object] samples = read_objects(sample_sheet)

  scatter (sample in samples) {
    
    # Run Plass
    call plass.plass {
      input:
        input_fastq=sample.sim_reads
    }
    # Make the headers unique
    call clean_headers.makeUniqueFastaHeaders as plass_clean {
      input:
        fasta_input=plass.contigs
    }
    # Cluster the Plass-detected proteins
    call clust.clusterProteins as plass_clust {
      input:
        fasta_in=plass_clean.fasta_output,
        identity=identity
    }

    # Cluster the reference proteins
    call clust.clusterProteinsAbunds as ref_clust {
      input:
        fasta_in=sample.ref_fasta,
        abund_in=sample.ref_abund,
        identity=identity
    }

    # Calculate the accuracy of Plass, with that clustering applied
    call calc_acc.calculateGeneAccuracy as plass_calc {
      input:
        ref_fasta=ref_clust.fasta_out,
        ref_abund=ref_clust.abund_out,
        detected_fasta=plass_clust.fasta_out,
        method_label="Plass"
    }

    # Now run DIAMOND on those clustered proteins
    # Make the database
    call famli.MakeDiamondDatabase {
      input:
        fasta=ref_clust.fasta_out
    }
    # Run the DIAMOND aligner
    call famli.DiamondBlastx {
      input:
        refdb=MakeDiamondDatabase.db,
        input_fastq=sample.sim_reads
    }

    # Filter the results with FAMLI
    call famli.FAMLI {
      input:
        input_aln=DiamondBlastx.aln
    }

    # Extract the genes detected by FAMLI
    call extract_famli_genes.extractDetectedGenesFAMLI as famli_genes {
      input:
        json_input=FAMLI.results,
        fasta_input=ref_clust.fasta_out
    }
    # Calculate the accuracy of the FAMLI results
    call calc_acc.calculateGeneAccuracy as famli_calc {
      input:
        ref_fasta=ref_clust.fasta_out,
        ref_abund=ref_clust.abund_out,
        detected_fasta=famli_genes.fasta_output,
        method_label="FAMLI"
    }

  }

  output {
    Array[File] plass=plass_calc.accuracy
    Array[File] famli=famli_calc.accuracy
  }

}