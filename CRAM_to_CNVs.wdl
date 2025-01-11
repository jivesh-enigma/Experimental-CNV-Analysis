# Workflow for running GATK GermlineCNVCaller on multiple case samples using a trained model (obtained from running 
# GATK GermlineCNVCaller in the cohort mode). Supports both WGS and WES.
# taken from: https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/germline/cnv_germline_case_workflow.wdl


# Notes:
#
# - The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
#   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   These intervals will be padded on both sides by the amount specified by padding (default 250)
#   and split into bins of length specified by bin_length (default 1000; specify 0 to skip binning,
#   e.g., for WES).  For WGS, the intervals should simply cover the chromosomes of interest.
#
# - Intervals can be blacklisted from coverage collection and all downstream steps by using the blacklist_intervals
#   argument, which accepts formats compatible with the GATK -XL argument
#   (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   This may be useful for excluding centromeric regions, etc. from analysis.  Alternatively, these regions may
#   be manually filtered from the final callset.
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_germline_case_workflow.wdl -i my_parameters.json
#
#############

version 1.0

import "https://raw.githubusercontent.com/broadinstitute/gatk/refs/heads/master/scripts/cnv_wdl/cnv_common_tasks.wdl" as CNVTasks

workflow CNVGermlineCaseWorkflow {

    input {
        ##################################
        #### required basic arguments ####
        ##################################
        File intervals
        File? blacklist_intervals
        File filtered_intervals
        Array[String]+ normal_bams
        Array[String]+ normal_bais
        File contig_ploidy_model_tar
        Array[File]+ gcnv_model_tars
        Int num_intervals_per_scatter
        File ref_fasta_dict
        File ref_fasta_fai
        File ref_fasta
        String gatk_docker

        String sample_id
        String permanent_bucket_path
        String clinical_bucket_path
        File vep_Rscript_file
        #dbnSFP
        File dbNSFPData
        File dbNSFPDataTbi
        File dbNSFPPlugin
        File dbNSFPReadme

        #   S-CAP
        File scap
        File scap_index

        #   CLINVAR
        File clinvar
        File clinvar_index

        # Gnomad SV
        File gnomad_sv_vcf
        File gnomad_sv_vcf_index

        # Variant filter input

        File master_gene_list
        
        File OMIM_genes
        File Clinvar_genes
        File GCD_genes

        File variant_filter_Rscript_file

        # ClassifyCNV input

        File Python_genotype_to_bed_script

        ##################################
        #### optional basic arguments ####
        ##################################
        File? gatk4_jar_override
        Int? preemptible_attempts

        # Required if BAM/CRAM is in a requester pays bucket
        String? gcs_project_for_requester_pays

        ####################################################
        #### optional arguments for PreprocessIntervals ####
        ####################################################
        Int? padding
        Int? bin_length

        ##############################################
        #### optional arguments for CollectCounts ####
        ##############################################
        Array[String]? disabled_read_filters_for_collect_counts
        String? collect_counts_format
        Boolean? collect_counts_enable_indexing
        Int? mem_gb_for_collect_counts

        ######################################################################
        #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
        ######################################################################
        Float? ploidy_mapping_error_rate
        Float? ploidy_sample_psi_scale
        Int? mem_gb_for_determine_germline_contig_ploidy
        Int? cpu_for_determine_germline_contig_ploidy
        Int? disk_for_determine_germline_contig_ploidy

        ##########################################################
        #### optional arguments for GermlineCNVCallerCaseMode ####
        ##########################################################
        Float? gcnv_p_alt
        Float? gcnv_cnv_coherence_length
        Int? gcnv_max_copy_number
        Int? mem_gb_for_germline_cnv_caller
        Int? cpu_for_germline_cnv_caller
        Int? disk_for_germline_cnv_caller

        # optional arguments for germline CNV denoising model
        Float? gcnv_mapping_error_rate
        Float? gcnv_sample_psi_scale
        Float? gcnv_depth_correction_tau
        String? gcnv_copy_number_posterior_expectation_mode
        Int? gcnv_active_class_padding_hybrid_mode

        # optional arguments for Hybrid ADVI
        Float? gcnv_learning_rate
        Float? gcnv_adamax_beta_1
        Float? gcnv_adamax_beta_2
        Int? gcnv_log_emission_samples_per_round
        Float? gcnv_log_emission_sampling_median_rel_error
        Int? gcnv_log_emission_sampling_rounds
        Int? gcnv_max_advi_iter_first_epoch
        Int? gcnv_max_advi_iter_subsequent_epochs
        Int? gcnv_min_training_epochs
        Int? gcnv_max_training_epochs
        Float? gcnv_initial_temperature
        Int? gcnv_num_thermal_advi_iters
        Int? gcnv_convergence_snr_averaging_window
        Float? gcnv_convergence_snr_trigger_threshold
        Int? gcnv_convergence_snr_countdown_window
        Int? gcnv_max_calling_iters
        Float? gcnv_caller_update_convergence_threshold
        Float? gcnv_caller_internal_admixing_rate
        Float? gcnv_caller_external_admixing_rate
        Boolean? gcnv_disable_annealing

        ###################################################
        #### arguments for PostprocessGermlineCNVCalls ####
        ###################################################
        Int ref_copy_number_autosomal_contigs
        Array[String]? allosomal_contigs
        Int? disk_space_gb_for_postprocess_germline_cnv_calls
        Int? mem_gb_for_postprocess_germline_cnv_calls

        ##########################
        #### arguments for QC ####
        ##########################
        Int maximum_number_events_per_sample
        Int maximum_number_pass_events_per_sample

        # Variant Annotation
        Int diskspace_combineOutputFiles
        Int GATK_diskGb
        Int GATK_memoryGb
        Int HDD_Unzip
        Int RAM_Unzip
        Int bootDiskSizeGb_VEP
        Int cpu_VEP
        Int diskGb_VEP
        Int fork_VEP
        Int memoryGb_VEP
        Int nearestGeneDistance_VEP
        Int HDD_VTRecal
        Int RAM_VTRecal
    }

    Array[Pair[String, String]] normal_bams_and_bais = zip(normal_bams, normal_bais)

    call CNVTasks.PreprocessIntervals {
        input:
            intervals = intervals,
            blacklist_intervals = blacklist_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            padding = padding,
            bin_length = bin_length,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    scatter (normal_bam_and_bai in normal_bams_and_bais) {
        call CNVTasks.CollectCounts {
            input:
                intervals = PreprocessIntervals.preprocessed_intervals,
                bam = normal_bam_and_bai.left,
                bam_idx = normal_bam_and_bai.right,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta_dict = ref_fasta_dict,
                format = collect_counts_format,
                enable_indexing = collect_counts_enable_indexing,
                disabled_read_filters = disabled_read_filters_for_collect_counts,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_collect_counts,
                preemptible_attempts = preemptible_attempts,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }
    }

    call DetermineGermlineContigPloidyCaseMode {
        input:
            read_count_files = CollectCounts.counts,
            contig_ploidy_model_tar = contig_ploidy_model_tar,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_determine_germline_contig_ploidy,
            cpu = cpu_for_determine_germline_contig_ploidy,
            disk_space_gb = disk_for_determine_germline_contig_ploidy,
            mapping_error_rate = ploidy_mapping_error_rate,
            sample_psi_scale = ploidy_sample_psi_scale,
            preemptible_attempts = preemptible_attempts
    }

    call CNVTasks.ScatterIntervals {
        input:
            interval_list = filtered_intervals,
            num_intervals_per_scatter = num_intervals_per_scatter,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    scatter (scatter_index in range(length(ScatterIntervals.scattered_interval_lists))) {
        call GermlineCNVCallerCaseMode {
            input:
                scatter_index = scatter_index,
                read_count_files = CollectCounts.counts,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                gcnv_model_tar = gcnv_model_tars[scatter_index],
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_germline_cnv_caller,
                cpu = cpu_for_germline_cnv_caller,
                p_alt = gcnv_p_alt,
                cnv_coherence_length = gcnv_cnv_coherence_length,
                max_copy_number = gcnv_max_copy_number,
                mapping_error_rate = gcnv_mapping_error_rate,
                sample_psi_scale = gcnv_sample_psi_scale,
                depth_correction_tau = gcnv_depth_correction_tau,
                copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                learning_rate = gcnv_learning_rate,
                adamax_beta_1 = gcnv_adamax_beta_1,
                adamax_beta_2 = gcnv_adamax_beta_2,
                log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                min_training_epochs = gcnv_min_training_epochs,
                max_training_epochs = gcnv_max_training_epochs,
                initial_temperature = gcnv_initial_temperature,
                num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                max_calling_iters = gcnv_max_calling_iters,
                caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                disable_annealing = gcnv_disable_annealing,
                preemptible_attempts = preemptible_attempts
        }
    }

    Array[Array[File]] call_tars_sample_by_shard = transpose(GermlineCNVCallerCaseMode.gcnv_call_tars)

    scatter (sample_index in range(length(normal_bams))) {
        call CNVTasks.PostprocessGermlineCNVCalls {
            input:
                entity_id = CollectCounts.entity_id[sample_index],
                gcnv_calls_tars = call_tars_sample_by_shard[sample_index],
                gcnv_model_tars = gcnv_model_tars,
                calling_configs = GermlineCNVCallerCaseMode.calling_config_json,
                denoising_configs = GermlineCNVCallerCaseMode.denoising_config_json,
                gcnvkernel_version = GermlineCNVCallerCaseMode.gcnvkernel_version_json,
                sharded_interval_lists = GermlineCNVCallerCaseMode.sharded_interval_list,
                allosomal_contigs = allosomal_contigs,
                ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                sample_index = sample_index,
                maximum_number_events = maximum_number_events_per_sample,
                maximum_number_pass_events = maximum_number_pass_events_per_sample,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts
        }
    }

    call CNVTasks.ScatterPloidyCallsBySample {
        input :
            contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
            samples = CollectCounts.entity_id,
            docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    call cnv_data_transfer {
        input:
            sample_id = sample_id,
            bucket_url = permanent_bucket_path,
            genotyped_intervals_vcfs = PostprocessGermlineCNVCalls.genotyped_intervals_vcf[0],
            genotyped_intervals_vcf_indexes = PostprocessGermlineCNVCalls.genotyped_intervals_vcf_index[0],
            genotyped_segments_vcfs = PostprocessGermlineCNVCalls.genotyped_segments_vcf[0],
            genotyped_segments_vcf_indexes = PostprocessGermlineCNVCalls.genotyped_segments_vcf_index[0],
    }

    scatter (segment_vcf in PostprocessGermlineCNVCalls.genotyped_segments_vcf) {

        call UnZip { 
            input:
                vcfFileGz = segment_vcf,
                sampleId = sample_id,
                HDD = HDD_Unzip,
                RAM = RAM_Unzip
        }
        
        call VTRecal {
            input:

                vcfFile=UnZip.vcfFile,
                refFasta=ref_fasta,
                refFastaDict=ref_fasta_dict,
                refFastaIdx=ref_fasta_fai,
                sampleId=sample_id,
                HDD = HDD_VTRecal,
                RAM = RAM_VTRecal
        }
        
        call GATKVariantsToTable {
            input:
                normalizedvcfFileGz=VTRecal.normalizedVCF,
                refFasta=ref_fasta,
                refFastaFai=ref_fasta_fai,
                refFastaDict=ref_fasta_dict,
                samplesetId=sample_id,
                GATK_diskGb = GATK_diskGb,
                GATK_memoryGb = GATK_memoryGb

            
        }
        call vep_task {
            input:
                dbNSFPData = dbNSFPData,
                dbNSFPDataTbi = dbNSFPDataTbi,
                dbNSFPPlugin = dbNSFPPlugin,
                dbNSFPReadme = dbNSFPReadme,
                scap = scap,
                scap_index = scap_index,
                clinvar = clinvar,
                clinvar_index = clinvar_index,
                gnomad_sv_vcf = gnomad_sv_vcf,
                gnomad_sv_vcf_index = gnomad_sv_vcf_index,
                refFasta=ref_fasta,
                refFastaFai=ref_fasta_fai,
                refFastaDict=ref_fasta_dict,
                samplesetId=sample_id,
                normalizedvcfFileGz=VTRecal.normalizedVCF,
                bootDiskSizeGb_VEP = bootDiskSizeGb_VEP,
                cpu_VEP = cpu_VEP,
                diskGb_VEP = diskGb_VEP,
                fork = fork_VEP,
                memoryGb_VEP = memoryGb_VEP,
                nearestGeneDistance = nearestGeneDistance_VEP

        }
        
        call combineOutputFiles {
            input:
                samplesetId=sample_id + ".segment",
                vepOutputFile=vep_task.VEP_Output,
                gatkOutputFile=GATKVariantsToTable.GATK_output,
                Rscript_file = vep_Rscript_file,
                diskSpace = diskspace_combineOutputFiles
        }
        
        File vepannotated_vcfs = combineOutputFiles.vepannotated_vcf
    }

    call genotype_to_bed {
        input:
            sample_id = sample_id,
            vep_genotypes = vepannotated_vcfs[0],
            genotype_to_bed_script = Python_genotype_to_bed_script
    }

    call ClassifyCNV {
        input:
            sample_id = sample_id,
            genotype_bed = genotype_to_bed.genotype_bed
    }

    call MergeCNVResults {
        input:
            sample_id = sample_id,
            vep_genotypes = vepannotated_vcfs[0],
            classifycnv_score = ClassifyCNV.classifycnv_score
    }

    call vep_data_transfer {
        input:
            samplename = sample_id,
            bucket_url = permanent_bucket_path,
            segment_calls_genotypes = vepannotated_vcfs[0],
            cnv_genotypes = MergeCNVResults.CNV_genotypes,
            classifycnv_score = ClassifyCNV.classifycnv_score
    }

    call VariantFilter {
        input:
            input_vcf = MergeCNVResults.CNV_genotypes,
            sample_id = sample_id,
            master_gene_list = master_gene_list,
            OMIM_genes = OMIM_genes,
            Clinvar_genes = Clinvar_genes,
            GCD_genes = GCD_genes,
            Rscript_file = variant_filter_Rscript_file
    }

    call variant_data_transfer_permanent {
        input:
            sample_id = sample_id,
            destination_bucket = permanent_bucket_path,
            gcnv_path_var_HQ = VariantFilter.gcnv_path_var_HQ,
            gcnv_path_var_HQ_non_clinical = VariantFilter.gcnv_path_var_HQ_non_clinical,
            gcnv_path_var_LQ = VariantFilter.gcnv_path_var_LQ,
            gcnv_path_var_LQ_non_clinical = VariantFilter.gcnv_path_var_LQ_non_clinical,
            gcnv_rare_homozygous_variants = VariantFilter.gcnv_rare_homozygous_variants
    }

    call variant_data_transfer_clinical {
        input:
            sample_id = sample_id,
            destination_bucket = clinical_bucket_path,
            gcnv_path_var_HQ = VariantFilter.gcnv_path_var_HQ,
            gcnv_path_var_HQ_non_clinical = VariantFilter.gcnv_path_var_HQ_non_clinical,
            gcnv_path_var_LQ = VariantFilter.gcnv_path_var_LQ,
            gcnv_path_var_LQ_non_clinical = VariantFilter.gcnv_path_var_LQ_non_clinical,
            gcnv_rare_homozygous_variants = VariantFilter.gcnv_rare_homozygous_variants
    }

    output {
        String genotyped_intervals_vcfs = cnv_data_transfer.genotyped_intervals_vcfs_path
        String genotyped_intervals_vcf_indexes = cnv_data_transfer.genotyped_intervals_vcf_indexes_path
        String genotyped_segments_vcfs = cnv_data_transfer.genotyped_segments_vcfs_path
        String genotyped_segments_vcf_indexes = cnv_data_transfer.genotyped_segments_vcf_indexes_path

        String vepannotated_vcf = vep_data_transfer.vepannotated_vcf_path
        String classifycnv_score = vep_data_transfer.classifycnv_score_path
        String cnv_genotypes = vep_data_transfer.cnv_genotypes_path

        String gcnv_path_var_HQ = variant_data_transfer_permanent.gcnv_path_var_HQ_path
        String gcnv_path_var_HQ_non_clinical = variant_data_transfer_permanent.gcnv_path_var_HQ_non_clinical_path
        String gcnv_path_var_LQ = variant_data_transfer_permanent.gcnv_path_var_LQ_path
        String gcnv_path_var_LQ_non_clinical = variant_data_transfer_permanent.gcnv_path_var_LQ_non_clinical_path
        String gcnv_rare_homozygous_variants = variant_data_transfer_permanent.gcnv_rare_homozygous_variants_path
    }
}

task DetermineGermlineContigPloidyCaseMode {
    input {
      Array[File] read_count_files
      File contig_ploidy_model_tar
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts

      # Model parameters
      Float? mapping_error_rate
      Float? sample_psi_scale
    }

    # We do not expose Hybrid ADVI parameters -- the default values are decent

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{default=8 cpu}
        export OMP_NUM_THREADS=~{default=8 cpu}

        mkdir contig-ploidy-model
        tar xzf ~{contig_ploidy_model_tar} -C contig-ploidy-model

        gatk --java-options "-Xmx~{command_mem_mb}m" DetermineGermlineContigPloidy \
            --input ~{sep=" --input " read_count_files} \
            --model contig-ploidy-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --mapping-error-rate ~{default="0.3" mapping_error_rate} \
            --sample-psi-scale ~{default="0.0001" sample_psi_scale}

        tar c -C ~{output_dir_}/case-calls . | gzip -1 > case-contig-ploidy-calls.tar.gz

        rm -rf contig-ploidy-model
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File contig_ploidy_calls_tar = "case-contig-ploidy-calls.tar.gz"
    }
}

task GermlineCNVCallerCaseMode {
    input {
      Int scatter_index
      Array[File] read_count_files
      File contig_ploidy_calls_tar
      File gcnv_model_tar
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts

      # Caller parameters
      Float? p_alt
      Float? cnv_coherence_length
      Int? max_copy_number

      # Denoising model parameters
      Float? mapping_error_rate
      Float? sample_psi_scale
      Float? depth_correction_tau
      String? copy_number_posterior_expectation_mode
      Int? active_class_padding_hybrid_mode

      # Hybrid ADVI parameters
      Float? learning_rate
      Float? adamax_beta_1
      Float? adamax_beta_2
      Int? log_emission_samples_per_round
      Float? log_emission_sampling_median_rel_error
      Int? log_emission_sampling_rounds
      Int? max_advi_iter_first_epoch
      Int? max_advi_iter_subsequent_epochs
      Int? min_training_epochs
      Int? max_training_epochs
      Float? initial_temperature
      Int? num_thermal_advi_iters
      Int? convergence_snr_averaging_window
      Float? convergence_snr_trigger_threshold
      Int? convergence_snr_countdown_window
      Int? max_calling_iters
      Float? caller_update_convergence_threshold
      Float? caller_internal_admixing_rate
      Float? caller_external_admixing_rate
      Boolean? disable_annealing
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])
    Int num_samples = length(read_count_files)

    command <<<
        set -eu
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{default=8 cpu}
        export OMP_NUM_THREADS=~{default=8 cpu}

        mkdir contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls

        mkdir gcnv-model
        tar xzf ~{gcnv_model_tar} -C gcnv-model

        gatk --java-options "-Xmx~{command_mem_mb}m"  GermlineCNVCaller \
            --run-mode CASE \
            --input ~{sep=" --input " read_count_files} \
            --contig-ploidy-calls contig-ploidy-calls \
            --model gcnv-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --p-alt ~{default="5e-4" p_alt} \
            --cnv-coherence-length ~{default="10000.0" cnv_coherence_length} \
            --max-copy-number ~{default="5" max_copy_number} \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --sample-psi-scale ~{default="0.01" sample_psi_scale} \
            --depth-correction-tau ~{default="10000.0" depth_correction_tau} \
            --copy-number-posterior-expectation-mode ~{default="HYBRID" copy_number_posterior_expectation_mode} \
            --active-class-padding-hybrid-mode ~{default="50000" active_class_padding_hybrid_mode} \
            --learning-rate ~{default="0.05" learning_rate} \
            --adamax-beta-1 ~{default="0.9" adamax_beta_1} \
            --adamax-beta-2 ~{default="0.99" adamax_beta_2} \
            --log-emission-samples-per-round ~{default="50" log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error ~{default="0.005" log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds ~{default="10" log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch ~{default="5000" max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs ~{default="100" max_advi_iter_subsequent_epochs} \
            --min-training-epochs ~{default="10" min_training_epochs} \
            --max-training-epochs ~{default="100" max_training_epochs} \
            --initial-temperature ~{default="2.0" initial_temperature} \
            --num-thermal-advi-iters ~{default="2500" num_thermal_advi_iters} \
            --convergence-snr-averaging-window ~{default="500" convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold ~{default="0.1" convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window ~{default="10" convergence_snr_countdown_window} \
            --max-calling-iters ~{default="10" max_calling_iters} \
            --caller-update-convergence-threshold ~{default="0.001" caller_update_convergence_threshold} \
            --caller-internal-admixing-rate ~{default="0.75" caller_internal_admixing_rate} \
            --caller-external-admixing-rate ~{default="1.00" caller_external_admixing_rate} \
            --disable-annealing ~{default="false" disable_annealing}

        tar czf case-gcnv-tracking-shard-~{scatter_index}.tar.gz -C ~{output_dir_}/case-tracking .

        CURRENT_SAMPLE=0
        NUM_SAMPLES=~{num_samples}
        NUM_DIGITS=${#NUM_SAMPLES}
        while [ $CURRENT_SAMPLE -lt $NUM_SAMPLES ]; do
            CURRENT_SAMPLE_WITH_LEADING_ZEROS=$(printf "%0${NUM_DIGITS}d" $CURRENT_SAMPLE)
            tar czf case-gcnv-calls-shard-~{scatter_index}-sample-$CURRENT_SAMPLE_WITH_LEADING_ZEROS.tar.gz -C ~{output_dir_}/case-calls/SAMPLE_$CURRENT_SAMPLE .
            let CURRENT_SAMPLE=CURRENT_SAMPLE+1
        done

        rm -rf contig-ploidy-calls
        rm -rf gcnv-model
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[File] gcnv_call_tars = glob("case-gcnv-calls-shard-~{scatter_index}-sample-*.tar.gz")
        File gcnv_tracking_tar = "case-gcnv-tracking-shard-~{scatter_index}.tar.gz"
        File calling_config_json = "~{output_dir_}/case-calls/calling_config.json"
        File denoising_config_json = "~{output_dir_}/case-calls/denoising_config.json"
        File gcnvkernel_version_json = "~{output_dir_}/case-calls/gcnvkernel_version.json"
        File sharded_interval_list = "~{output_dir_}/case-calls/interval_list.tsv"
    }
}

task cnv_data_transfer {
    input {
        String sample_id
        String bucket_url
        File genotyped_intervals_vcfs
        File genotyped_intervals_vcf_indexes
        File genotyped_segments_vcfs
        File genotyped_segments_vcf_indexes
        Int additional_disk_gb = 10
    }

    Int input_size = ceil(size(genotyped_intervals_vcfs,"GB") + size(genotyped_intervals_vcf_indexes,"GB") + size(genotyped_segments_vcfs,"GB") + size(genotyped_segments_vcf_indexes,"GB"))
    Int disk_gb = input_size + additional_disk_gb

    command <<<

        # echo -e "\n=== Sending Genotyped Intervals VCF ===\n"

        # for file in ~{sep=" " genotyped_intervals_vcfs}; do
        # gsutil -m cp ${file} ~{bucket_url}/~{sample_id}/gcnv_genotyped_intervals_vcfs/
        # echo "~{bucket_url}/~{sample_id}/gcnv_genotyped_intervals_vcfs/$(basename ${file})" >> genotyped_intervals_vcfs_gs_urls.txt
        # done

        # echo -e "\n=== Sending Genotyped Intervals VCF Indexes ===\n"

        # for file in ~{sep=" " genotyped_intervals_vcf_indexes}; do
        # gsutil -m cp ${file} ~{bucket_url}/~{sample_id}/gcnv_genotyped_intervals_vcf_indexes/
        # echo "~{bucket_url}/~{sample_id}/gcnv_genotyped_intervals_vcf_indexes/$(basename ${file})" >> genotyped_intervals_vcf_indexes_gs_urls.txt
        # done

        # echo -e "\n=== Sending Genotyped Segments VCF ===\n"

        # for file in ~{sep=" " genotyped_segments_vcfs}; do
        # gsutil -m cp ${file} ~{bucket_url}/~{sample_id}/gcnv_genotyped_segments_vcfs/
        # echo "~{bucket_url}/~{sample_id}/gcnv_genotyped_segments_vcfs/$(basename ${file})" >> genotyped_segments_vcfs_gs_urls.txt
        # done

        # echo -e "\n=== Sending Genotyped Segments VCF Indexes ===\n"

        # for file in ~{sep=" " genotyped_segments_vcf_indexes}; do
        # gsutil -m cp ${file} ~{bucket_url}/~{sample_id}/gcnv_genotyped_segments_vcf_indexes/
        # echo "~{bucket_url}/~{sample_id}/gcnv_genotyped_segments_vcf_indexes/$(basename ${file})" >> genotyped_segments_vcf_indexes_gs_urls.txt
        # done

        gsutil -m cp ~{genotyped_intervals_vcfs} ~{bucket_url}/~{sample_id}/~{sample_id}.GCNV_intervals_genotyped.vcf.gz
        gsutil -m cp ~{genotyped_intervals_vcf_indexes} ~{bucket_url}/~{sample_id}/~{sample_id}.GCNV_intervals_genotyped.vcf.gz.tbi
        gsutil -m cp ~{genotyped_segments_vcfs} ~{bucket_url}/~{sample_id}/~{sample_id}.GCNV_segments_genotyped.vcf.gz
        gsutil -m cp ~{genotyped_segments_vcf_indexes} ~{bucket_url}/~{sample_id}/~{sample_id}.GCNV_segments_genotyped.vcf.gz.tbi

    >>>

    output {
        String genotyped_intervals_vcfs_path = "~{bucket_url}/~{sample_id}/~{sample_id}.GCNV_intervals_genotyped.vcf.gz"
        String genotyped_intervals_vcf_indexes_path = "~{bucket_url}/~{sample_id}/~{sample_id}.GCNV_intervals_genotyped.vcf.gz.tbi"
        String genotyped_segments_vcfs_path = "~{bucket_url}/~{sample_id}/~{sample_id}.GCNV_segments_genotyped.vcf.gz"
        String genotyped_segments_vcf_indexes_path = "~{bucket_url}/~{sample_id}/~{sample_id}.GCNV_segments_genotyped.vcf.gz.tbi"
    }

    runtime {
        docker: "google/cloud-sdk"
        memory: "8GB"
        disks: 'local-disk ' + disk_gb + ' HDD' 
        preemptible: 5
    }

    meta {
        description: "This task transfers files to a Google Storage bucket and returns the gs URLs."
    }

    parameter_meta {
        sample_id: "The ID of the sample"
        bucket_url: "The Google Storage bucket URL"
    }
}

task UnZip {
    input {
        File vcfFileGz
        String sampleId
        Int RAM
        Int HDD
    }
    

    command {
    # Decompress bgzipped merged VCF file
    echo "bgzip -d -c ~{vcfFileGz} > vcfFile.vcf"
    bgzip -d -c ~{vcfFileGz} > ~{sampleId}.vcf
    }

    output {
        File vcfFile="~{sampleId}.vcf"
    }
    runtime {
        docker: "vanallenlab/vt:3.13.2018"
        memory: RAM + "GB"
        disks: "local-disk " + HDD + " HDD"
        preemptible: "3"
    }
}


task VTRecal {
    input {
        File vcfFile 
        File refFasta
        File refFastaIdx
        File refFastaDict
        String sampleId
        Int RAM
        Int HDD
    }

    command {

        echo "########### decompose VCF"
        /software/vt/./vt decompose -s \
        -o ~{sampleId}.vt1.vcf \
        ~{vcfFile}

       # echo "########### normalize VCF using ch38 genome build"
        #/software/vt/./vt normalize \
        #-r ~{refFasta} \
        #-o ~{sampleId}.vt2.vcf \
        #~{sampleId}.vt1.vcf
        
        bgzip ~{sampleId}.vt1.vcf
        
       # echo "########### normalizing the spanning alleles (*):"
       # sed 's/*/-/g' ~{sampleId}.vt2.vcf > ~{sampleId}.vt2_normalized_spanning_alleles.vcf
       # bgzip ~{sampleId}.vt2_normalized_spanning_alleles.vcf
        
        echo "########### creating an index for vcf.gz:"
        tabix -p vcf ~{sampleId}.vt1.vcf.gz 
    }

    output {
        File normalizedVCF="~{sampleId}.vt1.vcf.gz"
        File normalizedVCF_index="~{sampleId}.vt1.vcf.gz.tbi"
    }

    runtime {
        docker: "vanallenlab/vt:3.13.2018"
        memory: RAM + " GB"
        disks: "local-disk " + HDD + " HDD"
        preemptible: "3"
    }
}


task vep_task {
    input {
        File normalizedvcfFileGz
        String samplesetId
        # Customizations
        Int nearestGeneDistance
        # Optimizations
        Int fork

        File refFasta
        File refFastaDict
        File refFastaFai
        Int bootDiskSizeGb_VEP
        Int cpu_VEP
        Int diskGb_VEP
        Int memoryGb_VEP
        
        # Cache files
        File speciesCacheTarGzFile="gs://fc-268e3277-9518-4b28-9f5e-1664c9c5c093/ch38/annotation_files/VEP_files/homo_sapiens_merged_vep_109_GRCh38.tar.gz"
        String speciesCacheLabel="homo_sapiens_merged"
        String speciesCacheParameter="--merged"

        #dbnSFP
        File dbNSFPData
        File dbNSFPDataTbi
        File dbNSFPPlugin
        File dbNSFPReadme

        #   S-CAP
        File scap
        File scap_index

        #   CLINVAR
        File clinvar
        File clinvar_index

        # Gnomad SV
        File gnomad_sv_vcf
        File gnomad_sv_vcf_index
    }

    command <<<
        # Prepare directories for data
        echo "mkdir /opt/vep/.vep/~{speciesCacheLabel}/"
        mkdir /opt
        mkdir /opt/vep
        mkdir /opt/vep/.vep
        mkdir /opt/vep/.vep/~{speciesCacheLabel}/
        mkdir /opt/vep/.vep/Plugins
        mkdir /opt/vep/.vep/Plugins/data
        
        # Put Plugins in correct folder
        echo "symbolic links..."
       
      
        # dbNSFP
        ln -s ~{dbNSFPPlugin} /opt/vep/.vep/Plugins/dbNSFP.pm
        ln -s ~{dbNSFPData} /opt/vep/.vep/Plugins/data/dbNSFPv4.2a_modified_header.gz
        ln -s ~{dbNSFPDataTbi} /opt/vep/.vep/Plugins/data/dbNSFPv4.2a_modified_header.gz.tbi
        ln -s ~{dbNSFPReadme} /opt/vep/.vep/Plugins/data/dbNSFP4.1a.readme.txt
      
        
        # Uncompress the species cache to the .vep directory
        echo "# Uncompress the species cache to the .vep directory"
        echo "tar xzf ~{speciesCacheTarGzFile} -C ~"
        tar xzf ~{speciesCacheTarGzFile} -C .
        
        echo "ls -lh"
        ls -lh
        
        echo "mv ~{speciesCacheLabel}/* /opt/vep/.vep/~{speciesCacheLabel}/"
        mv ~{speciesCacheLabel}/* /opt/vep/.vep/~{speciesCacheLabel}
        echo "ls -lh /opt/vep/.vep/~{speciesCacheLabel}/109_GRCh38/*"
        ls -lh /opt/vep/.vep/~{speciesCacheLabel}/109_GRCh38/*
    
        # log progress to make sure that the VEP output is being generated
        set -xeuo pipefail
        function runtimeInfo() {            
            echo [$(date)]
            echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
            echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
            echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
            do runtimeInfo;
            sleep 30;
        done &
   
        
        
        # Run VEP
        echo "running VEP..."
        vep -v -i ~{normalizedvcfFileGz} -o ~{samplesetId}.vep.txt \
        --tab \
        --offline --cache ~{speciesCacheParameter} --dir /opt/vep/.vep --fasta ~{refFasta} \
        --force_overwrite --stats_text --symbol --everything \
        --regulatory --distance ~{nearestGeneDistance}  \
        --total_length --numbers --domains --pick --pick_order rank,mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,length,ensembl,refseq --variant_class --hgvs --hgvsg --ccds  --fork ~{fork} \
        --custom ~{scap},scap_v1.0,vcf,exact,0,Allele_region,rawscore,sensscore,rawscore_dom,sensscore_dom,rawscore_rec,senscore_rec \
        --plugin dbNSFP,/opt/vep/.vep/Plugins/data/dbNSFPv4.2a_modified_header.gz,hg19_chr,hg19_pos,FATHMM_score,FATHMM_pred,PROVEAN_score,MetaSVM_score,MetaLR_score,MetaLR_pred,MetaRNN_score,MetaRNN_pred,M-CAP_score,M-CAP_pred,REVEL_score,MutPred_score,MVP_score,Aloft_pred,LINSIGHT,CADD_raw,GenoCanyon_score,integrated_fitCons_score,Interpro_domain,gnomAD_genomes_MID_AC,gnomAD_genomes_MID_AN,gnomAD_genomes_MID_AF,gnomAD_genomes_MID_nhomalt \
        --custom ~{clinvar},ClinVar_updated_2023Feb,vcf,exact,0,ID,ALLELEID,CLNDN,CLNDISDB,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNVI,DBVARID \
        --custom ~{gnomad_sv_vcf},GnomadSV_V4,vcf,exact,0,ID,FILTER,ALGORITHMS,CHR2,AC,AN,AF,EVIDENCE 
    
        
        
        echo "ls -lh"
        ls -lh
        
        echo "ls -lh opt/vep/.vep/*"
        ls -lh /opt/vep/.vep/*
        
        echo "Number of VEP variants (grep -v # ~{samplesetId}.vep.txt | wc -l):"
        grep -v "#" ~{samplesetId}.vep.txt | wc -l
        
        gzip ~{samplesetId}.vep.txt
    >>>
    
    runtime {
        docker: "ensemblorg/ensembl-vep:release_109.2"    
        bootDiskSizeGb : "~{bootDiskSizeGb_VEP}"
        preemptible    : 0
        cpu            : "~{cpu_VEP}"
        disks          : "local-disk ~{diskGb_VEP} SSD"
        memory         : "~{memoryGb_VEP} GB"
    }

    output {        
        File VEP_Output="~{samplesetId}.vep.txt.gz"
        File VEP_Summary="~{samplesetId}.vep.txt_summary.txt"
    }
}


task GATKVariantsToTable {
    input {
        File normalizedvcfFileGz
        File refFasta
        File refFastaDict
        File refFastaFai
        String samplesetId
        Int GATK_diskGb
        Int GATK_memoryGb
    }
    
    
    command <<<      
      echo "ls -lh"
      ls -lh
      ls ~{normalizedvcfFileGz}
      
      mv ~{normalizedvcfFileGz} vcfFile.vcf.gz
      
      echo "ls -lh"
      ls -lh
      
      echo "bgzip decompressing vt recal VCF file"
      bgzip --decompress vcfFile.vcf.gz
    
      echo "ls -lh"
      ls -lh 
      
      echo "ls -lh vcfFile.vcf"
      ls -lh vcfFile.vcf
      
      echo "########### Using GATK to extract variants into a table format (GRCh38)"
      java -jar /usr/GenomeAnalysisTK.jar -R ~{refFasta} -T VariantsToTable \
      -V vcfFile.vcf -o ~{samplesetId}.vt2_GATK_annotations.txt \
      -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F END -GF GT -GF CN -GF NP -GF QA -GF QS -GF QSE -GF QSS --allowMissingData --showFiltered
      #-F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F ALGORITHMS -F CHR2  -F CN_COUNT -F CN_FREQ -F CN_NONREF_COUNT -F CN_NONREF_FREQ -F CN_NUMBER -F END -F EVIDENCE -F PREDICTED_INTERGENIC -F PREDICTED_NEAREST_TSS -F PREDICTED_NONCODING_SPAN -F SVLEN -F SVTYPE -F PROTEIN_CODING__NEAREST_TSS -F PROTEIN_CODING__INTERGENIC -F NONCODING_SPAN   -GF GT -GF GQ -GF RD_CN -GF RD_GQ -GF PE_GT -GF PE_GQ -GF SR_GT -GF SR_GQ -GF EV -GF CN -GF CNQ --allowMissingData --showFiltered 

      echo '########### Done extracting relevant fields using GATK'

      # count the number of GATK variants:
      echo "########### number of GATK variants: "
      cat ~{samplesetId}.vt2_GATK_annotations.txt | wc -l
      
      gzip ~{samplesetId}.vt2_GATK_annotations.txt
      
      echo "ls -lh"
      ls -lh
    >>>
    
    output {
        File GATK_output="~{samplesetId}.vt2_GATK_annotations.txt.gz"
    }
    
    runtime {
        docker: "vanallenlab/gatk3.7_with_htslib1.9:1.0"
        memory: "~{GATK_memoryGb} GB"
        cpu: "1"
        disks: "local-disk ~{GATK_diskGb} HDD"
        preemptible: 0
    }
}


task combineOutputFiles {
    input {
        File vepOutputFile
        File gatkOutputFile
        File Rscript_file
        String samplesetId
        Int diskSpace
    }

    command <<<
      cp ~{vepOutputFile} vepOutputFile.txt.gz
      cp ~{gatkOutputFile} gatkOutputFile.txt.gz
    
      gunzip vepOutputFile.txt.gz
      gunzip gatkOutputFile.txt.gz
    
      # remove the '#' from the first line before parsing:
      echo "########### removing the # sign from the first line of the VEP output file"
      sed  "s/#Uploaded_variation/Uploaded_variation/g" vepOutputFile.txt > ~{samplesetId}_vt2_VEP_temp2.txt

      # remove the excess header part:
      grep "#" ~{samplesetId}_vt2_VEP_temp2.txt > ~{samplesetId}_VEP_annotation_list.txt
      grep -v "##" ~{samplesetId}_vt2_VEP_temp2.txt > ~{samplesetId}_vt2_VEP.txt
      #tail -n +124 ~{samplesetId}_vt2_VEP_temp2.txt > ~{samplesetId}_vt2_VEP.txt

      # count the number of VEP variants:
      echo "########### number of VEP variants"
      cat ~{samplesetId}_vt2_VEP.txt | wc -l
      
      echo "########### Combining VEP  output files"
      #paste ~{samplesetId}_vt2_VEP.txt gatkOutputFile.txt | bgzip > ~{samplesetId}_gCNV_VEP_Genotypes.txt.gz
      Rscript ~{Rscript_file} ~{samplesetId}_vt2_VEP.txt ~{gatkOutputFile} ~{samplesetId}
      
    >>>
    
    output {
        File vepannotated_vcf="~{samplesetId}_gcnv_calls_genotypes.txt.gz"
        File filtered_calls="~{samplesetId}_gcnv_filtered_out_calls.txt"

    }
    
    runtime {
        docker: "lbwang/rocker-genome"    
        preemptible    : 0
        cpu            : "1"
        disks          : "local-disk ~{diskSpace} HDD"
        memory         : "14 GB"
    }
}

task vep_data_transfer {
    input {
        String samplename
        String bucket_url
        #Array[File] segment_calls_genotypes
        File segment_calls_genotypes
        File cnv_genotypes
        File classifycnv_score
        Int additional_disk_gb = 10
    }
    
    Int input_size = ceil(size(segment_calls_genotypes,"GB"))
    Int disk_gb = input_size + additional_disk_gb

    command <<<

        # echo -e "\n=== Sending VEP annotated VCFs ===\n"

        # for file in ~{sep=" " segment_calls_genotypes}; do
        # gsutil -m cp ${file} ~{bucket_url}/~{samplename}/segment_calls_genotypes/
        # echo "~{bucket_url}/~{samplename}/gcnv_segment_calls_genotypes/$(basename ${file})" >> segment_calls_genotypes_gs_urls.txt
        # done

        gsutil -m cp ~{segment_calls_genotypes} ~{bucket_url}/~{samplename}/~{samplename}.segment_gcnv_calls_genotypes.txt.gz

        gsutil -m cp ~{classifycnv_score} ~{bucket_url}/~{samplename}/~{samplename}.ClassifyCNV_Scoresheet.txt

        gsutil -m cp ~{cnv_genotypes} ~{bucket_url}/~{samplename}/~{samplename}.CNV_merged_segment_gcnv_calls_genotypes.txt.gz

    >>>

    output {
        String vepannotated_vcf_path = "~{bucket_url}/~{samplename}/~{samplename}.segment_gcnv_calls_genotypes.txt.gz"
        String classifycnv_score_path = "~{bucket_url}/~{samplename}/~{samplename}.ClassifyCNV_Scoresheet.txt"
        String cnv_genotypes_path = "~{bucket_url}/~{samplename}/~{samplename}.CNV_merged_segment_gcnv_calls_genotypes.txt.gz"
    }

    runtime {
        docker: "google/cloud-sdk"
        memory: "8GB"
        disks: 'local-disk ' + disk_gb + ' HDD' 
        preemptible: 5
    }

    meta {
        description: "This task transfers files to a Google Storage bucket and returns the gs URLs."
    }

    parameter_meta {
        samplename: "The ID of the sample"
        bucket_url: "The Google Storage bucket URL"
    }
}

task VariantFilter {
    input {
        File input_vcf
        File master_gene_list
        
        File OMIM_genes
        File Clinvar_genes
        File GCD_genes

        File Rscript_file
        
        Int ram_gb = 32
        Int boot_disk_gb = 24
        Int output_disk_gb = 20
        Int preemptible = 2
        String sample_id
}
   command <<<
       Rscript ~{Rscript_file} ~{master_gene_list} ~{OMIM_genes} ~{Clinvar_genes} ~{GCD_genes} ~{input_vcf} ~{sample_id}
   >>>
   
   runtime { 
     docker : "lbwang/rocker-genome"
     bootDiskSizeGb: "~{boot_disk_gb}"
     preemptible : preemptible
     disks: "local-disk ~{output_disk_gb} HDD"
     memory: "~{ram_gb}GB"
    }   
    output {
        File gcnv_path_var_HQ="~{sample_id}.pathogenic_CNV_calls_high_quality.csv"
        File gcnv_path_var_HQ_non_clinical="~{sample_id}.pathogenic_CNV_calls_high_quality_non_clinical.csv"        
        File gcnv_path_var_LQ="~{sample_id}.pathogenic_CNV_calls_low_quality.csv"
        File gcnv_path_var_LQ_non_clinical = "~{sample_id}.pathogenic_CNV_calls_low_quality_non_clinical.csv"
        File gcnv_rare_homozygous_variants = "~{sample_id}.rare_homozygous_CNV_calls.csv"
        
       }  
}

task genotype_to_bed {
    input {
        String sample_id
        File vep_genotypes
        File genotype_to_bed_script

        Int disk_gb = ceil(size(vep_genotypes,"GiB")) + 10
        Int mem_gb = 4
        Int preemptible = 2
        String genotype_to_bed_docker = "jiveshenigma/terra-pgs-env:v3"
    }

    command <<<
        python3 ~{genotype_to_bed_script} ~{vep_genotypes} ~{sample_id}.segment_gcnv_genotypes.bed
    >>>

    output {
        File genotype_bed = "~{sample_id}.segment_gcnv_genotypes.bed"
    }

    runtime {
        docker: genotype_to_bed_docker
        disks: "local-disk ~{disk_gb} HDD"
        memory: "~{mem_gb}GB"
        preemptible: preemptible
    }
}

task ClassifyCNV {
    input {
        String sample_id
        File genotype_bed
        
        String classifycnv_docker = "jiveshenigma/classifycnv:v1"
        Int disk_gb = ceil(size(genotype_bed,"GiB")) + 10
        Int mem_gb = 4
        Int preemptible = 2
    }

    command <<<
        bedtools --version
        python3 /app/ClassifyCNV/ClassifyCNV.py --help

        python3 /app/ClassifyCNV/ClassifyCNV.py --infile ~{genotype_bed} --GenomeBuild hg38 --precise --outdir ~{sample_id}.classifycnv_result

        # Renaming scorefile
        mv ~{sample_id}.classifycnv_result/Scoresheet.txt ~{sample_id}.classifycnv_result/~{sample_id}.ClassifyCNV_Scoresheet.txt
    >>>

    output {
        File classifycnv_score = "~{sample_id}.classifycnv_result/~{sample_id}.ClassifyCNV_Scoresheet.txt"
    }

    runtime {
        docker: classifycnv_docker
        disks: "local-disk ~{disk_gb} HDD"
        memory: "~{mem_gb}GB"
        preemptible: preemptible
    }
}

task MergeCNVResults {
    input {
        String sample_id
        File vep_genotypes
        File classifycnv_score

        Int disk_gb = ceil(size(vep_genotypes,"GiB")) + 10
        Int mem_gb = 4
        Int preemptible = 2
        String MergeCNVResults_docker = "jiveshenigma/terra-pgs-env:v3"
    }

    command <<<
        python3 <<CODE

        import pandas as pd

        # Load the TSV files into pandas DataFrames
        cnvscore = pd.read_csv("~{classifycnv_score}", sep="\t")
        genotypes = pd.read_csv("~{vep_genotypes}", sep="\t")

        # Function to extract the core identifier
        def get_core_id(id_str, id):
            """Extract core parts of the ID by removing known prefixes or suffixes."""
            # CNV_chr10_46390621_46391180
            # chr10_46390621_46391180_DEL
            if pd.isna(id_str):  # Handle NaN values
                return None
            # Split by common separators and join the core parts
            if id == 'Uploaded_variation':
                return "_".join(id_str.split("_")[1:])  # Skip the first part (like CNV)
            elif id == 'VariantID':
                return "_".join(id_str.split("_")[:-1])  # Skip the last part (like DEL/DUP/INS)
            else:
                print("==ERROR==\nUnknown Variant Identifier column detected. Please check the column names. Accepted values ('Uploaded_variation' from vep genotypes, 'VariantID' from classifycnv)\n")
                return None

        # Apply the function to extract core IDs
        genotypes['Core_ID'] = genotypes['Uploaded_variation'].apply(get_core_id, args=('Uploaded_variation',))
        cnvscore['Core_ID'] = cnvscore['VariantID'].apply(get_core_id, args=('VariantID',))

        # Specify the columns to append from cnvscore.txt
        columns_to_append = [
            'VariantID', 'Chromosome', 'Start', 'End', 'Type', 
            'Classification', 'Total score', 
            'Known or predicted dosage-sensitive genes', 
            'All protein coding genes'
        ]

        # Perform the join based on the extracted core IDs
        merged_df = genotypes.merge(
            cnvscore[columns_to_append + ['Core_ID']], 
            on='Core_ID', 
            how='left'
        )

        # Drop the temporary "Core_ID" column
        merged_df.drop(columns=['Core_ID'], inplace=True)

        # Save the merged data to a new file
        merged_df.to_csv("~{sample_id}.CNV_merged_segment_gcnv_calls_genotypes.txt.gz", sep="\t", index=False)

        CODE
    >>>

    output {
        File CNV_genotypes = "~{sample_id}.CNV_merged_segment_gcnv_calls_genotypes.txt.gz"
    }

    runtime {
        docker: MergeCNVResults_docker
        disks: "local-disk ~{disk_gb} HDD"
        memory: "~{mem_gb}GB"
        preemptible: preemptible
    }
}

task variant_data_transfer_permanent {
    input {
        String sample_id
        String destination_bucket
        File gcnv_path_var_HQ
        File gcnv_path_var_HQ_non_clinical
        File gcnv_path_var_LQ
        File gcnv_path_var_LQ_non_clinical
        File gcnv_rare_homozygous_variants

        Int disk_gb = 16
        Int memory_gb = 8
        Int preemptible = 2

    }

    command <<<
        gsutil -m cp ~{gcnv_path_var_HQ} ~{destination_bucket}/~{sample_id}/
        gsutil -m cp ~{gcnv_path_var_HQ_non_clinical} ~{destination_bucket}/~{sample_id}/
        gsutil -m cp ~{gcnv_path_var_LQ} ~{destination_bucket}/~{sample_id}/
        gsutil -m cp ~{gcnv_path_var_LQ_non_clinical} ~{destination_bucket}/~{sample_id}/
        gsutil -m cp ~{gcnv_rare_homozygous_variants} ~{destination_bucket}/~{sample_id}/
    >>>

    output {
        String gcnv_path_var_HQ_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.pathogenic_CNV_calls_high_quality.csv"
        String gcnv_path_var_HQ_non_clinical_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.pathogenic_CNV_calls_high_quality_non_clinical.csv"
        String gcnv_path_var_LQ_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.pathogenic_CNV_calls_low_quality.csv"
        String gcnv_path_var_LQ_non_clinical_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.pathogenic_CNV_calls_low_quality_non_clinical.csv"
        String gcnv_rare_homozygous_variants_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.rare_homozygous_CNV_calls.csv"
    }

    runtime {
        docker: "google/cloud-sdk"
        memory: "~{memory_gb}GB"
        disks: 'local-disk ~{disk_gb} HDD'
        preemptible: preemptible
    }
}

task variant_data_transfer_clinical {
    input {
        String sample_id
        String destination_bucket
        File gcnv_path_var_HQ
        File gcnv_path_var_HQ_non_clinical
        File gcnv_path_var_LQ
        File gcnv_path_var_LQ_non_clinical
        File gcnv_rare_homozygous_variants

        Int disk_gb = 16
        Int memory_gb = 8
        Int preemptible = 2

    }

    command <<<
        gsutil -m cp ~{gcnv_path_var_HQ} ~{destination_bucket}/
        gsutil -m cp ~{gcnv_path_var_HQ_non_clinical} ~{destination_bucket}/
        gsutil -m cp ~{gcnv_path_var_LQ} ~{destination_bucket}/
        gsutil -m cp ~{gcnv_path_var_LQ_non_clinical} ~{destination_bucket}/
        gsutil -m cp ~{gcnv_rare_homozygous_variants} ~{destination_bucket}/
    >>>

    output {
        String gcnv_path_var_HQ_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.pathogenic_CNV_calls_high_quality.csv"
        String gcnv_path_var_HQ_non_clinical_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.pathogenic_CNV_calls_high_quality_non_clinical.csv"
        String gcnv_path_var_LQ_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.pathogenic_CNV_calls_low_quality.csv"
        String gcnv_path_var_LQ_non_clinical_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.pathogenic_CNV_calls_low_quality_non_clinical.csv"
        String gcnv_rare_homozygous_variants_path = "~{destination_bucket}/~{sample_id}/~{sample_id}.rare_homozygous_CNV_calls.csv"
    }

    runtime {
        docker: "google/cloud-sdk"
        memory: "~{memory_gb}GB"
        disks: 'local-disk ~{disk_gb} HDD'
        preemptible: preemptible
    }
}