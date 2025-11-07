//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include "strelka_pos_processor.hh"
#include "strelka_streams.hh"
#include "strelka_run.hh"

#include "appstats/RunStatsManager.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "htsapi/bam_header_info.hh"
#include "htsapi/vcf_record_util.hh"
#include "starling_common/HtsMergeStreamerUtil.hh"
#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_pos_processor_util.hh"

// #undef USE_OPENMP

#if USE_OPENMP
#include <omp.h>
#endif



namespace INPUT_TYPE
{
enum index_t
{
    CANDIDATE_INDELS,
    FORCED_GT_VARIANTS,
    NOISE_VARIANTS,
    CALL_REGION
};
}



static
void
callRegion(
    const strelka_options& opt,
    const AnalysisRegionInfo& regionInfo,
    starling_read_counts& readCounts,
    reference_contig_segment& ref,
    HtsMergeStreamer& streamData,
    strelka_pos_processor& posProcessor)
{
    using namespace illumina::common;

    posProcessor.resetRegion(regionInfo.regionChrom, regionInfo.regionRange);
    streamData.resetRegion(regionInfo.streamerRegion.c_str());
    setRefSegment(opt, regionInfo.regionChrom, regionInfo.refRegionRange, ref);

    while (streamData.next())
    {
        const pos_t currentPos(streamData.getCurrentPos());
        const HTS_TYPE::index_t currentHtsType(streamData.getCurrentType());
        const unsigned currentIndex(streamData.getCurrentIndex());

        // wind posProcessor forward to position behind buffer head:
        posProcessor.set_head_pos(currentPos - 1);

        if (HTS_TYPE::BAM == currentHtsType)
        {
            // Remove the filter below because it's not valid for
            // RNA-Seq case, reads should be selected for the report
            // range by the bam reading functions
            //
            // /// get potential bounds of the read based only on current_pos:
            // const known_pos_range any_read_bounds(current_pos-maxIndelSize,current_pos+MAX_READ_SIZE+maxIndelSize);
            // if( posProcessor.is_range_outside_report_influence_zone(any_read_bounds) ) continue;

            // Approximate begin range filter: (removed for RNA-Seq)
            //if((current_pos+MAX_READ_SIZE+MAX_INDEL_SIZE) <= rlimit.begin_pos) continue;
            processInputReadAlignment(opt, ref, streamData.getCurrentBamStreamer(),
                                      streamData.getCurrentBam(), currentPos,
                                      readCounts, posProcessor, currentIndex);
        }
        else if (HTS_TYPE::VCF == currentHtsType)
        {
            assertExpectedVcfReference(ref, streamData.getCurrentVcfStreamer());
            const vcf_record& vcfRecord(streamData.getCurrentVcf());
            if (INPUT_TYPE::CANDIDATE_INDELS == currentIndex)     // process candidate indels input from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    process_candidate_indel(opt.maxIndelSize, vcfRecord, posProcessor);
                }
                else
                {
                    log_os << "WARNING: candidate indel vcf variant record cannot be categorized as indel:\n";
                    streamData.getCurrentVcfStreamer().report_state(log_os);
                }
            }
            else if (INPUT_TYPE::FORCED_GT_VARIANTS ==
                     currentIndex)     // process forced genotype tests from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    static const unsigned sample_no(0);
                    static const bool is_forced_output(true);
                    process_candidate_indel(opt.maxIndelSize, vcfRecord, posProcessor, sample_no, is_forced_output);
                }
                else if (vcfRecord.is_snv() or vcfRecord.is_ref_site())
                {
                    posProcessor.insert_forced_output_pos(vcfRecord.pos - 1);
                }
                else
                {
                    std::ostringstream oss;
                    oss << "forcedGT vcf variant record cannot be categorized as SNV or indel:\n";
                    streamData.getCurrentVcfStreamer().report_state(oss);
                    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
                }
            }
            else if (INPUT_TYPE::NOISE_VARIANTS == currentIndex)
            {
                if (vcfRecord.is_snv())
                {
                    SiteNoise sn;
                    set_noise_from_vcf(vcfRecord.line, sn);
                    posProcessor.insert_noise_pos(vcfRecord.pos - 1, sn);
                }
            }
            else
            {
                assert(false && "Unexpected hts index");
            }
        }
        else if (HTS_TYPE::BED == currentHtsType)
        {
            const bed_record& bedRecord(streamData.getCurrentBed());
            if (INPUT_TYPE::CALL_REGION == currentIndex)
            {
                known_pos_range2 range(bedRecord.begin,bedRecord.end);
                posProcessor.insertCallRegion(range);
            }
            else
            {
                assert(false && "Unexpected hts index");
            }
        }
        else
        {
            assert(false && "Invalid input condition");
        }
    }
}



void
strelka_run(
    const prog_info& pinfo,
    const strelka_options& opt)
{
    using namespace illumina::common;

    // ensure that this object is created first for runtime benchmark
    RunStatsManager statsManager(opt.segmentStatsFilename);

    opt.validate();

    const strelka_deriv_options dopt(opt);
    const StrelkaSampleSetSummary ssi;
    starling_read_counts readCounts;
    reference_contig_segment ref;

    ////////////////////////////////////////
    // setup streamData:
    //
    HtsMergeStreamer streamData(opt.referenceFilename);

    // additional data structures required in the region loop below, which are filled in as a side effect of
    // streamData initialization:
    std::vector<std::reference_wrapper<const bam_hdr_t>> bamHeaders;
    {
        std::vector<unsigned> registrationIndices;
        for (const bool isTumor : opt.alignFileOpt.isAlignmentTumor)
        {
            const unsigned rindex(isTumor ? STRELKA_SAMPLE_TYPE::TUMOR : STRELKA_SAMPLE_TYPE::NORMAL);
            registrationIndices.push_back(rindex);
        }

        bamHeaders = registerAlignments(opt.alignFileOpt.alignmentFilenames, registrationIndices, streamData);

        assert(not bamHeaders.empty());
        const bam_hdr_t& referenceHeader(bamHeaders.front());

        static const bool noRequireNormalized(false);
        registerVcfList(opt.input_candidate_indel_vcf, INPUT_TYPE::CANDIDATE_INDELS, referenceHeader, streamData,
                        noRequireNormalized);
        registerVcfList(opt.force_output_vcf, INPUT_TYPE::FORCED_GT_VARIANTS, referenceHeader, streamData);

        registerVcfList(opt.noise_vcf, INPUT_TYPE::NOISE_VARIANTS, referenceHeader, streamData);

        if (! opt.callRegionsBedFilename.empty())
        {
            streamData.registerBed(opt.callRegionsBedFilename.c_str(), INPUT_TYPE::CALL_REGION);
        }
    }

    const bam_hdr_t& referenceHeader(bamHeaders.front());
    const bam_header_info referenceHeaderInfo(referenceHeader);

    // parse and sanity check regions
    assert ((! opt.isHaplotypingEnabled) && "Region border size must be updated if haplotyping is enabled");
    const unsigned supplementalRegionBorderSize(opt.maxIndelSize);

    const auto& referenceAlignmentFilename(opt.alignFileOpt.alignmentFilenames.front());
    std::vector<AnalysisRegionInfo> regionInfoList;
    getStrelkaAnalysisRegions(opt, referenceAlignmentFilename, referenceHeaderInfo, supplementalRegionBorderSize,
                              regionInfoList);

    // Log region information
    log_os << "Total top-level regions to process: " << regionInfoList.size() << "\n";
    for (const auto& regionInfo : regionInfoList)
    {
        log_os << "Region: " << regionInfo.regionChrom 
               << ":" << (regionInfo.regionRange.begin_pos()+1) << "-" << regionInfo.regionRange.end_pos()
               << " (size: " << (regionInfo.regionRange.end_pos() - regionInfo.regionRange.begin_pos()) << " bp)\n";
    }

#if USE_OPENMP
    const int max_threads = omp_get_max_threads();
    const int num_regions = static_cast<int>(regionInfoList.size());
    if (max_threads > num_regions)
    {
        omp_set_num_threads(num_regions);
    }
    log_os << "Parallel processing: " << regionInfoList.size() 
           << " regions, max " << omp_get_max_threads() << " threads\n";
#else
    log_os << "Serial processing: " << regionInfoList.size() << " regions\n";
#endif

    // Process each top-level region in parallel (or serial if OpenMP not compiled)
    // Each iteration creates fresh data structures to avoid any state carryover
    #pragma omp parallel for schedule(dynamic)
    for (size_t regionIdx = 0; regionIdx < regionInfoList.size(); ++regionIdx)
    {
        const auto& regionInfo = regionInfoList[regionIdx];

        // Create local copy of options with region-specific output filenames
        // Use regionIdx instead of thread_id to avoid file overwriting when
        // a single thread processes multiple regions
        strelka_options opt_local(opt);
        
        // Add region-specific suffixes to output filenames
        const std::string region_suffix = ".region_" + std::to_string(regionIdx);
        
        if (! opt_local.somatic_snv_filename.empty())
        {
            opt_local.somatic_snv_filename += region_suffix;
        }
        if (! opt_local.somatic_indel_filename.empty())
        {
            opt_local.somatic_indel_filename += region_suffix;
        }
        if (! opt_local.somatic_callable_filename.empty())
        {
            opt_local.somatic_callable_filename += region_suffix;
        }
        if (! opt_local.segmentStatsFilename.empty())
        {
            opt_local.segmentStatsFilename += region_suffix;
        }
        if (! opt_local.realignedReadFilenamePrefix.empty())
        {
            opt_local.realignedReadFilenamePrefix += region_suffix;
        }

        // Create fresh data structures for this region (no reset, no state carryover)
        RunStatsManager statsManager_local(opt_local.segmentStatsFilename);
        starling_read_counts readCounts_local;
        reference_contig_segment ref_local;
        
        // Setup streamData
        HtsMergeStreamer streamData_local(opt_local.referenceFilename);
        std::vector<std::reference_wrapper<const bam_hdr_t>> bamHeaders_local;
        {
            std::vector<unsigned> registrationIndices;
            for (const bool isTumor : opt_local.alignFileOpt.isAlignmentTumor)
            {
                const unsigned rindex(isTumor ? STRELKA_SAMPLE_TYPE::TUMOR : STRELKA_SAMPLE_TYPE::NORMAL);
                registrationIndices.push_back(rindex);
            }

            bamHeaders_local = registerAlignments(opt_local.alignFileOpt.alignmentFilenames, registrationIndices, streamData_local);

            assert(not bamHeaders_local.empty());
            const bam_hdr_t& referenceHeader_local(bamHeaders_local.front());

            static const bool noRequireNormalized(false);
            registerVcfList(opt_local.input_candidate_indel_vcf, INPUT_TYPE::CANDIDATE_INDELS, referenceHeader_local, streamData_local,
                            noRequireNormalized);
            registerVcfList(opt_local.force_output_vcf, INPUT_TYPE::FORCED_GT_VARIANTS, referenceHeader_local, streamData_local);
            registerVcfList(opt_local.noise_vcf, INPUT_TYPE::NOISE_VARIANTS, referenceHeader_local, streamData_local);

            if (! opt_local.callRegionsBedFilename.empty())
            {
                streamData_local.registerBed(opt_local.callRegionsBedFilename.c_str(), INPUT_TYPE::CALL_REGION);
            }
        }

        const bam_hdr_t& referenceHeader_local(bamHeaders_local.front());
        const StrelkaSampleSetSummary ssi_local;
        strelka_streams fileStreams_local(opt_local, dopt, pinfo, referenceHeader_local, ssi_local);
        strelka_pos_processor posProcessor_local(opt_local, dopt, ref_local, fileStreams_local, statsManager_local);

        // Process this region (and its subregions if using BED file)
        if (not opt.isUseCallRegions())
        {
            // No BED file: process the single region
            callRegion(opt_local, regionInfo, readCounts_local, ref_local, streamData_local, posProcessor_local);
        }
        else
        {
            // With BED file: process subregions sequentially within this top-level region
            std::vector<known_pos_range2> subRegionRanges;
            getSubRegionsFromBedTrack(opt.callRegionsBedFilename, regionInfo.regionChrom, regionInfo.regionRange, subRegionRanges);

            for (const auto& subRegionRange : subRegionRanges)
            {
                AnalysisRegionInfo subRegionInfo;
                getStrelkaAnalysisRegionInfo(regionInfo.regionChrom, subRegionRange.begin_pos(), subRegionRange.end_pos(),
                                             supplementalRegionBorderSize, subRegionInfo);
                
                callRegion(opt_local, subRegionInfo, readCounts_local, ref_local, streamData_local, posProcessor_local);
            }
        }

        // Flush remaining variants for this region
        posProcessor_local.reset();
    }
}
