#!/usr/bin/env python3
"""
Create IGV session file from FLAIR evaluation outputs.
Based on nf-core/atacseq igv_files_to_session.py
"""

import argparse
import os
from pathlib import Path


def create_igv_session(output_xml, genome, bam_file=None, bed_files=None, bigwig_files=None, gtf_file=None):
    """
    Create IGV session XML file.
    
    Args:
        output_xml: Path to output XML file
        genome: Genome build (e.g., 'hg38')
        bam_file: BAM file path (optional)
        bed_files: List of BED file paths (optional)
        bigwig_files: List of BigWig file paths (optional)
        gtf_file: GTF file path (optional)
    """
    # Initialize lists
    bed_files = bed_files or []
    bigwig_files = bigwig_files or []
    
    # Filter out placeholder files and non-existent files
    if bam_file:
        if not os.path.exists(bam_file) or os.path.getsize(bam_file) == 0:
            bam_file = None
    
    bed_files = [f for f in bed_files if f and os.path.exists(f) and os.path.getsize(f) > 0]
    bigwig_files = [f for f in bigwig_files if f and os.path.exists(f) and os.path.getsize(f) > 0]
    
    if gtf_file and (not os.path.exists(gtf_file) or os.path.getsize(gtf_file) == 0):
        gtf_file = None
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_xml) or '.', exist_ok=True)
    
    # Start building XML string
    xml_str = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    xml_str += f'<Session genome="{genome}" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">\n'
    
    # ADD RESOURCES SECTION
    xml_str += '\t<Resources>\n'
    
    if bam_file:
        xml_str += f'\t\t<Resource path="{bam_file}"/>\n'
    
    for bw in bigwig_files:
        xml_str += f'\t\t<Resource path="{bw}"/>\n'
    
    for bed in bed_files:
        xml_str += f'\t\t<Resource path="{bed}"/>\n'
    
    if gtf_file:
        xml_str += f'\t\t<Resource path="{gtf_file}"/>\n'
    
    # Add reference track
    xml_str += '\t\t<Resource name="RefSeq All" path="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz" type="refgene"/>\n'
    
    xml_str += '\t</Resources>\n'
    
    # ADD DATA PANEL WITH TRACKS
    xml_str += '\t<Panel height="1160" name="DataPanel" width="1897">\n'
    
    # Add BAM tracks (coverage, junctions, and alignments)
    if bam_file:
        bam_name = os.path.basename(bam_file)
        
        # Coverage track
        xml_str += f'\t\t<Track attributeKey="{bam_name} Coverage" autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" fontSize="10" id="{bam_file}_coverage" name="{bam_name} Coverage" snpThreshold="0.2" visible="true">\n'
        xml_str += '\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>\n'
        xml_str += '\t\t</Track>\n'
        
        # Junction track
        xml_str += f'\t\t<Track attributeKey="{bam_name} Junctions" autoScale="true" clazz="org.broad.igv.sam.SpliceJunctionTrack" fontSize="10" groupByStrand="false" height="60" id="{bam_file}_junctions" maxdepth="1" name="{bam_name} Junctions" visible="true"/>\n'
        
        # Alignment track
        xml_str += f'\t\t<Track attributeKey="{bam_name}" clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" displayMode="SQUISHED" experimentType="RNA" fontSize="10" id="{bam_file}" name="{bam_name}" visible="true">\n'
        xml_str += '\t\t\t<RenderOptions colorOption="READ_STRAND" duplicatesOption="FILTER"/>\n'
        xml_str += '\t\t</Track>\n'
    
    # Add BigWig tracks
    for bw in bigwig_files:
        bw_name = os.path.basename(bw)
        xml_str += f'\t\t<Track attributeKey="{bw_name}" autoScale="true" clazz="org.broad.igv.track.DataSourceTrack" fontSize="10" id="{bw}" name="{bw_name}" renderer="BAR_CHART" visible="true" windowFunction="mean">\n'
        xml_str += '\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="1.0" minimum="0.0" type="LINEAR"/>\n'
        xml_str += '\t\t</Track>\n'
    
    xml_str += '\t</Panel>\n'
    
    # ADD BED PANEL if there are BED files
    if bed_files:
        xml_str += '\t<Panel height="100" name="BedPanel" width="1897">\n'
        for bed in bed_files:
            bed_name = os.path.basename(bed)
            xml_str += f'\t\t<Track attributeKey="{bed_name}" clazz="org.broad.igv.track.FeatureTrack" colorScale="ContinuousColorScale;0.0;1164.0;255,255,255;0,0,178" fontSize="10" groupByStrand="false" id="{bed}" name="{bed_name}" visible="true"/>\n'
        xml_str += '\t</Panel>\n'
    
    # ADD FEATURE PANEL with reference tracks
    xml_str += '\t<Panel height="200" name="FeaturePanel" width="1897">\n'
    xml_str += '\t\t<Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" sequenceTranslationStrandValue="+" shouldShowTranslation="false" visible="true"/>\n'
    
    if gtf_file:
        gtf_name = os.path.basename(gtf_file)
        xml_str += f'\t\t<Track attributeKey="{gtf_name}" clazz="org.broad.igv.track.FeatureTrack" displayMode="COLLAPSED" fontSize="10" groupByStrand="false" id="{gtf_file}" name="{gtf_name}" visible="true"/>\n'
    
    xml_str += '\t\t<Track attributeKey="RefSeq All" clazz="org.broad.igv.track.FeatureTrack" displayMode="COLLAPSED" fontSize="10" groupByStrand="false" id="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz" name="RefSeq All" visible="true"/>\n'
    xml_str += '\t</Panel>\n'
    
    # Panel layout fractions
    if bed_files:
        xml_str += '\t<PanelLayout dividerFractions="0.8,0.9"/>\n'
    else:
        xml_str += '\t<PanelLayout dividerFractions="0.85"/>\n'
    
    # Hidden attributes
    xml_str += '\t<HiddenAttributes>\n'
    xml_str += '\t\t<Attribute name="DATA FILE"/>\n'
    xml_str += '\t\t<Attribute name="DATA TYPE"/>\n'
    xml_str += '\t\t<Attribute name="NAME"/>\n'
    xml_str += '\t</HiddenAttributes>\n'
    
    xml_str += '</Session>'
    
    # Write to file
    with open(output_xml, 'w') as f:
        f.write(xml_str)
    
    print(f"Created IGV session file: {output_xml}")
    print(f"  BAM file: {'Yes' if bam_file else 'No'}")
    print(f"  BED files: {len(bed_files)}")
    print(f"  BigWig files: {len(bigwig_files)}")
    print(f"  GTF file: {'Yes' if gtf_file else 'No'}")


def main():
    parser = argparse.ArgumentParser(
        description='Create IGV session file from FLAIR evaluation outputs.',
        epilog='Example: python prepare_igv_files.py session.xml --genome hg38 --bam reads.bam --bed isoforms.bed --bed cage_peaks.bed'
    )
    parser.add_argument('output_xml', help='Output XML session file')
    parser.add_argument('--genome', default='hg38', help='Genome build (default: hg38)')
    parser.add_argument('--bam', help='BAM alignment file')
    parser.add_argument('--bed', action='append', dest='bed_files', help='BED file(s) to include (can specify multiple times)')
    parser.add_argument('--bigwig', action='append', dest='bigwig_files', help='BigWig file(s) to include (can specify multiple times)')
    parser.add_argument('--gtf', help='GTF annotation file')
    
    args = parser.parse_args()
    
    create_igv_session(
        output_xml=args.output_xml,
        genome=args.genome,
        bam_file=args.bam,
        bed_files=args.bed_files,
        bigwig_files=args.bigwig_files,
        gtf_file=args.gtf
    )


if __name__ == '__main__':
    main()
