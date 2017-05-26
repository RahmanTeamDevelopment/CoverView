import bamgen
import json
import pysam


def make_command_line_arguments(bam_file_name,
                                bed_file_name,
                                reference_file_name,
                                config_file_name):
    """
    Utility function to construct a list of command-line arguments in the form that is stored
    in sys.argv. This can then be fed to the main function and used to run CoverView.
    """
    arguments_list = [
        "-i",
        bam_file_name
    ]

    if bed_file_name is not None:
        arguments_list.extend([
            "-b",
            bed_file_name
        ])

    if reference_file_name is not None:
        arguments_list.extend([
            "-r",
            reference_file_name
        ])

    if config_file_name is not None:
        arguments_list.extend([
            "-c",
            config_file_name
        ])

    return arguments_list


def make_bam_file(file_name, read_sets):
    """
    Utility function to create a BAM file with the specified configuration. We create a mock
    reference with the union of all references in the read_sets and each reference sequence is
    at least as long as the max(2 * (start_position + read_length)), which is longer than we need but
    that's ok.
    """
    ref_file = bamgen.MockReferenceFile()

    bamgen.generate_bam_file(
        file_name,
        ref_file,
        read_sets
    )


def make_bed_file(file_name, regions):
    """
    Output a BED file containing a list of genomic intervals    
    """
    with open(file_name, 'w') as bed_file:
        for chrom, start, end, name in regions:
            bed_file.write("\t{}\t{}\t{}\t{}".format(
                chrom, start, end, name
            ))


def make_config_file(file_name, config_dict):
    with open(file_name, 'w') as config_file:
        json.dump(config_dict, config_file)


def load_bam_into_read_array(file_name):
    """
    Utility function for creating a read array from the contents of a sorted BAM file    
    """
    read_array = coverview.reads.pyReadArray()

    with pysam.AlignmentFile(file_name, 'rb') as bam_file:
        for read in bam_file:
            read_array.append(read)

    return read_array