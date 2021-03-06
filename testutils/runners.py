import bamgen.bamgen
import coverview_.main
import coverview_.reads
import coverview_.transcript
import json
import os
import pysam
import uuid


def remove_if_exists(file_name):
    if file_name is None:
        return

    if os.path.exists(file_name):
        os.remove(file_name)


def make_command_line_arguments(bam_file_name,
                                bed_file_name,
                                config_file_name,
                                transcript_file_name,
                                gui_output_file_name):
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

    if config_file_name is not None:
        arguments_list.extend([
            "-c",
            config_file_name
        ])

    if transcript_file_name is not None:
        arguments_list.extend([
            "-t",
            transcript_file_name
        ])

    if gui_output_file_name is not None:
        arguments_list.extend([
            "--gui_json_output_file",
            gui_output_file_name
        ])

    return arguments_list


def make_bam_file(file_name, read_sets):
    """
    Utility function to create a BAM file with the specified configuration. We create a mock
    reference with the union of all references in the read_sets and each reference sequence is
    at least as long as the max(2 * (start_position + read_length)), which is longer than we need but
    that's ok.
    """
    ref_file = bamgen.bamgen.MockReferenceFile()

    bamgen.bamgen.generate_bam_file(
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
            bed_file.write("{}\t{}\t{}\t{}\n".format(
                chrom, start, end, name
            ))


def make_transcripts_file(
        transcripts,
        transcript_file_name,
):
    coverview_.transcript.write_transcripts_to_indexed_tabix_file(
        transcripts,
        transcript_file_name
    )


def make_config_file(file_name, config_dict):
    with open(file_name, 'w') as config_file:
        json.dump(config_dict, config_file)


def load_bam_into_read_array(file_name):
    """
    Utility function for creating a read array from the contents of a sorted BAM file
    """
    read_array = coverview_.reads.pyReadArray()

    with pysam.AlignmentFile(file_name, 'rb') as bam_file:
        for read in bam_file:
            read_array.append(read)

    return read_array


class CoverViewTestRunner(object):
    """
    Utility class to wrap all the logic needed to run CoverView with simulated inputs
    """
    def __init__(self):
        self.bam_file_name = str(uuid.uuid4()) + ".bam"
        self.index_file_name = self.bam_file_name + ".bai"
        self.bed_file_name = self.bam_file_name.replace(".bam", ".bed")
        self.config_file_name = self.bam_file_name.replace(".bam", ".json")
        self.transcript_file_name = self.bam_file_name.replace(".bam", "_transcript_db.txt")
        self.compressed_transcript_file_name = self.bam_file_name.replace(".bam", "_transcript_db.txt.gz")
        self.transcript_file_index_name = self.bam_file_name.replace(".bam", "_transcript_db.txt.gz.tbi")
        self.gui_output_file_name = None
        self.read_sets = []
        self.regions = []
        self.transcripts = []
        self.config_data = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clean_up_input_files()
        self.clean_up_output_files()

    def add_reads(self, read_set):
        self.read_sets.append(read_set)

    def add_unmapped_reads(self, read_length, number_of_reads):
        self.read_sets.append((None, None, read_length, number_of_reads))

    def add_region(self, region):
        self.regions.append(region)

    def add_config_data(self, config_data):
        self.config_data.update(config_data)

    def add_transcript(self, transcript):
        self.transcripts.append(transcript)

    def add_gui_output_file(self, file_name):
        self.gui_output_file_name = file_name

    def generate_input_files(self):
        bamgen.bamgen.make_bam_file(
            self.bam_file_name,
            self.read_sets
        )

        if len(self.regions) > 0:
            make_bed_file(
                self.bed_file_name,
                self.regions
            )
        else:
            self.bed_file_name = None

        if len(self.transcripts) > 0:
            make_transcripts_file(
                self.transcripts,
                self.transcript_file_name,
            )
        else:
            self.transcript_file_name = None
            self.compressed_transcript_file_name = None
            self.transcript_file_index_name = None

        make_config_file(
            self.config_file_name,
            self.config_data
        )

    def clean_up_input_files(self):
        remove_if_exists(self.bam_file_name)
        remove_if_exists(self.bed_file_name)
        remove_if_exists(self.compressed_transcript_file_name)
        remove_if_exists(self.config_file_name)
        remove_if_exists(self.index_file_name)
        remove_if_exists(self.transcript_file_index_name)
        remove_if_exists(self.transcript_file_name)

    def clean_up_output_files(self):
        remove_if_exists("output_regions.txt")
        remove_if_exists("output_profiles.txt")
        remove_if_exists("output_summary.txt")
        remove_if_exists("output_poor.txt")
        remove_if_exists(self.gui_output_file_name)

    def run_coverview_and_get_exit_code(self):
        self.generate_input_files()

        command_line_args = make_command_line_arguments(
            bam_file_name=self.bam_file_name,
            bed_file_name=self.bed_file_name,
            config_file_name=self.config_file_name,
            transcript_file_name=self.compressed_transcript_file_name,
            gui_output_file_name=self.gui_output_file_name
        )

        return coverview_.main.main(command_line_args)
