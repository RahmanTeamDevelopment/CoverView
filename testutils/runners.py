import bamgen
import coverview.main
import coverview.reads
import json
import os
import pysam
import shutil
import uuid


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
            bed_file.write("{}\t{}\t{}\t{}\n".format(
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


class CoverViewTestRunner(object):
    """
    Utility class to wrap all the logic needed to run CoverView with simulated inputs
    """
    def __init__(self):
        self.bam_file_name = str(uuid.uuid4()) + ".bam"
        self.index_file_name = self.bam_file_name + ".bai"
        self.bed_file_name = self.bam_file_name.replace(".bam", ".bed")
        self.config_file_name = self.bam_file_name.replace(".bam", ".json")
        self.ref_file_name = "__MOCK__"
        self.read_sets = []
        self.regions = []
        self.config_data = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clean_up_input_files()
        self.clean_up_output_files()

    def add_reads(self, read_set):
        self.read_sets.append(read_set)

    def add_region(self, region):
        self.regions.append(region)

    def add_config_data(self, config_data):
        self.config_data.update(config_data)

    def generate_input_files(self):
        bamgen.make_bam_file(
            self.bam_file_name,
            self.read_sets
        )

        make_bed_file(
            self.bed_file_name,
            self.regions
        )

        make_config_file(
            self.config_file_name,
            self.config_data
        )

    def clean_up_input_files(self):
        os.remove(self.bam_file_name)
        os.remove(self.index_file_name)
        os.remove(self.bed_file_name)
        os.remove(self.config_file_name)

    def clean_up_output_files(self):
        os.remove("output_regions.txt")
        os.remove("output_profiles.txt")
        os.remove("output_summary.txt")

        output_config = self.config_data.get("outputs")

        if output_config is not None:
            has_gui_output = output_config.get(
                "gui",
                False
            )
        else:
            has_gui_output = False

        if has_gui_output:
            gui_output_dir = self.config_data.get("outputs").get(
                "gui_output_directory",
                os.path.join(os.getcwd(), "output_gui_data")
            )

            shutil.rmtree(gui_output_dir)

    def run_coverview_and_get_exit_code(self):
        self.generate_input_files()

        command_line_args = make_command_line_arguments(
            bam_file_name=self.bam_file_name,
            bed_file_name=self.bed_file_name,
            reference_file_name=self.ref_file_name,
            config_file_name=self.config_file_name
        )

        return coverview.main.main(command_line_args)
