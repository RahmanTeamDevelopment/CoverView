import bamgen
import coverview.main
import json
import os
import pysam
import shutil
import unittest
import uuid


def make_command_line_arguments(bam_file_name, bed_file_name, config_file_name):
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


class TestCoverViewWithGuiOutput(unittest.TestCase):
    """
    Testing that we can run CoverView and generate the GUI output correctly.
    """
    def setUp(self):
        self.unique_bam_file_name = str(uuid.uuid4()) + ".bam"
        self.unique_index_file_name = self.unique_bam_file_name + ".bai"
        self.unique_bed_file_name = self.unique_bam_file_name.replace(".bam", ".bed")
        self.unique_config_file_name = self.unique_bam_file_name.replace(".bam", ".json")

    def tearDown(self):
        os.remove(self.unique_bam_file_name)
        os.remove(self.unique_index_file_name)
        os.remove(self.unique_bed_file_name)
        os.remove(self.unique_config_file_name)
        os.remove("output_regions.txt")
        os.remove("output_profiles.txt")
        os.remove("output_summary.txt")
        shutil.rmtree("output_gui")

    def test_coverview_runs_and_returns_0_exit_code(self):
        read_sets = [
            ("1", 32, 100, 0)
        ]

        regions = [
            ( "1", 32, 132, "Region_1")
        ]

        coverview_dir = os.getcwd()
        gui_dir = os.path.join(coverview_dir, "gui")

        config = \
            {
                "outputs": {
                    "gui": True,
                    "gui_output_directory": os.path.join(coverview_dir, "output_gui")
                },

                "gui": {
                    "template_gui_html_file": os.path.join(gui_dir, "gui.html"),
                    "javascript_directory": os.path.join(gui_dir, "lib")
                },

                "reference_file": "__MOCK__"
            }


        bamgen.make_bam_file(self.unique_bam_file_name, read_sets)
        make_bed_file(self.unique_bed_file_name, regions)
        make_config_file(self.unique_config_file_name, config)

        command_line_args = make_command_line_arguments(
            self.unique_bam_file_name,
            self.unique_bed_file_name,
            self.unique_config_file_name
        )

        assert coverview.main.main(command_line_args) == 0


if __name__ == "__main__":
    unittest.main()