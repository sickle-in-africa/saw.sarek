#!/usr/bin/env nextflow

test_ch = Channel
	.of("hello world!")

process printMessge {
	echo true

	input:
		val(message) from test_ch

	script:
	aligner = "bwa"
    """
    echo "$message"
    echo "${params.sarekDir}"
    """
}
