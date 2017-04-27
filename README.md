Supporting Perl scripts for the SigSeeker analysis platform

The benchmarking data compiled for the SigSeeker Bioinformatics paper 2017 was produced using
benchmark_cmd.pl and benchmarking.pl

The intersection of peak calling results for the ensemble is generated using updated_intersectbed.pl and SIG_CORR_PIPELINE_3D.pl

The file is called in the following fashion:
perl updated_intersectbed.pl [Input File 1] [Input File 2 ... Input File n] [Minimum Overlap in Basepairs] [Quantitative Threshold in Percent] [Output Filename] [ID]
