How to run the example:
1. Compile ``find_repeats_chr2_232_4Mb.cpp`` and ``find_repeats_chr19_14_8Mb.cpp`` from the cpp_files directory (see comment at top of each file for information on the g++ command to run), and put the binary executables in the working directory.
2. Symlink the GRCh38 reference .fa and .fa.fai files to the working directory (or edit the path to the GRCh38 reference in ``compute_mid_length_mean_somatic_expansion_1000G_GBR.sh``).
3. Run the example analysis: ``bash compute_mid_length_mean_somatic_expansion_1000G_GBR.sh``
4. Reads extracted for analysis will be stored in the ``chr2_232_4Mb/reads/`` and ``chr19_14_8Mb/reads/`` directories, and the two example output files will be generated.
