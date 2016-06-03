## readquant

This package contains helper functions for parsing both expression values
and technical features from common RNA-seq quantification tools.

The goal is to simplify going from a collection of quantification results to
expression table, as well as sample meta data.

A minimal example could look something like this:

    In [1]: from readquant import read_quants, read_qcs

    In [2]: tpm = read_quants('salmon/*E4_salmon_out', version='0.4.0')

    In [3]: tpm.head()
    Out[3]:
            salmon/1771-026-195-E4_salmon_out  \
    Name
    ERCC-00158                                0.0
    ERCC-00154                                0.0
    ERCC-00150                                0.0
    ERCC-00143                                0.0
    ERCC-00142                                0.0

            salmon/1771-026-190-E4_salmon_out  \
    Name
    ERCC-00158                                0.0
    ERCC-00154                                0.0
    ERCC-00150                                0.0
    ERCC-00143                                0.0
    ERCC-00142                                0.0

            salmon/1771-026-196-E4_salmon_out  \
    Name
    ERCC-00158                              0.000
    ERCC-00154                            906.236
    ERCC-00150                              0.000
    ERCC-00143                              0.000
    ERCC-00142                              0.000

            salmon/1771-026-198-E4_salmon_out  \
    Name
    ERCC-00158                                0.0
    ERCC-00154                                0.0
    ERCC-00150                                0.0
    ERCC-00143                                0.0
    ERCC-00142                                0.0

            salmon/1771-023-118-E4_salmon_out  \
    Name
    ERCC-00158                            0.00000
    ERCC-00154                            9.35988
    ERCC-00150                            0.00000
    ERCC-00143                          724.44800
    ERCC-00142                            0.00000

            salmon/1771-026-193-E4_salmon_out
    Name
    ERCC-00158                                0.0
    ERCC-00154                                0.0
    ERCC-00150                                0.0
    ERCC-00143                                0.0
    ERCC-00142                                0.0

    In [4]: sample_info = read_qcs('salmon/*_salmon_out', version='0.4.0', flen_lim=(10, 100))

    In [5]: sample_info.head()
    Out[5]:
                                   percent_mapped  num_processed  \
    salmon/1771-026-197-G4_salmon_out         71.8270      2665988.0
    salmon/1771-026-198-D2_salmon_out         67.7853      3841894.0
    salmon/1771-026-195-H7_salmon_out         73.4492      3875822.0
    salmon/1771-026-194-E9_salmon_out         47.5425       979446.0
    salmon/1771-026-195-E4_salmon_out         56.3955      5051277.0

                                   global_fl_mode  robust_fl_mode
    salmon/1771-026-197-G4_salmon_out           103.0           103.0
    salmon/1771-026-198-D2_salmon_out           122.0           122.0
    salmon/1771-026-195-H7_salmon_out           110.0           110.0
    salmon/1771-026-194-E9_salmon_out           111.0           111.0
    salmon/1771-026-195-E4_salmon_out           111.0           111.0
