mkdir: cannot create directory ‘/mnt/causes-vnx1/Pipelines/AnnotateVariants/Test/tmpdir/’: File exists
normalize v0.5

options:     input VCF file                                  -
         [o] output VCF file                                 -
         [w] sorting window size                             10000
         [m] no fail on masked reference inconsistency       false
         [n] no fail on reference inconsistency              false
         [q] quiet                                           false
         [d] debug                                           false
         [r] reference FASTA file                            /mnt/causes-data01/data/GENOMES/hg19/FASTA/hg19.fa

decompose v0.5

options:     input VCF file        -
         [s] smart decomposition   true (experimental)
         [o] output VCF file       -


stats: no. variants                 : 146217
       no. biallelic variants       : 143319
       no. multiallelic variants    : 2898

       no. additional biallelics    : 3379
       total no. of biallelics      : 149596

Time elapsed: 1.53s


stats: biallelic
          no. left trimmed                      : 0
          no. right trimmed                     : 1967
          no. left and right trimmed            : 0
          no. right trimmed and left aligned    : 0
          no. left aligned                      : 0

       total no. biallelic normalized           : 1967

       multiallelic
          no. left trimmed                      : 0
          no. right trimmed                     : 0
          no. left and right trimmed            : 0
          no. right trimmed and left aligned    : 0
          no. left aligned                      : 0

       total no. multiallelic normalized        : 0

       total no. variants normalized            : 1967
       total no. variants observed              : 149596
       total no. reference observed             : 0

Time elapsed: 2.31s



NEW VERSION!
	There is a new SnpEff version available: 
		Version      : 4.3Q
		Release date : 2018-01-04
		Download URL : http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip


=============================================
vcfanno version 0.2.8 [built with go1.8]

see: https://github.com/brentp/vcfanno
=============================================
vcfanno.go:115: found 72 sources from 25 files
vcfanno.go:241: annotated 149596 variants in 278.81 seconds (536.5 / second)
/opt/tools/miniconda/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:219: SAWarning: Unicode type received non-unicode bind param value 'NA12878_Trio'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
/opt/tools/miniconda/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:219: SAWarning: Unicode type received non-unicode bind param value '-9'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
/opt/tools/miniconda/lib/python2.7/site-packages/sqlalchemy/sql/sqltypes.py:219: SAWarning: Unicode type received non-unicode bind param value '2'. (this warning may be suppressed after 10 occurrences)
  (util.ellipses_string(value),))
10000 variant_impacts:31898	effects time: 6.4	chunk time:20.8	481.68 variants/second	expanded columns:2.62
20000 variant_impacts:20041	effects time: 4.5	chunk time:18.5	509.23 variants/second	expanded columns:4.50
30000 variant_impacts:28359	effects time: 5.5	chunk time:21.6	492.46 variants/second	expanded columns:6.06
40000 variant_impacts:31212	effects time: 5.8	chunk time:23.6	473.03 variants/second	expanded columns:7.60
50000 variant_impacts:30696	effects time: 6.0	chunk time:25.4	454.81 variants/second	expanded columns:9.18
60000 variant_impacts:18603	effects time: 4.3	chunk time:24.2	447.32 variants/second	expanded columns:10.59
70000 variant_impacts:28286	effects time: 5.4	chunk time:26.2	436.69 variants/second	expanded columns:10.86
80000 variant_impacts:19010	effects time: 4.3	chunk time:25.2	431.24 variants/second	expanded columns:11.97
90000 variant_impacts:39563	effects time: 7.4	chunk time:33.6	410.82 variants/second	expanded columns:14.48
100000 variant_impacts:31915	effects time: 6.0	chunk time:30.6	400.45 variants/second	expanded columns:14.15
indexing ... finished in 5.5 seconds...
total time: in 416.0 seconds...
