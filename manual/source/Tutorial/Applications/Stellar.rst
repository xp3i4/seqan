.. sidebar:: ToC

    .. contents::

.. _tutorial-apps-mason:

STELLAR
======

Overview
--------

STELLAR is a tool for pairwise local alignments. It has full sensitivity for :math:`{\epsilon}`-alignements with a maximal error rate and a minimal length. It implements SWIFT filter algorithm and is very practical and fast on very long sequences.   

Algorithm
---------
**Filtration**

    Filter in STELLAR applies SWIFT algorithm (Rasmussen et al. 2006) to find out potential regions.The SWIFT algorithm provides a necessary condition that for any given error rate :math:`\epsilon` and alignment length :math:`n_0` there exist :math:`\omega, q, e` and :math:`\tau` such that every :math:`\epsilon`-match contains :math:`\tau` q-hits that reside in a :math:`\omega × e` parallelogram. A :math:`\omega × e` parallelogram is the intersection of :math:`e + 1` consecutive diagonals and :math:`w + 1` consecutive columns in the dotplot. Where :math:`q,\tau,\omega,e` meets the following two requirments

.. math:: 
   q < \lceil 1 / \epsilon \rceil \quad and \quad \tau \leq min\lbrace U(n_0, q, \epsilon), U(n_1, q, \epsilon)\rbrace 

   \omega = (\tau - 1) + q(e + 1) \quad and \quad e = \lfloor \frac{2(\tau - 1) + (q -1 )}{1/\epsilon - q} \rfloor\quad 
..

   where :math:`{n_1=\lceil(\lfloor\epsilon n_0\rfloor + 1) / \epsilon \rceil \quad and \quad U(n, q, \epsilon) = (n+1)-q(\lfloor \epsilon n \rfloor + 1))}`

   Based on that STELLAR detects :math:`\omega × e` parallelograms with :math:`\tau` q-hits in the dotplot. It slides from left to right over one sequence and searches overlapping q-grams in a q-gram index of the other sequence. Parallelograms whose q-gram counter has reached :math:`\tau` is output as a potential alignment region.

.. image:: ./stellar1.jpg  
   :width: 300px
   :align: right
.. 

  The following figure shows an example for SWIFT hits containing either a subalignment of an :math:`\epsilon`-match, whole :math:`\epsilon`-matches, no :math:`\epsilon`-match or an :math:`\epsilon`-match with an X-drop, where :math:`n_0` = 20, :math:`\epsilon` = 0.1, and :math:`q` = 6. Accordingly, :math:`\omega` =20, :math:`\tau` = 3, and :math:`e` = 2. SWIFT searches parallelograms that contain at least τ = 3 q-gram hits by streaming over sequence 1 and searching common q-grams in sequence 2. Subfigure (a) shows an ε-match that results in two SWIFT hits and the :math:`\epsilon`-match is longer than both of the two hits. (b) shows a SWIFT hit that contains two :math:`epsilon`-matches and (c) shows a false positive SWIFT hit induced by three separated q-gram hits. (d) shows a SWIFT hit that contains an :math:`\epsilon`-match with a 3-drop.

**Verification**

    In the verification part, STELLAR will first find all possible :math:`{\epsilon}`-core by applying a banded version of the Waterman-Eggert local alignment algorithm [31], where :math:`{\epsilon}`-core is a segment of an :math:`\epsilon`-match and its match score is bigger than a given minimal value. Then it filters out the ones containing :math:`{\epsilon}`-X-drops. Then it spans those :math:`{\epsilon}`-core from the last step to maximal length. After that it finds the longest extension as a possible :math:`{\epsilon}`-match. At last it remove overlapped :math:`{\epsilon}`-match and choose the maximal-:math:`{\epsilon}`-match as the final :math:`{\epsilon}`-match.


Program options
---------------

**Usage**

.. code-block:: bash

  stellar [Options] <FASTA sequence file 1> <FASTA sequence file 2>

..

  <FASTA sequence file 1>
 
  Set the name of the file containing the database sequence(s) in
  (multi-)Fasta format.  Note that STELLAR expects Fasta identifiers to
  be unique before the first whitespace.
 
  <FASTA sequence file 2>
 
  Set the name of the file containing the query sequence(s) in (multi-)
  Fasta format. All query sequences will be compared to all database
  sequences.  Note that STELLAR expects Fasta identifiers to be unique
  before the first whitespace.

  Without any additional parameters STELLAR compares the query sequence(s)
  to both strands of the database sequence(s) with an error rate of 0.05
  (i.e. 5% errors, an identity of 95%) and a minimal length of 100, and
  dumps all local alignments in an output file named "stellar.gff".

**Main options**
  [ -e NUM ],  [ --epsilon NUM ]
  
  Set the maximal error rate for local alignments. NUM must be a floating
  point number between 0 and 0.25. The default value is 0.05. NUM is the
  number of edit operations needed to transform the aligned query substring
  to the aligned database substring divided by the length of the local
  alignment. For example specify '-e 0.1' for a maximal error rate of 10%
  (90% identity).

  [ -l NUM ],  [ --minLength NUM ]
  
  Set the minimal length of a local alignment. The default value is 100.
  NUM is the minimal length of an alignment, not the length of the
  aligned substrings.

  [ -f ],  [ --forward ]

  Only compare query sequence(s) to the positive/forward strand of the
  database sequence(s). By default, both strands are scanned.

  [ -r ],  [ --reverse ]

  Only compare query sequence(s) to the negative/reverse-complemented 
  database sequence(s). By default, both strands are scanned.

  [ -a ],  [ --alphabet STR ]

  Select the alphabet type of your input sequences. Valid values are dna,
  dna5, rna, rna5, protein, and char. Choose dna5 and rna5 for input
  sequences containing N characters, choose protein for amino acid
  sequences. By default, dna5 is assumed.

  [ -v ],  [ --verbose ]
  
  Verbose. Print extra information and running times.

Usage examples
--------------
In this example we selected chromosome of D. melanogaster and D. pseudoobscura for our test runs. 

Input files:
    Sequences data came from FlyBase

    1. Chromosome arm 2L from D. melanogaster: `dmel-2L-chromosome-r5.26.fasta dpse-4_group3-r2.14.fasta`_
.. _dmel-2L-chromosome-r5.26.fasta dpse-4_group3-r2.14.fasta: ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.16_FB2009_03/fasta/

..    

    2. Group 3 from the chromosome 4 assembly of D. pseudoobscura: `dpse-4_group3-r2.14.fasta`_ extracted from dpse-all-chromosome-r2.14.fasta
.. _dpse-4_group3-r2.14.fasta: ftp://ftp.flybase.net/12_species_analysis/genomes/Drosophila_pseudoobscura/dpse_r2.14_FB2010_08/fasta/

..    

Test: 
    with minimal length = 200 and default maximal error rate = 0.05

.. code-block:: bash 

    ./stellar -l 200 dmel-2L-chromosome-r5.26.fasta dpse-4_group3-r2.14.fasta 

..
    
Output:
    The output looks like the following console. There are altogether 44 matches found by STELLAR in test one. The alignment result is written into a .gff file as follow

.. code-block:: console

    I/O options:
      database file   : dmel-2L-chromosome-r5.26.fasta
      query file      : dpse-4_group3-r2.14.fasta
      alphabet        : dna5
      output file     : stellar.gff
      output format   : gff

    User specified parameters:
      minimal match length             : 200
      maximal error rate (epsilon)     : 0.05
      maximal x-drop                   : 5
      search forward strand            : yes
      search reverse complement        : yes

      verification strategy            : exact
      maximal number of matches        : 50
      duplicate removal every          : 500

    Calculated parameters:
      k-mer length : 18
      s^min        : 18
      threshold    : 3
      distance cut : 200
      delta        : 16
      overlap      : 10

    Loaded 1 query sequence.
    Loaded 1 database sequence.

    All matches matches resulting from your search have an E-value of: 
            1.26922e-74 or smaller  (match score = 1, error penalty = -2)

    Constructing index...

    Aligning all query sequences to database sequence...
      2L type=chromosome_arm; loc=2L:1..23011544; ID=2L; dbxref=REFSEQ:NT_033779,GB:AE014134; MD5=bfdfb99d39fa5174dae1e2ecd8a231cd; length=23011544; release=r5.26; species=Dmel;
      2L type=chromosome_arm; loc=2L:1..23011544; ID=2L; dbxref=REFSEQ:NT_033779,GB:AE014134; MD5=bfdfb99d39fa5174dae1e2ecd8a231cd; length=23011544; release=r5.26; species=Dmel;, complement

    # Eps-matches     : 44
..

.. includefrags:: manual/source/Tutorial/Applications/stellar.gff 
..

References
----------
