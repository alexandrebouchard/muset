package muset;

import tutorialj.Tutorial;

public class Doc
{
  /**
   * 
   * Summary
   * -------
   * 
   * Muset is a library for manipulating and inferring multiple sequence alignments.
   * 
   * The main features are:
   * 
   * 1. MSAPoset, an efficient poset representation of MSAs, based on http://bioinformatics.oxfordjournals.org/content/23/2/e24.full
   * 2. Efficient and flexible pairwise alignment (max, sum, and sampling), based on alignment random fields
   *   http://www.eecs.berkeley.edu/Pubs/TechRpts/2010/EECS-2010-153.pdf 
   * 3. Variational inference algorithms for MSAs 
   *   http://papers.nips.cc/paper/4036-variational-inference-over-combinatorial-spaces.pdf
   *  
   * Note: as of Apr 6 2015, 1 and 2 are operational, while code for 3 is being transferred from a legacy SVN library.
   * 
   * Muset stands for MUltiple SEquence Toolkit.
   * 
   * 
   * Installation
   * ------------
   * 
   * There are two ways to install:
   * 
   * ### Integrate to a gradle script
   * 
   * Simply add the following lines (replacing 1.0.0 by the current version (see git tags)):
   * 
   * ```groovy
   * repositories {
   *  mavenCentral()
   *  jcenter()
   *  maven {
   *     url "http://www.stat.ubc.ca/~bouchard/maven/"
   *   }
   * }
   * 
   * dependencies {
   *   compile group: 'ca.ubc.stat', name: 'muset', version: '1.0.0'
   * }
   * ```
   * 
   * ### Compile using the provided gradle script
   * 
   * - Check out the source ``git clone git@github.com:alexandrebouchard/muset.git``
   * - Compile using ``gradle installApp``
   * - Add the jars in ``build/install/muset/lib/`` into your classpath
   * 
   * ### Use in eclipse
   * 
   * - Check out the source ``git clone git@github.com:alexandrebouchard/muset.git``
   * - Type ``gradle eclipse`` from the root of the repository
   * - From eclipse:
   *   - ``Import`` in ``File`` menu
   *   - ``Import existing projects into workspace``
   *   - Select the root
   *   - Deselect ``Copy projects into workspace`` to avoid having duplicates
   */
  @Tutorial(startTutorial = "README.md", showSource = false)
  public void installInstructions()
  {
  }
  
  /**
   * Using ``MSAPoset``
   * ------------------
   */
  @Tutorial(startTutorial = "README.md", showSource = false, nextStep = MSAPosetTest.class)
  public void msaPosetTest() {}
  
  /**
   * Using ``Aligner``
   * -----------------
   * 
   * ``Aligner`` is a main function with two purposes:
   * 
   * 1. calling the machinery in ``muset.pef`` to train alignment weights 
   * based on ``ExponentialFamily`` (i.e. GLM of maxent model for alignments, as in Chapter 
   * Chapter 5 of
   * [http://www.eecs.berkeley.edu/Pubs/TechRpts/2010/EECS-2010-153.pdf](http://www.eecs.berkeley.edu/Pubs/TechRpts/2010/EECS-2010-153.pdf)
   * 2. creating alignments by using ``MSAPoset.consensus`` computed from the weights in (1)
   * 
   * #### Preparing the input
   * 
   * The program takes as input a dataset formatted as follows (see ``src/test/resources/testdataset.csv``):
   * 
   * ```
   * homologous group 1,taxon1,firstLetter secondLetter
   * homologous group 1,taxon2,thirdLetter firstLetter
   * "homologous group, number 2",taxon3,secondLetter secondLetter
   * "homologous group, number 2",taxon1,secondLetter firstLetter
   * ```
   * 
   * Explanations of the format:
   * 
   * - The first column is an identifier for a group of homologous sequences (cognate id or protein family)
   * - The second column is an identifier for the taxon (language, species, dialect, operational taxonomic unit)
   * - The third column contains space separated letters (note that these 'letters' or character do not need to be
   *   one letter long, or even have the same length.
   * - Note that the format is fairly rigid beyond what is demonstrated here (avoid empty lines,
   * headers, trailing spaces, etc)
   *   
   * #### Usage  
   *   
   * The aligner can then be invoked by calling from the root of the project (assuming you already built the
   * application using ``build installApp``):
   * 
   * ```
   * build/install/muset/bin/aligner -csvFile src/test/resources/testdataset.csv
   * ```
   * 
   * #### Output
   * 
   * You can see the results in ``results/latest`` (and ``results/all`` for previous experiments), 
   * where you will find, for each iteration:
   * 
   * - Weights (which can be used as initialization later on, see below)
   * - Alignments, in the sub-directories ``alignment-txt`` (for visual inspection) and ``alignment-fasta`` 
   *   (in the standard fasta format).
   *   
   * #### Options and customization
   *   
   * To get help on other functionalities, use:
   * 
   * ```
   * build/install/muset/bin/aligner -help
   * ```
   * 
   * In particular, the help contains information on the following additional options/functionalities:
   * 
   * - ``-initParams``: to load parameters from a previous run as initialization
   * - ``-rocGridSize``: to create alignment with different trade-offs between precision and recall (by 
   *   default, alignments produced maximize recall).
   * - ``-useLongGaps``, ``-addPairSpecific``, ``-featuresFile``, ``-useLetterPairs`` to control the set of 
   *   features (sufficient statistics) used for the alignment exponential family. More details on 
   *   ``-featuresFile`` below.
   *   
   * #### Using features
   * 
   * An attractive aspect of muset's ``Aligner`` is that it is very simple to build letter-specific features used to share
   * statistical strengths across different substitutions and indel operations. 
   * 
   * To do this, set the options ``-featuresFile`` to the path to a file such as the one below:
   * 
   * ```
   * input,features
   * firstLetter,[ f1v1 f2v1 f3v1 ]
   * secondLetter,[ f1v1 f2v2 f3v1 ]
   * thirdLetter,[ someFeatName someOtherFeatName ]
   * fourthLetter,[ someFeatName yetAnotherOne ]
   * ```
   * 
   * Note the mandatory header, the space after the opening bracket and before the closing bracket. The name of a
   * letter-specific feature is arbitrary but should not contain spaces. 
   * 
   * This file is provided in ``src/test/resources/testfeatures.csv`` so its effect can be seen by looking that 
   * the ``weight.txt`` file in ``results/latest/iteration_0/`` after calling:
   * 
   * ```
   * build/install/muset/bin/aligner -csvFile src/test/resources/testdataset.csv -featuresFile src/test/resources/testfeatures.csv
   * ```
   * 
   * A letter-specific feature are viewed 
   * as a vector. Given a substitution between two letters (x,y), sufficient statistics are created by indicating which 
   * letter-specific features are 
   * changed between x and y. For example, and indicator ``featChange(1/2)`` means that  letter-specific feature number 1 
   * (zero indexed) was changed between two letter specific vectors of length 2. When two 
   * letter-specific vectors have different lengths, all letter-specific features of both are considered to be changed, 
   * and a special extra feature is activated 
   * to indicate that a change of feature dimensionality occurred. 
   * 
   */
  @Tutorial(showSource = false)
  public void alignerTest() {}
  
}
