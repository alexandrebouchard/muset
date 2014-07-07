package muset;

import java.io.File;
import java.util.Map;

import muset.util.Edge;

import org.junit.Test;

import briefj.BriefFiles;
import briefj.collections.Counter;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import tutorialj.Tutorial;

import static org.junit.Assert.*;

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
   *   [http://www.eecs.berkeley.edu/Pubs/TechRpts/2010/EECS-2010-153.pdf]. 
   * 3. Variational inference algorithms for MSAs 
   *   [http://papers.nips.cc/paper/4036-variational-inference-over-combinatorial-spaces.pdf]
   *  
   * Note: as of Jul 7 2014, 1. is operational, while code for 2 and 3 is being transferred from a legacy SVN library.
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
   */
  @Tutorial(startTutorial = "README.md", showSource = false)
  public void installInstructions()
  {
  }
  
  /**
   * Using MSAPoset
   * --------------
   */
  @Tutorial(showLink = true, linkPrefix = "src/test/java/")
  @Test
  public void msaPosetTest()
  {
    // creating an alignment
    Map<SequenceId,String> sequences = Maps.newLinkedHashMap();
    SequenceId 
      seq0 = new SequenceId("seq-0"),
      seq1 = new SequenceId("seq-1"),
      seq2 = new SequenceId("seq-2");
    sequences.put(seq0, "ACA");
    sequences.put(seq1, "CC");
    sequences.put(seq2, "ACC");
    MSAPoset align = new MSAPoset(sequences);
    
    // add links
    // NB: tryAdding returns whether adding the link is allowed by MSA constraints
    Edge 
      e0 = new Edge(1, 0, seq0, seq1),
      e1 = new Edge(0, 1, seq1, seq2);
    assertTrue(align.tryAdding(e0));
    assertTrue(align.tryAdding(e1));
    
    System.out.println(align);
    /*
      A-C-A-|seq-0
      --CC--|seq-1
      -AC--C|seq-2
     */
    
    // check for the presence of a link
    // note: transitive links are taken care automatically
    Edge
      e2 = new Edge(1, 1, seq0, seq2);
    assertTrue(align.containsEdge(e2));
    
    // edges are not added if they violate MSA constraints
    Edge
      e3 = new Edge(0, 2, seq0, seq2);
    assertFalse(align.tryAdding(e3));
    
    // deep clone
    MSAPoset clone = new MSAPoset(align);
    Edge 
      e4 = new Edge(1, 2, seq1, seq2);
    assertTrue(align.tryAdding(e4));
    assertTrue(align.nEdges() == 4);
    assertTrue(clone.nEdges() == 3);
    
    // create an alignment from edge scores
    Counter<Edge> scores = new Counter<Edge>();
    scores.setCount(e0, 100.0);
    scores.setCount(e1, -5.0);
    scores.setCount(e2, 10.5);
    scores.setCount(e3, 1000.0);
    scores.setCount(e4, -10.0);
    MSAPoset msa = MSAPoset.maxRecallMSA(sequences, scores);
    assertTrue(msa.nEdges() == 2);
    
    System.out.println(msa);
    /*
      --AC-A|seq-0
      ---CC-|seq-1
      ACC---|seq-2
     */
    
    // compute alignment metrics
    assertEquals(Metrics.edgeF1(msa, align), 1.0/3.0, 0.0);
    
    // remove links
    align.split(align.column(seq0, 1), Sets.newHashSet(seq0));
    assertTrue(align.nEdges() == 2);
    
    System.out.println(align);
    /*
      A-C--A|seq-0
      ---CC-|seq-1
      -A-CC-|seq-2
     */
    
    // input and output to and from FASTA format
    File tempFile = BriefFiles.createTempFile();
    align.toFASTA(tempFile);
    MSAPoset readMSA = MSAPoset.parseFASTA(tempFile);
    assertTrue(MSAPoset.deepEquals(readMSA, align));
  }
}
