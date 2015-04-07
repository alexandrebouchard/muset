package muset;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.Map;

import muset.util.Edge;

import org.junit.Test;

import tutorialj.Tutorial;
import briefj.BriefFiles;
import briefj.collections.Counter;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;



public class MSAPosetTest
{
  @Tutorial(showLink = true, linkPrefix = "src/test/java/")
  @Test
  public void msaPosetTest()
  {
    // creating an alignment
    Map<SequenceId,Sequence> sequences = Maps.newLinkedHashMap();
    Alphabet alphabet = new Alphabet();
    SequenceId 
      seq0 = new SequenceId("seq-0"),
      seq1 = new SequenceId("seq-1"),
      seq2 = new SequenceId("seq-2");
    sequences.put(seq0, Sequence.buildSimpleSequence(alphabet, "ACA"));
    sequences.put(seq1,  Sequence.buildSimpleSequence(alphabet, "CC"));
    sequences.put(seq2,  Sequence.buildSimpleSequence(alphabet, "ACC"));
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
    Alphabet alpha = new Alphabet();
    MSAPoset readMSA = MSAPoset.parseFASTA(alpha, tempFile);
    assertTrue(MSAPoset.deepEquals(readMSA, align));
  }
}
