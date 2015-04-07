package muset;

import java.io.File;
import java.util.Map;
import java.util.Set;

import muset.util.Edge;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;


/**
 * Metrics to compare and summarize MSAPoset instances.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class Metrics
{

  /**
   * 
   * @return The average over pairs of sequences and character of identical linked characters.
   */
  public static double getIdentityStatistic(MSAPoset msa) { return basicStat(msa, true); }

  /**
   * 
   * @return The average over pairs of sequences and character of linked characters
   */
  public static double getAlignedStatistics(MSAPoset msa) { return basicStat(msa, false);}
  
  private static double basicStat(MSAPoset msa, boolean isIdentity)
  {
    SummaryStatistics stat = new SummaryStatistics();
    for (SequenceId t1 : msa.allSequenceIds())
      for (SequenceId t2 : msa.allSequenceIds())
        if (!t1.equals(t2))
          stat.addValue(msa.basicStat(t1,t2,isIdentity));
    return stat.getMean();
  }
  
  public static double getMeanSequenceLength(MSAPoset msa) 
  {
    SummaryStatistics result = new SummaryStatistics();
    for (Sequence seq : msa.sequences().values())
      result.addValue(seq.length());
    return result.getMean();
  }

  /**
   * @param gold
   * @param guess
   * @return
   */
  public static double columnRecall(MSAPoset gold, MSAPoset guess)
  {
    Set<Map<SequenceId,Integer>> goldColumns = gold.points(),
      guessColumns = guess.points();
    double num = 0.0, denom = 0.0;
    for (Map<SequenceId,Integer> goldColumn : goldColumns)
      if (goldColumn.size() > 1) // in bali references, this will have the effect of only counting core blocks
      {
        denom++;
        if (guessColumns.contains(goldColumn))
          num++;
      }
    return num/denom;
  }

  /**
   * 
   * @param gold
   * @param guess
   * @return
   */
  public static double edgeRecall(MSAPoset gold, MSAPoset guess)
  {
    double num = 0.0, denom = 0.0;
    for (Edge e : gold.edges())
    {
      denom++;
      if (guess.containsEdge(e))
        num++;
    }
    return num/denom;
  }

  /**
   * 
   * @param gold
   * @param guess
   * @return
   */
  public static double edgePrecision(MSAPoset gold, MSAPoset guess)
  {
    return edgeRecall(guess, gold);
  }

  /**
   * 
   * @param gold
   * @param guess
   * @return
   */
  public static double edgeF1(MSAPoset gold, MSAPoset guess)
  {
    return Metrics.f1Score(edgePrecision(gold, guess), edgeRecall(gold, guess));
  }

  /**
   * 
   * @param precision
   * @param recall
   * @return
   */
  public static double f1Score(double precision, double recall)
  {
    if (precision + recall == 0.0) return 0.0;
    return 2 * (precision * recall) / (precision + recall);
  }
  
  public static void main(String [] args)
  {
    if (args.length < 1 || args.length > 2)
    {
      System.err.println("One or two arguments: path to alignments in fasta format");
      System.err.println("If two are provided, the first is truth, second is guess");
      System.exit(1);
    }
    Alphabet alpha = new Alphabet();
    MSAPoset 
      first = MSAPoset.parseFASTA(alpha, new File(args[0])),
      second = args.length == 2 ? MSAPoset.parseFASTA(alpha, new File(args[1])) : null;
    printStats(first);
    if (second != null)
      printStats(second);
    printComparison(first, second);
  }

  private static void printComparison(MSAPoset first, MSAPoset second)
  {
    System.out.println("columnRecall=" + columnRecall(first, second));
    System.out.println("edgeRecall=" + edgeRecall(first, second));
    System.out.println("edgePrecision=" + edgePrecision(first, second));
    System.out.println("f1Score=" + edgeF1(first, second));
  }

  private static void printStats(MSAPoset msa)
  {
    System.out.println("identityStat=" + getIdentityStatistic(msa));
    System.out.println("alignedStat=" + getAlignedStatistics(msa));
    System.out.println("meanSeqLen=" + getMeanSequenceLength(msa));
    System.out.println("---");
    System.out.println(msa);
    System.out.println("---");
  }
  
  private Metrics() {}
}
