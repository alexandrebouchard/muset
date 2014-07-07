package muset;

import java.util.Map;
import java.util.Set;

import muset.util.Edge;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import com.google.common.collect.Maps;



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
    for (String seq : msa.sequences().values())
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

  /**
   * keeps only capitalized links, then capitalize everything
   * @param msa
   * @return
   */
  public static MSAPoset processBenchmarkReference(MSAPoset msa)
  {  
    MSAPoset result = new MSAPoset(msa.sequences);
    result.disableLinearization();
    for (MSAPoset.Column c : msa.linearizedColumns.keySet())
    {
      Map<SequenceId,Integer> processed = Maps.newLinkedHashMap();
      for (SequenceId item : c.points.keySet())
        if (Character.isUpperCase(msa.charAt(c, item)))
          processed.put(item, c.points.get(item));
        result.tryAdding(processed);
    }
    for (SequenceId sequenceId : result.allSequenceIds())
      result.sequences.put(sequenceId, result.sequences.get(sequenceId).toUpperCase());
    result.enableLinearization();
    return result;
  }

}
