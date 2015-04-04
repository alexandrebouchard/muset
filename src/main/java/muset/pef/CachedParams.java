package muset.pef;
import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;

import bayonet.regression.MaxentClassifier;
import muset.SequenceId;
import muset.hmm.HetPairHMMSpecification;



public final class CachedParams
{
  private final double [][][][][] cachedLogPrs; // state 1 -> state 2 -> f-suff-stat -> (symbol1 U EPSILON) -> (symbol2 U EPSILON)
  private final Model model;
  
  /**
   * Weights coming for example from a message passing algorithm (see MessageComputations for weights[][][] semantics---called r messages there)
   */
  public HetPairHMMSpecification getReweightedHMM(final double [][][] logWeights, final String top, final String bot, final SequenceId topTaxon, final SequenceId botTaxon)
  {
    final StrTaxonSuffStat.StrTaxonSuffStatExtractor extractor = 
      model.stSuffStat.getExtractor(top, bot,topTaxon, botTaxon);
    return new HetPairHMMSpecification() {
      @Override public final int startState() { return model.startState; }
      @Override public final int endState() { return model.endState; }
      @Override public final int nStates() { return model.nStates; }
      @Override public final double logWeight(int prevState, int currentState, int xpos, int ypos,
          int deltaX, int deltaY)
      {
        final int 
          xid = model.charIdAt(top, xpos, deltaX),
          yid = model.charIdAt(bot, ypos, deltaY);
        
        final int stss = extractor.extract(xpos, ypos);
        // reweight:
        final boolean isAligned = deltaX == 1 && deltaY == 1;
        final double logWeight = (isAligned ? logWeights[xpos][ypos][1] - logWeights[xpos][ypos][0] : 0.0);
        return logWeight + getLogPr(prevState, currentState, stss, xid, yid);
      }
    };
  }
  
  public HetPairHMMSpecification getUnsupPairHMM(final String top, final String bot, final SequenceId topTaxon, final SequenceId botTaxon)
  {
    final StrTaxonSuffStat.StrTaxonSuffStatExtractor extractor = 
      model.stSuffStat.getExtractor(top, bot,topTaxon, botTaxon);
    
    return new HetPairHMMSpecification() {
      @Override public final int startState() { return model.startState; }
      @Override public final int endState() { return model.endState; }
      @Override public final int nStates() { return model.nStates; }
      @Override public final double logWeight(int prevState, int currentState, int xpos, int ypos,
          int deltaX, int deltaY)
      {
        final int 
          xid = model.charIdAt(top, xpos, deltaX),
          yid = model.charIdAt(bot, ypos, deltaY);
        final int stss = extractor.extract(xpos, ypos);
        return getLogPr(prevState, currentState, stss, xid, yid);
      }
    };
  }
  
  public CachedParams(Model model)
  {
    this.model = model;
    cachedLogPrs = new double[model.nStates][model.nStates][model.stSuffStat.valuesIndexer.size()][model.enc.size()+1][model.enc.size()+1];
    deepFill(cachedLogPrs, Double.NEGATIVE_INFINITY);
  }
  
  private static void deepFill(Object o, double v)
  {
    if (o instanceof double[])
    {
      double [] a = (double []) o;
      Arrays.fill(a, v);
    }
    else if (o instanceof Object[])
    {
      Object [] a = (Object []) o;
      for (Object cur : a)
        deepFill(cur, v);
    }
    else throw new RuntimeException();
  }
  
  private void set(Input in, Output out, double logpr)
  {
    if (logpr > 0.0001) throw new RuntimeException("Invalid log pr:" + logpr);
    cachedLogPrs[in.state1][out.state2][in.strTaxSuffStat][out.topSymbol][out.botSymbol] = logpr;
  }
  public double getLogPr(Input in, Output out)
  {
    return cachedLogPrs[in.state1][out.state2][in.strTaxSuffStat][out.topSymbol][out.botSymbol];
  }
  public final double getLogPr(int state1, int state2, int strTaxSuffStat, int topSymbol, int botSymbol)
  {
    return cachedLogPrs[state1][state2][strTaxSuffStat][topSymbol][botSymbol];
  }

  // this assumes we are working in prs, not logprs, because leaves untouched (to zero) entries not in support
  public CachedParams(Model model,
      MaxentClassifier<Input, Output, Object> maxentClassifier)
  {
    this(model); // this initializes all entries to -infinity
    // for each input, populate outputs
    for (Input in : model.allInputs())
    {
      SortedSet<Output> outs = maxentClassifier.getLabels(in); 
      double [] logprs = maxentClassifier.logProb(in);
      int i = 0;
      for (Output out : outs)
      {
        set(in, out, logprs[i]);
        i++;
      }
    }
  }
  
  @Override
  public String toString()
  {
    StringBuilder result = new StringBuilder();
    List<Output> allOuts = model.allOutputs();
    for (Input in : model.allInputs())
      for (Output out : allOuts)
      {
        double val = getLogPr(in, out);
        if (!Double.isInfinite(val))
          result.append("" + in + "\t" + out + "\t" + Math.exp(val) + "\n");
      }
    return result.toString();
  }

}
