package muset.hmm;

import java.util.List;

public class SPMinBayesRiskDecoder
{
  public static HetPairHMM getDecoder(final HetPairHMM originalHMM)
  {
    HetPairHMMSpecification specs = new HetPairHMMSpecification() {
      @Override public int endState() { return originalHMM.hmm.endState(); }
      @Override public double logWeight(int prevState, int currentState, int x, int y,
          int deltaX, int deltaY)
      {
        final double logPostAligned = originalHMM.logSumProduct(prevState, currentState, x, y, deltaX, deltaY);
        if (logPostAligned == Double.NEGATIVE_INFINITY)
          return Double.NEGATIVE_INFINITY;
        final boolean isAlignmentPoint = (deltaX == 1 && deltaY == 1);
        if (isAlignmentPoint)
          return Math.exp(logPostAligned);// - Math.log(1.0 - Math.min(1.0,Math.exp(logPostAligned)));
        else
          return 0;
      }
      @Override public int nStates() { return originalHMM.hmm.nStates(); }
      @Override public int startState() { return originalHMM.hmm.startState(); }
    };
    return new HetPairHMM(originalHMM.str1, originalHMM.str2, specs);
  }
  
  public static Derivation decode(final HetPairHMM originalHMM)
  {
    final HetPairHMM decoder = getDecoder(originalHMM);
    List<Integer> list = null; //new ArrayList<Integer>();
    final Derivation result = decoder.viterbi(list);
    return result;
  }
}
