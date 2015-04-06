package muset.hmm;

import java.util.List;
import java.util.Random;

import muset.Alphabet.Letter;
import muset.Sequence;
import bayonet.distributions.Multinomial;
import bayonet.math.CoordinatePacker;
import bayonet.math.NumericalUtils;



/**
 * Implements Heterogeneous pair HMMs, meaning that the index in the top and bottom
 * strings can be used to determine the cost of the current transition
 * 
 * A Mealey automata
 * @author bouchard
 *
 */
public final class HetPairHMM
{
  private int startState() { return hmm.startState(); }
  private int endState() { return hmm.endState(); }
  
  private double [][][] prefix, suffix, maxSuffix; // state -> x pos -> y pos
  public final HetPairHMMSpecification hmm;
  private final int nStates;
  public final Sequence str1, str2;
  private boolean fwdInitialized = false, bwdInitialized = false, bwdMaxInitialized = false;
  
  public HetPairHMM(Sequence str1, Sequence str2, HetPairHMMSpecification pairHMM)
  {
    this.str1 = str1;
    this.str2 = str2;
    this.hmm = pairHMM;
    this.nStates = pairHMM.nStates();
  }
  
  public double logSumProduct()
  {
    return prefixLogSumProduct(endState(), str1.length(), str2.length()); 
  }
  
  // warning: not subtracting logSumProduct (i.e. unnormalized!)
  public double logSumProduct(int state1, int state2, int x, int y, int deltaX, int deltaY)
  {
    return prefixLogSumProduct(state1, x, y) + hmm.logWeight(state1, state2, x, y, deltaX, deltaY) + 
           suffixLogSumProduct(state2, str1.length()-x-deltaX, str2.length()-y-deltaY);
  }
  
  /**
   * The log posterior pr that these two points are aligned
   * @param x
   * @param y
   * @return
   */
  public double logPosteriorAlignment(int x, int y)
  {
    double logSum = Double.NEGATIVE_INFINITY;
    for (int s1 = 0; s1 < nStates; s1++)
      for (int s2 = 0; s2 < nStates; s2++)
        logSum = NumericalUtils.logAdd(logSum, logSumProduct(s1, s2, x, y, 1, 1));
    return logSum - logSumProduct();
  }
  
  /**
   * Forward recursion
   * The log of the sum of the weight of all the paths starting at the startState() and ending in state finalState,
   * and emitting the x first symbols of str1, and the y first symbols of str2
   */
  public double prefixLogSumProduct(int finalState, int x, int y)
  {
    if (!fwdInitialized) computeForward();
    return prefix[finalState][x][y];
  }
  
  private void computeForward()
  {
    final int 
      len1 = str1.length(),
      len2 = str2.length();
    if (startState() != 0 || endState() != 0) throw new RuntimeException();
    this.prefix = new double[hmm.nStates()][len1+1][len2+1];
    final double [][][] prefix = this.prefix;
    for (int x = 0; x <= len1; x++)
      for (int y = 0; y <= len2; y++)
        for (int finalState = 0; finalState < nStates; finalState++)
        {
          double result = Double.NEGATIVE_INFINITY;
          if (x == 0 && y == 0) 
          {
            if (finalState == startState()) result = 0.0; // log(1.0)
            else                            result = Double.NEGATIVE_INFINITY; // log(0.0)
          }
          else
          {
            if (x > 0)
              for (int previousState = 0; previousState < nStates; previousState++)
                result = NumericalUtils.logAdd(result, prefix[previousState][x-1][y]   + hmm.logWeight(previousState, finalState, x-1,y,1,0));
            if (y > 0)
              for (int previousState = 0; previousState < nStates; previousState++)
                result = NumericalUtils.logAdd(result, prefix[previousState][x][y-1]   + hmm.logWeight(previousState, finalState, x,y-1,0,1));
            if (x > 0 && y > 0)
              for (int previousState = 0; previousState < nStates; previousState++)
                result = NumericalUtils.logAdd(result, prefix[previousState][x-1][y-1] + hmm.logWeight(previousState, finalState, x-1,y-1,1,1));
          }
          prefix[finalState][x][y] = result;
        }
    fwdInitialized = true;
  }

  /**
   * Backward recursion:
   * The log of the sum of the weight of all the paths starting in state firstState ending in endState()
   * and emitting the x last symbols of str1, i.e [str1.length-x, str1.length()) and the y last symbols of str1
   * Note: actually computed in order of increasing x, y but represent suffix being grown from right to left 
   */
  public double suffixLogSumProduct(int firstState, int x, int y)
  {
    if (!bwdInitialized) computeBackward();
    return suffix[firstState][x][y];
  }
  
  private void computeBackward()
  {
    final int 
      len1 = str1.length(),
      len2 = str2.length();
    this.suffix = new double[hmm.nStates()][len1+1][len2+1];
    final double [][][] suffix = this.suffix;
    for (int x = 0; x <= len1; x++)
      for (int y = 0; y <= len2; y++)
        for (int firstState = 0; firstState <nStates; firstState++)
        {
          double result = Double.NEGATIVE_INFINITY;
          if (x == 0 && y == 0) 
          {
            if (firstState == endState())     result = 0.0; // log(1.0)
            else                              result = Double.NEGATIVE_INFINITY; // log(0.0)
          }
          else
          {
            if (x > 0)
              for (int nextState = 0; nextState < nStates; nextState++)
                result = NumericalUtils.logAdd(result, suffix[nextState][x-1][y]   + hmm.logWeight(firstState, nextState, len1-x,len2-y,1,0));
            if (y > 0)
              for (int nextState = 0; nextState < nStates; nextState++)
                result = NumericalUtils.logAdd(result, suffix[nextState][x][y-1]   + hmm.logWeight(firstState, nextState, len1-x,len2-y,0,1));
            if (x > 0 && y > 0)
              for (int nextState = 0; nextState < nStates; nextState++)
                result = NumericalUtils.logAdd(result, suffix[nextState][x-1][y-1] + hmm.logWeight(firstState, nextState, len1-x,len2-y,1,1));
          }
          suffix[firstState][x][y] = result;
        }
    bwdInitialized = true;
  }
  
  public double suffixLogMaxProduct(int firstState, int x, int y)
  {
    if (!bwdMaxInitialized) computeMaxBackward();
    return maxSuffix[firstState][x][y];
  }
  
  private void computeMaxBackward()
  {
    this.maxSuffix = new double[hmm.nStates()][str1.length()+1][str2.length()+1];
    final double [][][] maxSuffix = this.maxSuffix;
    final int 
      len1 = str1.length(),
      len2 = str2.length();
    for (int x = 0; x <= len1; x++)
      for (int y = 0; y <= len2; y++)
        for (int firstState = 0; firstState <nStates; firstState++)
        {
          double result = Double.NEGATIVE_INFINITY;
          if (x == 0 && y == 0) 
          {
            if (firstState == endState())   result = 0.0; // log(1.0)
            else                            result = Double.NEGATIVE_INFINITY; // log(0.0)
          }
          else
          {
            if (x > 0)
              for (int nextState = 0; nextState < nStates; nextState++)
                result = Math.max(result, maxSuffix[nextState][x-1][y]   + hmm.logWeight(firstState, nextState, len1-x,len2-y,1,0));
            if (y > 0)
              for (int nextState = 0; nextState < nStates; nextState++)
                result = Math.max(result, maxSuffix[nextState][x][y-1]   + hmm.logWeight(firstState, nextState, len1-x,len2-y,0,1));
            if (x > 0 && y > 0)
              for (int nextState = 0; nextState < nStates; nextState++)
                result = Math.max(result, maxSuffix[nextState][x-1][y-1] + hmm.logWeight(firstState, nextState, len1-x,len2-y,1,1));
          }
          maxSuffix[firstState][x][y] = result;
        }
    bwdMaxInitialized = true;
  }
  
  public static Derivation removeBoundary(Derivation d, Letter bound)
  {
    final int 
      oldLastBotIndex = d.getCurrentWord().length() - 1,
      oldLastTopIndex = d.getAncestorWord().length() - 1;
    if (!d.getCurrentWord().letterAt(oldLastBotIndex).equals(bound) ||
        !d.getAncestorWord().letterAt(oldLastTopIndex).equals(bound) ||
        d.ancestor(oldLastBotIndex) != oldLastTopIndex)
      throw new RuntimeException();
    int [] newAnc = new int[oldLastBotIndex];
    for (int i = 0; i < newAnc.length ; i++)
      if (d.hasAncestor(i))
        newAnc[i] = d.ancestor(i);
      else
        newAnc[i] = Derivation.INSERTED;
    return new Derivation(newAnc, 
        d.getAncestorWord().subsequence(0, oldLastTopIndex),
        d.getCurrentWord(). subsequence(0, oldLastBotIndex));
  }
  
  public Derivation viterbi(List<Integer> stateSequence)
  {
    return _viterbiOrSample(stateSequence, false, null);
  }
  
  public Derivation viterbi()
  {
    return _viterbiOrSample(null, false, null);
  }
  
  public Derivation sample(Random rand)
  {
    return _viterbiOrSample(null, true, rand);
  }
  
  public Derivation sample(Random rand, List<Integer> stateSequence)
  {
    return _viterbiOrSample(stateSequence, true, rand);
  }
  
  private Derivation _viterbiOrSample(List<Integer> stateSequence, boolean sample, Random rand)
  {
    if (stateSequence != null && stateSequence.size() > 0) throw new RuntimeException();
    final int [] ancestors = new int[str2.length()];
    for (int i = 0; i < ancestors.length; i++)
      ancestors[i] = Derivation.INSERTED - 1;
    int x = str1.length(), y = str2.length(); // number of characters left in each sequence (i.e. not processed yet)
    int previousState = startState();
    if (stateSequence != null) stateSequence.add(startState());
    CoordinatePacker cp = getPacker();
    double [] choices = (sample ? new double[2*2*nStates] : null);
    while (x > 0 || y > 0)
    {
      // need to pick a next state and which characters to consume
      int argmax = -1;
      double max = Double.NEGATIVE_INFINITY;
      if (sample)
        for (int i =0 ; i < choices.length; i++)
          choices[i] = Double.NEGATIVE_INFINITY;
      for (int currentState = 0; currentState < nStates; currentState++)
      {
        if (x > 0)
        {
          final double current = (sample ? suffixLogSumProduct(currentState,x-1,y) : suffixLogMaxProduct(currentState,x-1,y))   + hmm.logWeight(previousState, currentState, str1.length()-x,str2.length()-y,1,0);
          final int currentCoord = cp.coord2int(1,0,currentState);
          if (sample)
            choices[currentCoord] = current;
          if (current > max)
          {
            max = current;
            argmax = currentCoord;
          }
        }
        if (y > 0)
        {
          final double current = (sample ? suffixLogSumProduct(currentState,x,y-1) : suffixLogMaxProduct(currentState,x,y-1))   + hmm.logWeight(previousState, currentState, str1.length()-x,str2.length()-y,0,1);
          final int currentCoord = cp.coord2int(0,1,currentState);
          if (sample)
            choices[currentCoord] = current;
          if (current > max)
          {
            max = current;
            argmax = currentCoord;
          }
        }
        if (x > 0 && y > 0)
        {
          final double current = (sample ? suffixLogSumProduct(currentState, x-1,y-1) : suffixLogMaxProduct(currentState, x-1,y-1)) + hmm.logWeight(previousState, currentState, str1.length()-x,str2.length()-y,1,1);
          final int currentCoord = cp.coord2int(1,1,currentState);
          if (sample)
            choices[currentCoord] = current;
          if (current > max)
          {
            max = current;
            argmax = currentCoord;
          }
        }
      }
      //
      if (Double.isInfinite(max))
        throw new RuntimeException();
      if (sample)
      {
        // redefine argmax to be taken from the sampled distribution
        Multinomial.expNormalize(choices);
        argmax = Multinomial.sampleMultinomial(rand, choices);
      }
      int [] coords = cp.int2coord(argmax);
      if (stateSequence != null) stateSequence.add(coords[2]);
      final int positionInTop = str1.length() - x;
      final int positionInBot = str2.length() - y;
      if (coords[0] == 1 && coords[1] == 1)
        ancestors[positionInBot] = positionInTop;
      else if (coords[0] == 0 && coords[1] == 1)
        ancestors[positionInBot] = Derivation.INSERTED;
      else
        ;
      x = x - coords[0];
      y = y - coords[1];
      previousState = coords[2];
    }
    return new Derivation(ancestors, str1, str2);
  }
  
  private CoordinatePacker getPacker()
  {
    int [] sizes = new int[3];
    sizes[0] = sizes[1] = 2;
    sizes[2] = nStates;
    return new CoordinatePacker(sizes);
  }

}
