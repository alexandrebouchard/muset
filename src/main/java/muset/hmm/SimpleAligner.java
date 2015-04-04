package muset.hmm;

import java.util.List;
import java.util.Random;


/**
 * this type of aligner can be contextual in the sense of knowning where in the string
 * an indel/sub occurs, but does not know the previous indel/sub operation
 * 
 * Internally, it does keep track of it to avoid overcounting, and the restriction on 
 * not looking on previous op is just to simplify parameter specification.
 * 
 * Use HetPairHMM for full generality
 */
public final class SimpleAligner
{
  private final double [] insLogWeight, delLogWeight; // indexed by pos in the string
  private final double [][] subLogWeigths;
  
  private final int modLastTop, modLastBot;
  public final String topS, botS;
  
  public SimpleAligner(double [] ins, double [] del, double [][] sub)
  {
    int topL = del.length;
    int botL = ins.length;
    this.topS = briefj.opt.StrUtils.repeat("*",topL); this.botS = briefj.opt.StrUtils.repeat("*",botL);
    if (del.length != sub.length || ins.length != sub[0].length)
      throw new RuntimeException();
    this.insLogWeight = ins;
    this.delLogWeight = del;
    this.subLogWeigths = sub;
    this.modLastTop = topS.length(); // i.e. length is topL + 1
    this.modLastBot = botS.length(); // i.e. length is botL + 1
  }
  
  private HetPairHMM _pairHmm = null;
  private double norm = Double.NaN;
  
  public HetPairHMM getHMM()
  {
    ensureHMM();
    return _pairHmm;
  }
  
  public double logSumProduct()
  {
    if (Double.isNaN(norm))
    {
      norm  = 0;
      for (double x : insLogWeight)
        norm += x;
      for (double x : delLogWeight)
        norm += x;
    }
    return getHMM().logSumProduct() + norm;
  }
  
  public Derivation viterbi()
  {
    return HetPairHMM.removeBoundary(getHMM().viterbi(), BOUND);
  }
  
  public Derivation sample(Random rand)
  {
    return HetPairHMM.removeBoundary(getHMM().sample(rand), BOUND);
  }
  
  public Derivation sample(Random rand, List<Integer> fullDerivation)
  {
    return HetPairHMM.removeBoundary(getHMM().sample(rand, fullDerivation), BOUND);
  }
  
  public double pathLogProbability(Derivation d)
  {
    double num = 0.0;
    for (int i = 0; i < d.getCurrentWord().length(); i++)
      if (d.hasAncestor(i))
        num += subLogWeigths[d.ancestor(i)][i] - delLogWeight[d.ancestor(i)] - insLogWeight[i];
    return num - getHMM().logSumProduct();
  }
  
  private void ensureHMM()
  {
    if (_pairHmm != null) return;
    _pairHmm = new HetPairHMM(topS+BOUND, botS+BOUND, new HmmAdaptor());
  }
  
  public final Character BOUND = '#';

  private final class HmmAdaptor implements HetPairHMMSpecification
  {
    @Override public double logWeight(int prevState, int currentState, int x, int y,
        int deltaX, int deltaY)
    {
      final boolean isSub = deltaX == 1 && deltaY == 1,
                    isIns = deltaX == 0 && deltaY == 1,
                    isDel = deltaX == 1 && deltaY == 0;
        
      if (prevState == 2 && currentState == 1) 
        return Double.NEGATIVE_INFINITY;
      if (currentState == 2 && isDel && x != modLastTop)
        return 0.0;
      if (currentState == 1 && isIns && y != modLastBot)
        return 0.0;
      if (currentState == 0 && isSub)
      {
        if (x == modLastTop && y == modLastBot)
          return 0.0;
        else if (x == modLastTop || y == modLastBot)
          return Double.NEGATIVE_INFINITY;
        else
        {
//          return 0.0;
//          System.out.println("" + (idx++) + "\t" + Math.exp(subLogWeigths[x][y] - delLogWeight[x] - insLogWeight[y]));
          return subLogWeigths[x][y] - delLogWeight[x] - insLogWeight[y];
        }
      }
      return Double.NEGATIVE_INFINITY;
    }
    @Override public int endState() { return 0; }
    @Override public int startState() { return 0; }
    @Override public int nStates() { return 3; }
  }
  int idx = 0;
  
//  public static void main(String [] args)
//  {
////    String top = "aa", bot = "a";
//    
//    double [] insLogWeight = new double[]{0.0001,1,2};
//    MtxUtils.logInPlace(insLogWeight);
//    double [] delLogWeight = new double[]{1,1};
//    MtxUtils.logInPlace(delLogWeight);
//    double [][] subLogWeigths = new double[][]{{1,1,3},{1,4,1}};
//    for (double [] mtx : subLogWeigths)
//      MtxUtils.logInPlace(mtx);
//    
//    SimpleAligner aligner = new SimpleAligner(insLogWeight, delLogWeight, subLogWeigths);
//    
//
//    
//    System.out.println(Math.exp(aligner.logSumProduct()));
//    System.out.println(aligner.viterbi());
//    
//    Random rand = new Random(1);
//    Counter<Derivation> derivs = new Counter<Derivation>();
//    CounterMap<Derivation,List<Integer>> fullDerivations = new CounterMap<Derivation,List<Integer>>();
//    for (int i =0 ; i < 100000; i ++)
//    {
//      List<Integer> fullDeriv = list();
//      Derivation deriv = aligner.sample(rand, fullDeriv);
//      derivs.incrementCount(deriv,1.0);
//      fullDerivations.incrementCount(deriv, fullDeriv, 1.0);
//    }
//    derivs.normalize();
//    double sum = 0.0;
//    System.out.println("-------------------------------------------------------");
//    for (Derivation d : derivs)
//    {
//      System.out.println(d);
//      System.out.println("N full derivs:" + fullDerivations.getCounter(d).size());
//      final double curAnalytic = Math.exp(aligner.pathLogProbability(d));
//      sum += curAnalytic;
//      System.out.println("MC:" + derivs.getCount(d) + "\tAnalytic:" + curAnalytic);
//      System.out.println("-------------------------------------------------------");
//    }
//    System.out.println("\t\tSum:" + sum);
//  }

}
