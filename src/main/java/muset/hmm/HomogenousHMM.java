package muset.hmm;

import muset.Alphabet;
import muset.Alphabet.Letter;
import muset.Sequence;

/**
 * Homogenous is a special case of Heterogenous
 * Note: for convenience, the costs are not in log space here,
 * they are converted in createPairHMM()
 * 
 * All the weights are initialized to zero 
 * @author bouchard
 */
public final class HomogenousHMM
{
  private final Alphabet alphabet;
  private final double [][][][] sub;   // prev state -> cur state -> top char -> bot char
  private final double [][][]   del;   // prev state -> cur state -> top char
  private final double [][][]   ins;   // prev state -> cur state -> bot char
  private final int startState, endState;
  private int c2i(Letter c) { return alphabet.indexer.o2i(c); }
  public HomogenousHMM(Alphabet alphabet, int nStates, int startState, int endState)
  {
    this.alphabet = alphabet;
    this.startState = startState;
    this.endState = endState;
    this.sub = new double[nStates][nStates][alphabet.indexer.size()][alphabet.indexer.size()];
    this.del = new double[nStates][nStates][alphabet.indexer.size()];
    this.ins = new double[nStates][nStates][alphabet.indexer.size()];
  }
  public int nStates() { return sub.length; }
  public void setSub(int s1, int s2, Letter top, Letter bot, double value)
  {
    this.sub[s1][s2][c2i(top)][c2i(bot)] = value;
  }
  public void setIns(int s1, int s2, Letter bot, double value)
  {
    this.ins[s1][s2][c2i(bot)] = value;
  }
  public void setDel(int s1, int s2, Letter top, double value)
  {
    this.del[s1][s2][c2i(top)] = value;
  }
  public HetPairHMM createPairHMM(Sequence top, Sequence bot)
  {
    final double [][][][] logsub = new double[nStates()][nStates()][top.length()][bot.length()];   // prev state -> cur state -> top idx -> bot idx
    final double [][][]   logdel = new double[nStates()][nStates()][top.length()];   // prev state -> cur state -> top idx
    final double [][][]   logins = new double[nStates()][nStates()][bot.length()];
    for (int s1 = 0; s1 < nStates(); s1++)
      for (int s2 = 0; s2 < nStates(); s2++)
      {
        for (int t = 0; t < top.length(); t++)
          for (int b = 0; b < bot.length(); b++)
            logsub[s1][s2][t][b] = Math.log(sub[s1][s2][c2i(top.letterAt(t))][c2i(bot.letterAt(b))]);
        for (int t = 0; t < top.length(); t++)
          logdel[s1][s2][t] = Math.log(del[s1][s2][c2i(top.letterAt(t))]);
        for (int b = 0; b < bot.length(); b++)
          logins[s1][s2][b] = Math.log(ins[s1][s2][c2i(bot.letterAt(b))]);
      }
    final HetPairHMMSpecification specs = new HetPairHMMSpecification() {
      @Override public int nStates() { return logsub.length; }
      @Override
      public double logWeight(int prevState, int currentState, int x, int y,
          int deltaX, int deltaY)
      {
        if (deltaX == 1 && deltaY == 1)
         return logsub[prevState][currentState][x][y];
        else if (deltaX == 1)
          return logdel[prevState][currentState][x];
        else if (deltaY == 1)
          return logins[prevState][currentState][y];
        else
          throw new RuntimeException();
      }
      @Override public int endState() { return endState; }
      @Override public int startState() { return startState; }
    };
    return new HetPairHMM(top, bot, specs);
  }
}