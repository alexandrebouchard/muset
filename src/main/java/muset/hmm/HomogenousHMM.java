package muset.hmm;

import briefj.Indexer;

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
  private final Indexer<Character> enc;
  private final double [][][][] sub;   // prev state -> cur state -> top char -> bot char
  private final double [][][]   del;   // prev state -> cur state -> top char
  private final double [][][]   ins;   // prev state -> cur state -> bot char
  private final int startState, endState;
  private int c2i(char c) { return enc.o2i(c); }
  public HomogenousHMM(Indexer<Character> enc, int nStates, int startState, int endState)
  {
    this.enc = enc;
    this.startState = startState;
    this.endState = endState;
    this.sub = new double[nStates][nStates][enc.size()][enc.size()];
    this.del = new double[nStates][nStates][enc.size()];
    this.ins = new double[nStates][nStates][enc.size()];
  }
  public int nStates() { return sub.length; }
  public void setSub(int s1, int s2, char top, char bot, double value)
  {
    this.sub[s1][s2][c2i(top)][c2i(bot)] = value;
  }
  public void setIns(int s1, int s2, char bot, double value)
  {
    this.ins[s1][s2][c2i(bot)] = value;
  }
  public void setDel(int s1, int s2, char top, double value)
  {
    this.del[s1][s2][c2i(top)] = value;
  }
  public HetPairHMM createPairHMM(String top, String bot)
  {
    final double [][][][] logsub = new double[nStates()][nStates()][top.length()][bot.length()];   // prev state -> cur state -> top idx -> bot idx
    final double [][][]   logdel = new double[nStates()][nStates()][top.length()];   // prev state -> cur state -> top idx
    final double [][][]   logins = new double[nStates()][nStates()][bot.length()];
    for (int s1 = 0; s1 < nStates(); s1++)
      for (int s2 = 0; s2 < nStates(); s2++)
      {
        for (int t = 0; t < top.length(); t++)
          for (int b = 0; b < bot.length(); b++)
            logsub[s1][s2][t][b] = Math.log(sub[s1][s2][c2i(top.charAt(t))][c2i(bot.charAt(b))]);
        for (int t = 0; t < top.length(); t++)
          logdel[s1][s2][t] = Math.log(del[s1][s2][c2i(top.charAt(t))]);
        for (int b = 0; b < bot.length(); b++)
          logins[s1][s2][b] = Math.log(ins[s1][s2][c2i(bot.charAt(b))]);
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