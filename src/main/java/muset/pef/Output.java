package muset.pef;

import muset.Alphabet.Letter;

public final class  Output  implements Comparable<Output>
{
  public final int state2, topSymbol, botSymbol;
  private final Model model;
  public Output(int state2, int topSymbol, int botSymbol, Model model)
  {
    this.state2 = state2;
    this.topSymbol = topSymbol;
    this.botSymbol = botSymbol;
    this.model = model;
    if (topSymbol == model.epsilon() && botSymbol == model.epsilon()) // we do not allow pure epsilon for now
      throw new RuntimeException();
  }
  @Override public int compareTo(Output other)
  {
    if (this.state2 < other.state2) return -1;
    if (this.state2 == other.state2) 
    {
      if (this.topSymbol < other.topSymbol) return -1;
      if (this.topSymbol == other.topSymbol) 
      {
        if (this.botSymbol < other.botSymbol) return -1;
        if (this.botSymbol == other.botSymbol) return 0;
        else return 1;
      }
      else return 1;
    }
    else return 1;
  }
  @Override
  public int hashCode()
  {
    final int prime = 31;
    int result = 1;
    result = prime * result + botSymbol;
    result = prime * result + ((model == null) ? 0 : model.hashCode());
    result = prime * result + state2;
    result = prime * result + topSymbol;
    return result;
  }
  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    Output other = (Output) obj;
    if (botSymbol != other.botSymbol)
      return false;
    if (model == null)
    {
      if (other.model != null)
        return false;
    } else if (!model.equals(other.model))
      return false;
    if (state2 != other.state2)
      return false;
    if (topSymbol != other.topSymbol)
      return false;
    return true;
  }
  @Override
  public String toString()
  {
    return "[state2=" + model.stateIndexer.i2o(state2) + "," +
    		"topSym=" + topToChar() + "," +
    		"botSym=" + botToChar() + "]";
  }
  public Letter toLetter(int symbolId)
  {
    if (symbolId != model.epsilon())
      return model.enc.indexer.i2o(symbolId);
    return null;
  }
  public Letter topToChar() { return toLetter(topSymbol); }
  public Letter botToChar() { return toLetter(botSymbol);   }
}
