package muset.pef;

public final class Input
{
  public final int state1, strTaxSuffStat;
  public final Model model;
  public Input(int state1, int strTaxSuffStat, Model model)
  {
    this.state1 = state1;
    this.strTaxSuffStat = strTaxSuffStat;
    this.model = model;
  }

  @Override
  public String toString()
  {
    return "[state1=" + model.stateIndexer.i2o(state1) + "," +
    		"topSuffStat="+model.stSuffStat.valuesIndexer.i2o(strTaxSuffStat)   +"]";
  }

  @Override
  public int hashCode()
  {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((model == null) ? 0 : model.hashCode());
    result = prime * result + state1;
    result = prime * result + strTaxSuffStat;
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
    Input other = (Input) obj;
    if (model == null)
    {
      if (other.model != null)
        return false;
    } else if (!model.equals(other.model))
      return false;
    if (state1 != other.state1)
      return false;
    if (strTaxSuffStat != other.strTaxSuffStat)
      return false;
    return true;
  }
  
}
