package muset;

import java.io.Serializable;



public class SequenceId implements Serializable, Comparable<SequenceId>
{
  private static final long serialVersionUID = 1L;
  protected final String string;
  public SequenceId(String string) 
  { 
    this.string = string; 
  }
  @Override
  public boolean equals(Object obj)
  {
    if (obj == null) return false;
    return this.string.equals(((SequenceId) obj).string);
  }
  @Override
  public int hashCode()
  {
    return string.hashCode();
  }
  @Override
  public String toString()
  {
    return string;
  }
  public int compareTo(SequenceId arg0)
  {
    return this.string.compareTo(arg0.string);
  }
}
