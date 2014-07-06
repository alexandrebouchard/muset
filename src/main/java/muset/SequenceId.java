package muset;

import java.io.Serializable;



public class SequenceId implements Serializable, Comparable<SequenceId>
{
  private static final long serialVersionUID = 1L;
//  public static final Taxon dummy = new Taxon("DUMMY_LANGUAGE");
  protected final String string;
  public SequenceId(String string) 
  { 
//    if (string == null || string.equals(""))
//      throw new RuntimeException("Name of lang should be nontrivial");
    this.string = string; 
  }
  @Override
  public boolean equals(Object obj)
  {
//    if (!(obj instanceof Language)) throw new RuntimeException();
    if (obj == null) return false;
//    if (!(obj instanceof SequenceId))
//      throw new RuntimeException("Only Taxa can be compared to Taxa");
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

  
//  public static class AnnotatedTaxon extends Taxon
//  {
//    public final Serializable annotations;
//    public AnnotatedTaxon(String name, Serializable annotation)
//    {
//      super(name);
//      this.annotations = annotation;
//    }
//  }
  
//  public static class LanguageUtils
//  {
//    public static Arbre<Taxon> convert(Arbre<String> a)
//    {
//      return a.preOrderMap(new ArbreMap<String,Taxon>() {
//        @Override public Taxon map(Arbre<String> currentDomainNode)
//        {
//          return new Taxon(currentDomainNode.getContents());
//        }
//      });
//    }
//  }
  
}
