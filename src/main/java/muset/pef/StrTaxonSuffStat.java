package muset.pef;
import briefj.Indexer;
import muset.SequenceId;

public abstract class StrTaxonSuffStat
{
  public abstract StrTaxonSuffStatExtractor getExtractor(
      String str1, 
      String str2, 
      SequenceId taxon1, 
      SequenceId taxon2);
  
  @SuppressWarnings("rawtypes")
  public Indexer valuesIndexer = new Indexer(); // all
  
  public static interface StrTaxonSuffStatExtractor
  {
    public int extract(int position1, int position2);
  }
}
