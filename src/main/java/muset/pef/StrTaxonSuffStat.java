package muset.pef;
import briefj.Indexer;
import muset.Sequence;
import muset.SequenceId;

public abstract class StrTaxonSuffStat
{
  public abstract StrTaxonSuffStatExtractor getExtractor(
      Sequence str1, 
      Sequence str2, 
      SequenceId taxon1, 
      SequenceId taxon2);
  
  @SuppressWarnings("rawtypes")
  public Indexer valuesIndexer = new Indexer(); // all
  
  public static interface StrTaxonSuffStatExtractor
  {
    public int extract(int position1, int position2);
  }
}
