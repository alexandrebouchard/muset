package muset.pef;
import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Splitter;

import muset.Alphabet;
import muset.Alphabet.Letter;
import muset.Sequence;
import muset.SequenceId;
import muset.pef.StrTaxonSuffStat.StrTaxonSuffStatExtractor;
import bayonet.regression.LabeledInstance;
import briefj.BriefIO;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;



public final class FeatureExtractor implements
    bayonet.regression.FeatureExtractor<LabeledInstance<Input, Output>,Object>
{
  private static final long serialVersionUID = 1L;
  public final FeatureOptions options;
  public final Indexer<UnorderedPair<SequenceId, SequenceId>> sequenceIdPairsIndex;
  private final Alphabet alphabet;
  
  public static class FeatureOptions
  {
    @Option(gloss = "Use a feature that gives a different penalty to indels of one vs. more than one letters")
    public boolean useLongGaps = true;
    
    @Option(gloss = "Create one set of feature that is specific to the pair of taxa being aligned in addition to the default behavior which ignores the identity of the pair of taxa being aligned")
    public boolean addPairSpecific = false;
    
    @Option(gloss = "Path to a feature file. See README.md") 
    public String featuresFile = "";
    
    @Option(gloss = "Use features that are specific to the letters being subsituted/deleted/inserted") 
    public boolean useLetterPairs = true;
  }
  
  // phoneme index -> list of features
  private List<String> [] _features = null;
  
  @SuppressWarnings("unchecked")
  private List<String> [] features()
  {
    if (_features != null)
      return _features;
    
    if (StringUtils.isEmpty(options.featuresFile))
      return null;
    
    _features = new List[alphabet.indexer.size()];
    
    boolean atLeastOneFeatMatch = false;
    loop : for (Map<String,String> line : BriefIO.readLines(new File(options.featuresFile)).indexCSV())
    {
      String phonemeStr = line.get("input");
      if (!alphabet.containsLetter(phonemeStr))
        continue loop;
      Letter phoneme = alphabet.getExistingLetter(phonemeStr);
      String features = line.get("features").replace("[ ", "").replace(" ]", "");
      List<String> parsedFeatures = Splitter.onPattern("\\s+").splitToList(features);
      _features[alphabet.indexer.o2i(phoneme)] = parsedFeatures;
      atLeastOneFeatMatch = true;
    }
    if (!atLeastOneFeatMatch)
      throw new RuntimeException("No feature created. Make sure the alphabets match: " + alphabet.indexer + " vs " + options.featuresFile);
    
    return _features;
  }
  
  @Override
  public Counter<Object> extractFeatures(LabeledInstance<Input, Output> instance)
  {
    List<String> [] features = features();
    final Counter<Object> result = new Counter<Object>();
    
    for (String prefix : (options.addPairSpecific ? new String[]{"", sequenceIdPairsIndex.i2o(instance.getInput().strTaxSuffStat).toString() + ","} : new String[]{""}))
    {
      // pairs
      if (options.useLetterPairs)
      {
        final String emission = prefix + (instance.getLabel().topSymbol < instance.getLabel().botSymbol ?
            "pair(" +  toString(instance.getLabel().topToChar()) + "," + toString(instance.getLabel().botToChar()) + ")" :
            "pair(" +  toString(instance.getLabel().botToChar()) + "," + toString(instance.getLabel().topToChar()) + ")" ) ;
        result.incrementCount(emission, 1.0);
      }
      
      if (features != null && 
          instance.getLabel().topSymbol < features.length && 
          instance.getLabel().botSymbol < features.length)
      {
        List<String> 
          feat1 = features[instance.getLabel().topSymbol],
          feat2 = features[instance.getLabel().botSymbol];
        
        if (feat1 != null && feat2 != null) // check needed because of boundary symbols
        {
          if (feat1.size() != feat2.size())
          {
            result.incrementCount(prefix + "featDimChange(" + Math.abs(feat1.size() - feat2.size()) + ")", 1.0);
            for (int d = 0; d < feat1.size(); d++) result.incrementCount("featChange(" + d + "/" + feat1.size() + ")", 1.0);
            for (int d = 0; d < feat2.size(); d++) result.incrementCount("featChange(" + d + "/" + feat2.size() + ")", 1.0);
          }
          else
          {
            final int size = feat1.size();
            for (int d = 0; d < size; d++)
              if (!feat1.get(d).equals(feat2.get(d)))
                result.incrementCount(prefix + "featChange(" + d + "/" + size + ")", 1.0);
          }
        }
      }
      
      // self-sub
      if (instance.getLabel().topSymbol == instance.getLabel().botSymbol) // self sub
        result.incrementCount(prefix + "selfsub", 1.0);
      
      // is it del, sub or ins?
      result.incrementCount(prefix + "state(" + collapseInDelStates(instance.getInput().state1) + ")", 1.0); 
      
      // subsumes long ins, del
      if (options.useLongGaps)
        result.incrementCount(prefix + "statePair(" + collapseInDelStates(instance.getInput().state1) + 
            "," + collapseInDelStates(instance.getLabel().state2) + ")", 1.0); 
    }
    
    return result;
  }
  
  private String toString(Letter letter)
  {
    if (letter == null)
      return "-";
    return letter.toString();
  }

  public static int collapseInDelStates(int s)
  {
    if (s == 0) return 0;
    else if (s == 1 || s == 2) return 1;
    else throw new RuntimeException();
  }

  @Override 
  public double regularizationFactor(Object feature) { return 1.0; }

  public StrTaxonSuffStat getStrTaxonSuffStat()
  {
    return new SuffStat();
  }
  
  public final class SuffStat extends StrTaxonSuffStat
  {
    @SuppressWarnings("unchecked")
    public SuffStat()
    {
      if (options.addPairSpecific)
        this.valuesIndexer.addAllToIndex(sequenceIdPairsIndex.objectsList());
      else
        this.valuesIndexer.addToIndex("NONE");
    } 
    @Override
    public StrTaxonSuffStatExtractor getExtractor(
        Sequence str1, 
        Sequence str2,
        SequenceId taxon1, 
        SequenceId taxon2)
    {
      return new Extractor(options.addPairSpecific ?
          sequenceIdPairsIndex.o2i(UnorderedPair.of(taxon1, taxon2)) :
          0);
    }
  }
  
  public static class Extractor implements StrTaxonSuffStatExtractor
  {
    private final int pairId;
    public Extractor(int pairId) 
    { 
      this.pairId = pairId;
    }

    @Override
    public int extract(int position1, int position2)
    {
      return pairId;
    }
  }
  
  public FeatureExtractor(
      Alphabet alphabet, Set<UnorderedPair<SequenceId, SequenceId>> sequenceIdPairs,  
      FeatureOptions options)
  {
    this.alphabet = alphabet;
    this.options = options;
    this.sequenceIdPairsIndex = new Indexer<UnorderedPair<SequenceId,SequenceId>>(sequenceIdPairs);
  }
}
