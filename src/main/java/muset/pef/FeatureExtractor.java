package muset.pef;
import java.util.Set;

import muset.SequenceId;
import muset.pef.StrTaxonSuffStat.StrTaxonSuffStatExtractor;
import bayonet.regression.LabeledInstance;
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
  
  public static class FeatureOptions
  {
    @Option public boolean useLongGaps = true;
    @Option public boolean addPairSpecific = true;
  }
  
  @Override
  public Counter<Object> extractFeatures(LabeledInstance<Input, Output> instance)
  {
    final Counter<Object> result = new Counter<Object>();
    
//    final int [] coord = cp.int2coord(instance.getInput().strTaxSuffStat);
//    final int distQuantile = coord[Index.DIST.ordinal()];
//    final int hIndex = coord[Index.HC.ordinal()];
    for (String prefix : (options.addPairSpecific ? new String[]{"", sequenceIdPairsIndex.i2o(instance.getInput().strTaxSuffStat).toString() + ","} : new String[]{""}))
    {
    
      // pairs
      final String emission = prefix + (instance.getLabel().topSymbol < instance.getLabel().botSymbol ?
          "(" +  instance.getLabel().topToChar() + instance.getLabel().botToChar() + ")" :
          "(" +  instance.getLabel().botToChar() + instance.getLabel().topToChar() + ")" ) ;
      result.incrementCount(emission, 1.0);
      
      // self-sub
      if (instance.getLabel().topSymbol == instance.getLabel().botSymbol) // self sub
        result.incrementCount(prefix + "selfsub", 1.0);
    
      
      // is it del, sub or ins?
      result.incrementCount(prefix + "state1=" + collapseInDelStates(instance.getInput().state1), 1.0); 
      
      // subsumes long ins, del
      if (options.useLongGaps)
        result.incrementCount(prefix + "," + "state1=" + collapseInDelStates(instance.getInput().state1) + 
            "&state2=" + collapseInDelStates(instance.getLabel().state2), 1.0); 
    }
    
    return result;
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
//      for (int i = 0; i < cp.max; i++)
//      {
//        StringBuilder current = new StringBuilder();
//        int [] c = cp.int2coord(i);
//        for (int j = 0; j < indices.length; j++)
//        {
//          current.append(indices[j].toString() + "=" + c[j]);
//          if (j != indices.length -1)
//            current.append(",");
//        }
//        valuesIndexer.addToIndex(current.toString()); 
//      }
    } 
    @Override
    public StrTaxonSuffStatExtractor getExtractor(
        String str1, 
        String str2,
        SequenceId taxon1, 
        SequenceId taxon2)
    {
      return new Extractor(options.addPairSpecific ?
          sequenceIdPairsIndex.o2i(UnorderedPair.of(taxon1, taxon2)) :
          0);
//      int [] c = new int[indices.length];
//      // find largest threshold that has current distance >= to it
//      if (options.nQuantiles > 1)
//      {
//        final double d = distances.get(new UnorderedPair<Taxon,Taxon>(taxon1,taxon2));
//        loopThresholds : for (int i = thresholds.length - 1; i >= 0; i--)
//          if (d >= thresholds[i])
//          {
//            c[Index.DIST.ordinal()] = i;
//            break loopThresholds;
//          }
//      }
//      return new Extractor(c, 
//          options.hydrophobicModeling ? hydrophobicModelingHeuristic.compute(str1) : null, 
//          options.hydrophobicModeling ? hydrophobicModelingHeuristic.compute(str2) : null);
    }
  }
  
  public static class Extractor implements StrTaxonSuffStatExtractor
  {
    private final int pairId;
    public Extractor(int pairId) 
    { 
      this.pairId = pairId;
//      this.c = c; 
//      this.topHydro = topHydro;
//      this.botHydro = botHydro;
    }
//    private final int [] c;
//    private final int [] topHydro, botHydro;
    @Override
    public int extract(int position1, int position2)
    {
      return pairId;
//      if (topHydro != null && botHydro != null)
//      {
//        if (position1 < topHydro.length && position2 < botHydro.length)
//          c[Index.HC.ordinal()] = topHydro[position1] + botHydro[position2];
//        else
//          c[Index.HC.ordinal()] = 0;
//      }
//      return cp.coord2int(c);
    }
  }
  
  public FeatureExtractor(
      Set<UnorderedPair<SequenceId, SequenceId>> sequenceIdPairs,  FeatureOptions options)
  {
    this.options = options;
    this.sequenceIdPairsIndex = new Indexer<UnorderedPair<SequenceId,SequenceId>>(sequenceIdPairs);
//    this.distances = distances;
//    // find distance quantiles
//    thresholds = new double[options.nQuantiles];
//    if (options.nQuantiles > 1)
//    {
//      DescriptiveStatistics ds = new DescriptiveStatistics();
//      for (UnorderedPair<Taxon,Taxon> key : distances.keySet())
//        if (!key.getFirst().equals(key.getSecond()))
//          ds.addValue(distances.get(key));
//      // first coordinate is a hack
//      double [] quantiles = new double[options.nQuantiles]; //{0.5/3.0*100, 1.0/3.0*100, 2.0/3.0*100};
//      for (int i = 1; i < quantiles.length; i++)
//        quantiles[i] = 100.0 * ((double) i)/((double) quantiles.length);
//      quantiles[0] = Double.NaN;
//      
//      for (int i = 0; i < quantiles.length; i++)
//        thresholds[i] = (i == 0 ? Double.NEGATIVE_INFINITY : ds.getPercentile(quantiles[i]));
//      LogInfo.logsForce("Thresholds:" + Arrays.toString(thresholds));
//    }
//    //
//    indices = Index.values();
//    final int [] sizes = new int[indices.length];
//    for (Index idx : indices)
//      sizes[idx.ordinal()] = idx.size(FeatureExtractor.this);
//    cp = new MSCoordinatePacker(sizes);
    //
//    hydroDB = (options.hydrophobicModeling ? new HydropathyDB() : null);
  }
}
