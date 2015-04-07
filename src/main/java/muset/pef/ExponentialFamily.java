package muset.pef;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import muset.Alphabet;
import muset.MSAPoset;
import muset.Sequence;
import muset.SequenceId;
import muset.hmm.HetPairHMM;
import muset.pef.FeatureExtractor.FeatureOptions;
import muset.pef.Model.ThreeStatesBaseMeasure;
import muset.util.Edge;
import bayonet.regression.BaseMeasures;
import bayonet.regression.LabeledInstance;
import bayonet.regression.MaxentClassifier;
import bayonet.regression.MaxentClassifier.MaxentOptions;
import briefj.BriefIO;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.run.Results;


/**
 * Pairwise pseudo-likelihood training of Multiple alignments.
 * 
 * Implements the learning method of Chapter 5 of
 * http://www.eecs.berkeley.edu/Pubs/TechRpts/2010/EECS-2010-153.pdf
 * 
 * In order to make this suitable to unsupervised learning, we
 * use normalized decision where the output is the pair of phoneme 
 * (or epsilon) emitted at each step of the alignment process.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public final class ExponentialFamily
{
  private CachedParams cachedParams;
  private Counter<Object> naturalParams, regularizationCenters;
  
  public Counter<LabeledInstance<Input,Output>> suffStats
    = new Counter<LabeledInstance<Input,Output>>();
  
  public final MaxentOptions<Object> learningOptions;
  public final Model model;
  public final BaseMeasures<Input,Output> bm;
  public final FeatureExtractor featureExtractor;
  
  @Override
  public String toString()
  {
    return cachedParams.toString();
  }
  
  public ExponentialFamily(
      Counter<Object> naturalParams,
      Counter<Object> regularizationCenters,
      MaxentOptions<Object> learningOptions, 
      Model model,
      BaseMeasures<Input, Output> bm, FeatureExtractor featureExtractor)
  {
    this.naturalParams = naturalParams;
    this.regularizationCenters = regularizationCenters;
    this.learningOptions = learningOptions;
    this.model = model;
    this.bm = bm;
    this.featureExtractor = featureExtractor;
    initParameters();
  }
  
  public static class ExponentialFamilyOptions
  {
    @Option(gloss = "Path to an init parameters produced in a file called 'weights.txt' in the results folder. "
        + "Use the string 'ZERO' to use a zero vector instead. Note that if a file is used, you should probably "
        + "make sure that the set of feature options are the same (not checked explicitly).")
    public String initParams = ZERO;
    
    @Option(gloss = "Center where the regularization parameters should be centered (e.g., output of simpler model)")
    public String reguCenterParams = SAME;
    
    public static final String SAME = "SAME";
    public static final String ZERO = "ZERO";
    public static final String INTERNAL = "INTERNAL";
    public Counter<Object> internal = new Counter<Object>();
    
    @SuppressWarnings({ "unchecked", "rawtypes" })
    public Counter<Object> getInitCounter()
    {
      if (SAME.equals(initParams)) throw new RuntimeException();
      if (INTERNAL.equals(initParams)) return new Counter<Object>(internal);
      else if (ZERO.equals(initParams)) return new Counter<Object>();
      else return (Counter) restoreCounter(new File(initParams));
    }
    @SuppressWarnings({ "unchecked", "rawtypes" })
    public Counter<Object> getCenterCounter()
    {
      if (SAME.equals(reguCenterParams)) return getInitCounter();
      if (INTERNAL.equals(reguCenterParams)) return new Counter<Object>(internal);
      else if (ZERO.equals(reguCenterParams)) return new Counter<Object>();
      else return (Counter) restoreCounter(new File(reguCenterParams));
    }
  }
  
  public static <T> void saveCounter(Counter<T> weights, File file) 
  {
    PrintWriter out = BriefIO.output(file); 
    for (T key : weights)
      out.append(key.toString() + "\t" + weights.getCount(key) + "\n");
    out.close();
  }
  public static Counter<String> restoreCounter(File file) 
  {
    Counter<String> result = new Counter<String>();
    try { for (String line : BriefIO.readLines(file))
        if (!line.equals(""))
        {
          String [] fields = line.split("\\t+");
          String f = fields[0];
          double w = Double.parseDouble(fields[1]);
          if (result.containsKey(f))
            throw new RuntimeException("Duplicate entries for " + f + " in " + file);
          result.setCount(f,w);
        }
    } catch (Exception e) { 
      throw new RuntimeException(e); 
      }
    return result;
  }
  
  public static ExponentialFamily createExpfam(
      MaxentOptions<Object> learningOptions,
      ExponentialFamilyOptions options, 
      FeatureOptions fo, 
      Set<UnorderedPair<SequenceId, SequenceId>> taxaPairs,
      Alphabet alphabet)
  {
    Counter<Object> initParams = 
      options.getInitCounter(),
      centerParams = options.getCenterCounter();
    
    FeatureExtractor fe = new FeatureExtractor(alphabet, taxaPairs, fo);
    Model model = Model.stdBranchSpecificModel(alphabet, fe.getStrTaxonSuffStat());
    ThreeStatesBaseMeasure tsmb = new ThreeStatesBaseMeasure(model);
    
    return new ExponentialFamily(initParams, centerParams, learningOptions, model, tsmb, fe);
  }

  /**
   * Create new natural parameters using penalized likelihood (or other optimization/sampling algorithms)
   * Also: flush the suff stats
   */
  @SuppressWarnings("unchecked")
  public void updateParameters()
  {
    // learn new natural params
    MaxentOptions<Object> currentLearningOptions = MaxentOptions.cloneWithWeights(learningOptions, naturalParams);
    MaxentClassifier<Input,Output,Object> maxentClassifier =
      MaxentClassifier.learnMaxentClassifier(bm, suffStats, featureExtractor,
          currentLearningOptions, regularizationCenters);
    naturalParams = maxentClassifier.weights();
    // create the cached version
    cachedParams = new CachedParams(model, maxentClassifier);
    // flush suff stats
    suffStats = new Counter<LabeledInstance<Input,Output>>();
  }
  public void saveWeightsInExec(String name)
  {
    saveWeights(Results.getFileInResultFolder(name)); 
  }
  public void saveWeights(File f)
  {
    ExponentialFamily.saveCounter(naturalParams, f); 
  }
  
  private void initParameters()
  {
    MaxentClassifier<Input,Output,Object> maxentClassifier =
      MaxentClassifier.createMaxentClassifierFromWeights(bm, naturalParams, featureExtractor);
    // create the cached version
    cachedParams = new CachedParams(model, maxentClassifier);
    // flush suff stats
    suffStats = new Counter<LabeledInstance<Input,Output>>();
  }

  public void addSufficientStatistics( 
      final HetPairHMM pairHMM, final SequenceId topTaxon, final SequenceId botTaxon)
  {
    addSufficientStatistics(this.suffStats, pairHMM, topTaxon, botTaxon);
  }
  public void addSufficientStatistics(final Counter<LabeledInstance<Input,Output>> suffStats, 
      final HetPairHMM pairHMM, final SequenceId topTaxon, final SequenceId botTaxon)
  {
    StrTaxonSuffStat.StrTaxonSuffStatExtractor extractor = model.stSuffStat.getExtractor(pairHMM.str1, pairHMM.str2,topTaxon, botTaxon);
    final double logSumPr = pairHMM.logSumProduct();
    for (int xpos = 0; xpos <= pairHMM.str1.length(); xpos++)
      for (int dx = 0; dx < 2 && xpos + dx <= pairHMM.str1.length(); dx++)
      {
        final int xid = model.charIdAt(pairHMM.str1, xpos, dx);
        for (int ypos = 0; ypos <= pairHMM.str2.length(); ypos++)
          for (int dy = 0; dy < 2 && ypos + dy <= pairHMM.str2.length(); dy++)
            if (dx > 0 || dy > 0)
            {
              final int yid = model.charIdAt(pairHMM.str2, ypos, dy);
              final int gss = extractor.extract(xpos, ypos);
              for (int s1 = 0; s1 < model.nStates; s1++)
                for (int s2 = 0; s2 < model.nStates; s2++)
                {
                  final Input in = new Input(s1, gss, model);
                  final Output out = new Output(s2, xid, yid, model);
                  final double value = Math.exp(pairHMM.logSumProduct(s1, s2, xpos, ypos, dx, dy) - logSumPr);
                  if (value > 0.0)
                    suffStats.incrementCount(new LabeledInstance<Input,Output>(out,in), 
                        value);
                }
            }
      }
  }
  
  public HetPairHMM getReweightedHMM(double [][][] logWeights, Sequence top, Sequence bot, SequenceId topL, SequenceId botL)
  {
    top = top.append(model.BOUNDARY_SYMBOL); 
    bot = bot.append(model.BOUNDARY_SYMBOL); 
    return new HetPairHMM(top, bot, cachedParams.getReweightedHMM(logWeights, top, bot, topL, botL));
  }
  
  public HetPairHMM getHMM(Sequence top, Sequence bot, SequenceId topL, SequenceId botL)
  {
    top = top.append(model.BOUNDARY_SYMBOL); 
    bot = bot.append(model.BOUNDARY_SYMBOL); 
    return new HetPairHMM(top, bot, cachedParams.getUnsupPairHMM(top, bot, topL, botL));
  }
  
  public Counter<Edge> allPairsPosterior(Map<SequenceId,Sequence> sequences)
  {
    Counter<Edge> edgePosteriors = new Counter<Edge>();
    List<SequenceId> langs = new ArrayList<SequenceId>(sequences.keySet());
    for (int i = 0; i < langs.size(); i++)
    {
      final SequenceId l1 = langs.get(i);
      for (int j = i+1; j < langs.size(); j++)
      {
        final SequenceId l2 = langs.get(j);
        final Sequence 
          s1 = sequences.get(l1),
          s2 = sequences.get(l2);
        HetPairHMM hmm = getHMM(s1, s2, l1, l2);
        for (int p1 = 0; p1 < s1.length(); p1++)
          for (int p2 = 0; p2 < s2.length(); p2++)
          {
            final Edge current = new Edge(p1, p2, l1, l2);
            edgePosteriors.setCount(current, Math.exp(hmm.logPosteriorAlignment(p1, p2)));
          }
      }
    }
    return edgePosteriors;
  }
  
  public MSAPoset maxRecallAlignFromAllPairs(Map<SequenceId,Sequence> sequences)
  {
    return MSAPoset.maxRecallMSA(sequences, allPairsPosterior(sequences));
  }
}
