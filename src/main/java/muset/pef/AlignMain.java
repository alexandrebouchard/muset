package muset.pef;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import bayonet.regression.LabeledInstance;
import bayonet.regression.MaxentClassifier.MaxentOptions;
import briefj.BriefIO;
import briefj.BriefMaps;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;
import muset.MSAPoset;
import muset.SequenceId;
import muset.hmm.HetPairHMM;
import muset.pef.ExponentialFamily.ExponentialFamilyOptions;
import muset.pef.FeatureExtractor.FeatureOptions;
import muset.util.Edge;



public class AlignMain implements Runnable
{
  @OptionSet(name = "data")
  public SequenceDataset dataset = new SequenceDataset();
  
  @Option
  public int nIterations = 10;

  @OptionSet(name = "maxent")
  @SuppressWarnings("rawtypes")
  public MaxentOptions learningOptions = new MaxentOptions();

  @OptionSet(name = "learn")
  public ExponentialFamilyOptions learnOptions = new ExponentialFamilyOptions();
  
  @OptionSet(name = "features")
  public FeatureOptions fo = new FeatureOptions();
  
  private ExponentialFamily learnedModel;
  private File iterationSpecificOutput;

  public static class SequenceDataset
  {
    @Option(required = true)
    public File csvFile;
    
    private Map<GroupId,Map<SequenceId,String>> data = null;
    private Indexer<Character> indexer = null;
    
    public Collection<GroupId> groupIds()
    {
      ensureLoaded();
      return data.keySet();
    }
    
    public Indexer<Character> getIndexer()
    {
      ensureLoaded();
      return indexer;
    }
    
    public Map<SequenceId,String> getSequences(GroupId id)
    {
      ensureLoaded();
      return data.get(id);
    }

    private void ensureLoaded()
    {
      if (data != null)
        return;
      data = new LinkedHashMap<AlignMain.GroupId, Map<SequenceId,String>>();
      indexer = new Indexer<Character>();
      for (List<String> line : BriefIO.readLines(csvFile).splitCSV())
      {
        if (line.size() != 3)
          throw new RuntimeException("There should be 3 fields: the group, the taxon, the string");
        
        String seq = line.get(2);
        BriefMaps.getOrPutMap(data, new GroupId(line.get(0))).put(new SequenceId(line.get(1)), seq);
        for (char c : seq.toCharArray())
          indexer.addToIndex(c);
      }
    }

    public Set<UnorderedPair<SequenceId, SequenceId>> taxaPairs()
    {
      ensureLoaded();
      Set<UnorderedPair<SequenceId, SequenceId>> result = new LinkedHashSet<UnorderedPair<SequenceId,SequenceId>>();
      for (Map<SequenceId, String> datum : data.values())
        for (SequenceId key1 : datum.keySet())
          for (SequenceId key2 : datum.keySet())
            if (!key1.equals(key2))
              result.add(UnorderedPair.of(key1,key2));
      return result;
    }
  }
  
  public static class GroupId
  {
    private final String groupName;

    @Override
    public String toString()
    {
      return "GroupId [groupName=" + groupName + "]";
    }

    private GroupId(String groupName)
    {
      super();
      this.groupName = groupName;
    }

    @Override
    public int hashCode()
    {
      final int prime = 31;
      int result = 1;
      result = prime * result
          + ((groupName == null) ? 0 : groupName.hashCode());
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
      GroupId other = (GroupId) obj;
      if (groupName == null)
      {
        if (other.groupName != null)
          return false;
      } else if (!groupName.equals(other.groupName))
        return false;
      return true;
    }
    
    
  }
  
  private void doIteration(int iterationNumber)
  {
    iterationSpecificOutput = (Results.getFolderInResultFolder("iteration_" + iterationNumber));
    
    for (GroupId groupId : dataset.groupIds())
      align(groupId, dataset.getSequences(groupId));
    learnedModel.updateParameters();
  }

  private void align(GroupId groupId, Map<SequenceId, String> datum)
  {
    Counter<Edge> edgePosteriors = new Counter<Edge>();
    // get all pairs
    List<Map<SequenceId,String>> pairs = pairs(datum);
    Counter<LabeledInstance<Input,Output>> suffStats = new Counter<LabeledInstance<Input,Output>>();
    // align them
    for (Map<SequenceId,String> pair : pairs)
    {
      List<SequenceId> ids = new ArrayList<SequenceId>(pair.keySet());
      SequenceId topId = ids.get(0), botId = ids.get(1);
      String top = pair.get(topId), bot = pair.get(botId);
      HetPairHMM hmm = learnedModel.getHMM(top, bot, topId, botId);
      learnedModel.addSufficientStatistics(suffStats, hmm, topId, botId);
      for (int botPos = 0; botPos < bot.length(); botPos++)
        for (int topPos = 0; topPos < top.length(); topPos++)
          edgePosteriors.setCount(
              new Edge(topPos, botPos, topId, botId), 
              Math.exp(hmm.logPosteriorAlignment(topPos, botPos)));
    }
    
    // update suff stats
    learnedModel.suffStats.incrementAll(suffStats);
    
    // create consensus align
    MSAPoset maxRecallMSA = MSAPoset.maxRecallMSA(datum, edgePosteriors);
    
    BriefIO.write(new File(iterationSpecificOutput, "" + groupId + ".align.txt") , maxRecallMSA.toString());
  }

  private List<Map<SequenceId, String>> pairs(Map<SequenceId, String> datum)
  {
    List<SequenceId> ids = new ArrayList<SequenceId>(datum.keySet());
    List<Map<SequenceId, String>> result = new ArrayList<Map<SequenceId,String>>();
    for (int i = 0; i < ids.size(); i++)
      for (int j = i + 1; j < ids.size(); j++)
      {
        Map<SequenceId,String> pair = new LinkedHashMap<SequenceId, String>();
        pair.put(ids.get(i), datum.get(ids.get(i)));
        pair.put(ids.get(j), datum.get(ids.get(j)));
      }
    return result;
  }

  @SuppressWarnings("unchecked")
  @Override
  public void run()
  {
    learnedModel = ExponentialFamily.createExpfam(learningOptions , learnOptions, fo, dataset.taxaPairs(), dataset.getIndexer());
    
    for (int iter = 0; iter < nIterations; iter++)
      doIteration(iter);
  }
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new AlignMain());
  }
  
}
