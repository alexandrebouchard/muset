package muset;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

import briefj.BriefIO;
import briefj.CSV;
import briefj.opt.Option;
import briefj.run.Mains;
import briefj.run.Results;



public class ConvertCognateTable implements Runnable
{
  @Option(required = true)
  public File cognates;
  
  @Option(required = true)
  public File features;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new ConvertCognateTable());
  }

  @Override
  public void run()
  {
    Map<String, String> readFeatures = readFeatures();
    
    List<String> result = new ArrayList<String>();
    for (Map<String,String> line : BriefIO.readLines(cognates).indexCSV())
    {
      String languageName = line.get("");
      for (String cognateName : line.keySet())
        if (!"".equals(cognateName) && !"".equals(line.get(cognateName)))
        {
          String word = line.get(cognateName).replaceAll("\\]\\s+.*", "").replace("(", "").replace(")", "");
          word = word.replace("-[", "[").replace("]-", "]");
          List<String> phonemes = new ArrayList<String>();
          for (String phoneme : Splitter.on("]").split(word))
            if (!"".equals(phoneme))
            {
              phoneme = phoneme + "]";
              if (readFeatures.containsKey(phoneme))
                phonemes.add(readFeatures.get(phoneme));
              else
                throw new RuntimeException("Phoneme not found: " + phoneme);
            }
          result.add(CSV.toCSV(cognateName, languageName, Joiner.on(" ").join(phonemes)));
        }
    }
    Collections.sort(result);
    BriefIO.write(Results.getFileInResultFolder("output.csv"), Joiner.on("\n").join(result));
    
  }

  private Map<String,String> readFeatures()
  {
    Map<String,String> result = new LinkedHashMap<String,String>();
    for (Map<String,String> line : BriefIO.readLines(features).indexCSV())
    {
      String features = line.get("features");
      if (result.containsKey(features))
        throw new RuntimeException("Key already there:" + features);
      result.put(features, line.get("input"));
    }
    return result;
  }

}
