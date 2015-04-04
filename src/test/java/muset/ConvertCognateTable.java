package muset;

import briefj.opt.Option;
import briefj.run.Mains;



public class ConvertCognateTable implements Runnable
{
  @Option(required = true)
  public File inputCsv;

  public static void main(String[] args)
  {
    Mains.instrumentedRun(args, new ConvertCognateTable());
  }

  @Override
  public void run()
  {
    
  }

}
