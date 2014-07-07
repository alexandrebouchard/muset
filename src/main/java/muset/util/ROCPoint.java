package muset.util;



public class ROCPoint
{
  public final double precision, recall, posterior;

  public ROCPoint(double precision, double recall, double posterior)
  {
    this.precision = precision;
    this.recall = recall;
    this.posterior = posterior;
  }
  
}