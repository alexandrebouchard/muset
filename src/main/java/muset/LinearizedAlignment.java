package muset;


import java.util.HashMap;
import java.util.List;
import java.util.Map;

import briefj.Indexer;

import muset.MSAPoset.Column;




/**
 * An alignment and a linearization of the columns.
 * 
 * Technically, MSAPoset contains a linearization internally, but it has no facilities to set this linearization to
 * a new value.
 * @author bouchard
 *
 */
public class LinearizedAlignment
{
  private final List<Column> columns;
  private final MSAPoset msa;
  @SuppressWarnings({ "rawtypes" })
  public final Map cache = new HashMap();
  
  public LinearizedAlignment(MSAPoset msa)
  {
    this.columns = msa.linearizedColumns();
    this.msa= msa;
  }
  public LinearizedAlignment(MSAPoset msa, List<Column> columns)
  {
    this.columns = columns;
    this.msa = msa;
  }
  public double [][] indicators(SequenceId t, Indexer<Character> indexer, int gapIndex)
  {
    String curString = msa.sequences().get(t);
    double [][] result = new double[nColumns()][indexer.size()+1];
    for (int s = 0; s < nColumns(); s++)
    {
      Map<SequenceId,Integer> points = columns.get(s).getPoints();
      if (points.containsKey(t))
        result[s][indexer.o2i(curString.charAt(points.get(t)))] = 1.0;
      else
        result[s][gapIndex] = 1.0;
    }
    return result;
  }
  public List<Column> getColumns() { return columns; }
  public int nColumns() { return columns.size(); }
  public MSAPoset getMsa()
  {
    return msa;
  }
}