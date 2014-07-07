package muset.util;


import java.io.Serializable;
import java.util.List;
import java.util.regex.Pattern;

import muset.SequenceId;

import org.apache.commons.lang3.tuple.Pair;

import briefj.BriefStrings;
import briefj.collections.UnorderedPair;





/**
 * An edge connects two alignment points
 * An alignment point is characterized by a sequenceId and an index in the sequence
 */
public final class Edge  implements Serializable
{
  private static final long serialVersionUID = 1L;
  public int index1()     { return data.getFirst ().getLeft (); }
  public int index2()     { return data.getSecond().getLeft (); }
  public SequenceId sequenceId1() { return data.getFirst ().getRight(); }
  public SequenceId sequenceId2() { return data.getSecond().getRight(); }
  private final UnorderedPair<Pair<Integer,SequenceId>,Pair<Integer,SequenceId>> data;
  public Edge(int index1, int index2, SequenceId sequenceId1,
      SequenceId sequenceId2)
  {
    data = new UnorderedPair<Pair<Integer,SequenceId>,Pair<Integer,SequenceId>>(
        Pair.of(index1,sequenceId1),Pair.of(index2,sequenceId2));
  }
  @Override public String toString() { return data.toString(); }
  public static Pattern p = Pattern.compile("\\(\\(([^,]*), ([^)]*)\\), \\(([^,]*), ([^)]*)\\)\\)");
  public static Edge fromString(String str)
  {
    List<String> matches = BriefStrings.allGroupsFromFirstMatch(p, str);
    int 
      i1 = Integer.parseInt(matches.get(0)),
      i2 = Integer.parseInt(matches.get(2));
    SequenceId 
      l1 = new SequenceId(matches.get(1)),
      l2 = new SequenceId(matches.get(3));
    return new Edge(i1,i2,l1,l2);
  }
  @Override public int hashCode() { return data.hashCode(); }
  @Override public boolean equals(Object obj) 
  { return ((Edge) obj).data.equals(this.data); }
}
