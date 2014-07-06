package muset;


import java.io.File;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import muset.TopoSort.PartialOrder;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.CombinatoricsUtils;

import bayonet.math.NumericalUtils;
import bayonet.math.SpecialFunctions;
import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefMaps;
import briefj.BriefMath;
import briefj.BriefStrings;
import briefj.collections.Counter;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;


/**
 * A set of column (a column is a map from taxon to sequence position) together with
 * a linearization of the columns
 * @author bouchard
 *
 */
public class MSAPoset implements Serializable
{
  private static final long serialVersionUID = 1L;

  private final PartialOrder<Column> poset;
  private TreeSet<Double> linearizedLocations = new TreeSet<Double>();
  private Map<Column,Double> linearizedColumns = Maps.newLinkedHashMap();
  private final Map<SequenceId,Column[]> columnMaps = Maps.newLinkedHashMap();
  private final Map<SequenceId,String> sequences;
  private final List<SequenceId> taxa;
  private boolean linearizationEnabled = true;
  
  public MSAPoset(Map<SequenceId,String> sequences)
  {
    this.taxa = Lists.newArrayList(sequences.keySet());
    Collections.sort(taxa);
    this.sequences = sequences;
    this.poset = new ImplicitPoset(); 
    Counter<Column> fraction = new Counter<Column>();
    
    for (SequenceId lang : taxa)
    {
      Column[] currentMap = new Column[sequences.get(lang).length()];
      columnMaps.put(lang, currentMap);
      for (int i =0 ; i < sequences.get(lang).length(); i++)
      {
        Column currentCol = new Column(lang, i);
        currentMap[i] = (currentCol);
        fraction.setCount(currentCol, -((double) i) / (1.0 + ((double) sequences.get(lang).length())));
      }
    }
    int curIdx = 0;
    for (Column c : fraction)
    {
      Double key = (double) curIdx;
      linearizedColumns.put(c, key);
      linearizedLocations.add(key);
      curIdx++;
    }
  }
  
  public void disableLinearization()
  {
    this.linearizationEnabled = false;
    this.linearizedLocations = null;
  }
  
  public void enableLinearization()
  {
    if (linearizationEnabled)
      return;
    List<Column> linearized = TopoSort.topologicalSort(poset);
    if (linearized == null)
      throw new RuntimeException();
    this.linearizedLocations = new TreeSet<Double>();
    this.linearizedColumns = Maps.newLinkedHashMap(); 
    double idx = 0.0;
    for (Column c : linearized)
    {
      linearizedLocations.add(idx);
      linearizedColumns.put(c, idx);
      idx++;
    }
    this.linearizationEnabled = true;
  }
  
  public final List<SequenceId> taxa() { return Collections.unmodifiableList(taxa); }
  
  public MSAPoset(MSAPoset base)
  {
    this(base.sequences);
    for (Column c : base.linearizedColumns.keySet())
      if (!this.tryAdding(c))
        throw new RuntimeException();
  }
  
  public static MSAPoset maxRecallMSA(Map<SequenceId,String> sequences, Counter<Edge> edgeCounter)
  {
    MSAPoset result = new MSAPoset(sequences);
    for (Edge e : edgeCounter)
      result.tryAdding(e);
    return result;
  }
  
  public static class ROCPoint
  {
    public final double precision, recall, posterior;

    public ROCPoint(double precision, double recall, double posterior)
    {
      this.precision = precision;
      this.recall = recall;
      this.posterior = posterior;
    }
    
  }
  
  public static List<ROCPoint> ROC(
      Map<SequenceId,String> sequences, 
      Counter<Edge> edgeCounter, 
      MSAPoset ref,
      int nPoints)
  {
    MSAPoset currentGuess = new MSAPoset(sequences);
    List<ROCPoint> result = Lists.newArrayList(); 
    // check n points that will be successful
    int totalNPoints = 0;
    for (Edge e : edgeCounter)
      if (currentGuess.tryAdding(e))
        totalNPoints++;
    final int interval = totalNPoints / nPoints;
    currentGuess = new MSAPoset(sequences);
    int current = 0;
    for (Edge e : edgeCounter)
      if (currentGuess.tryAdding(e))
      {
        current++;
        if (current == totalNPoints || current % interval == 0)
          result.add(new ROCPoint(
              edgePrecision(ref, currentGuess),
              edgeRecall(ref, currentGuess),
              edgeCounter.getCount(e)));
      }
    return result;
  }
  
  public int nTaxa() { return sequences.size(); }
  
  public int nEdges()
  {
    int result = 0;
    for (Column c : linearizedColumns.keySet())
      result += CombinatoricsUtils.binomialCoefficient(c.points.size(), 2); 
    return result;
  }
  
  public Collection<Edge> edges()
  {
    List<Edge> result = Lists.newArrayList();
    for (Column c : linearizedColumns.keySet())
      for (int l1i = 0; l1i < taxa.size(); l1i++)
      {
        final SequenceId l1 = taxa.get(l1i);
        if (c.points.containsKey(l1))
          for (int l2i = l1i+1; l2i < taxa.size(); l2i++)
          {
            final SequenceId l2 = taxa.get(l2i);
            if (!l1.equals(l2) && c.points.containsKey(l2))
              result.add(new Edge(c.points.get(l1), c.points.get(l2), l1, l2));
          }
      }
    return result;
  }
  
//  public static MSAPoset parseAlnOrMsfFormats(File f)
//  {
//    return fromMultiAlignmentObject(MultiAlignment.parse(f.getAbsolutePath()));
//  }
//  
//  public static MSAPoset fromMultiAlignmentObject(MultiAlignment ma)
//  {
//    return _fromMultiAlignmentObject(ma, false);
//  }
//  
//  public static MSAPoset coreBlocksFromMultiAlignmentObject(MultiAlignment ma)
//  {
//    return _fromMultiAlignmentObject(ma, true);
//  }
//  
//  private static MSAPoset _fromMultiAlignmentObject(MultiAlignment ma, boolean keepOnlyRef)
//  {
//    MSAPoset result = new MSAPoset(ma.getSequences());
//    for (SequenceCoordinate sc : ma.eqClasses().representatives())
//      if (!keepOnlyRef || sc.isCoreBlock())
//        for (SequenceCoordinate other : ma.eqClasses().eqClass(sc))
//          if (!other.equals(sc))
//          {
//            if (keepOnlyRef && sc.isCoreBlock() != other.isCoreBlock())
//              throw new RuntimeException("The old format, based on annotation files," +
//              		" assumes that columns are either all or all not core block links.");
//            final Edge currentEdge = new Edge(
//                sc.indexInSequence(),  other.indexInSequence(),
//                sc.getNodeIdentifier(),other.getNodeIdentifier());
//            if (!result.tryAdding(currentEdge))
//              throw new RuntimeException();
//          }
//    return result;
//  }
  
  
  
  public Map<SequenceId,String> sequences() 
  {
    return Collections.unmodifiableMap(sequences);
  }
  
  public char charAt(Column c, SequenceId lang)
  {
    return sequences.get(lang).charAt(c.points.get(lang));
  }
  
  public char charAt(Edge e, boolean first)
  {
    return sequences.get(first ? e.lang1() : e.lang2())
      .charAt(first ? e.index1() : e.index2());
  }
  
//  public static void save(MSAPoset msa, File file)
//  {
//    ObjectOutputStream out = IOUtils.openBinOutHard(file);
//    try  {
//      out.writeObject(msa);
//      out.close();
//    } catch (Exception e) { throw new RuntimeException(e); }
//  }
//  
//  public static MSAPoset restore(File filePath) 
//  {
//    try{
//    ObjectInputStream ois = IOUtils.openBinIn(filePath);
//    return (MSAPoset) ois.readObject();
//    } catch (Exception e) { throw new RuntimeException(e); }
//  }
  
  public static boolean isValidSplit(Column c, Set<SequenceId> keepInCurrent)
  {
    if (keepInCurrent.size() == 0 || keepInCurrent.size() == c.points.size())
      return false; // there must be at leas one point in each cc
    if (!c.points.keySet().containsAll(keepInCurrent))
      return false; // malformed request!
    return true;
  }
  
  public void setString(SequenceId t, String s)
  {
    String cur = sequences.get(t);
    if (cur.length() != s.length())
      throw new RuntimeException();
    sequences.put(t, s);
  }
  
  public void fixMSAUsingRandomCharacters(Set<Character> allowed, Random rand)
  {
    List<Character> list = Lists.newArrayList(allowed);
    for (SequenceId t : taxa())
    {
      StringBuilder replacement = new StringBuilder();
      for (char c : sequences().get(t).toCharArray())
        if (allowed.contains(c))
          replacement.append(c);
        else
          replacement.append(list.get(rand.nextInt(list.size())));
      setString(t, replacement.toString());
    }
  }
  
  private Double findNextLocation(Column c)
  {
    // find the new location
    Double currentLoc = linearizedColumns.get(c);
    Iterator<Double> locIterator = linearizedLocations.tailSet(currentLoc).iterator();
    if (!currentLoc.equals(locIterator.next()))
      throw new RuntimeException();
    Double insertLoc = null;
    if (locIterator.hasNext())
    {
      Double nextExistingLoc = locIterator.next();
      insertLoc = currentLoc + (nextExistingLoc - currentLoc) / 2.0;
      if (insertLoc.equals(nextExistingLoc) || insertLoc.equals(currentLoc))
        return null;
    }
    else
      insertLoc = currentLoc + 1.0;
    return insertLoc;
  }
  
  private void recreateLocations()
  {
    final TreeSet<Double> newLinearizedLocations = new TreeSet<Double>();
    final Map<Column,Double> newLinearizedColumns = Maps.newLinkedHashMap();
    
    double curIdx = 0.0;
    for (Column c : linearizedColumns())
    {
      newLinearizedLocations.add(curIdx);
      newLinearizedColumns.put(c, curIdx);
      curIdx++;
    }
    
    this.linearizedColumns = newLinearizedColumns;
    this.linearizedLocations = newLinearizedLocations;
  }
  
  public void split(Column c, Set<SequenceId> keepInCurrent)
  {
    Double insertLoc = null;
    
    if (linearizationEnabled)
    {
      insertLoc = findNextLocation(c);
      if (insertLoc == null)
      {
        recreateLocations();
        insertLoc = findNextLocation(c);
      }
    }
    
    if (!isValidSplit(c, keepInCurrent))
      throw new RuntimeException(); 
    Set<SequenceId> putInNew = Sets.newLinkedHashSet(c.points.keySet());
    putInNew.removeAll(keepInCurrent);
    Column resurrected = new Column(); //resurrect();
    for (SequenceId taxon : putInNew)
      resurrected.points.put(taxon, c.points.get(taxon));
    c.points.keySet().removeAll(putInNew);

    if (linearizationEnabled)
      linearizedLocations.add(insertLoc);
    linearizedColumns.put(resurrected, insertLoc);
    
    for (SequenceId t : resurrected.points.keySet())
    {
      final int strIdx = resurrected.points.get(t);
      columnMaps.get(t)[strIdx] = resurrected;
    }
    
  }
  

  


//  public static int _testArcs(MSAPoset msa)
//  {
//    int nProblems = 0;
//    for (Column c : msa.poset.nodes())
//      for (Column c2 : msa.poset.next(c))
//        if (!intersects(c.getPoints().keySet(), c2.getPoints().keySet()))
//        {
//          System.out.println("Problem:" + c.getPoints() + " and " + c2.getPoints());
//          nProblems++;
//        }
////    System.out.println("N problems: " + nProblems);
//    return nProblems;
//  }
  
  public boolean containsEdge(Edge e)
  {
    return getColumn(e,true) == getColumn(e,false);
  }
  
  public static MSAPoset restrict(MSAPoset msa, Set<SequenceId> taxa)
  {
    Map<SequenceId,String> sequences = new HashMap<SequenceId, String>();
    Set<SequenceId> inter = Sets.intersection(taxa, msa.sequences().keySet());
    for (SequenceId t : inter)
      sequences.put(t, msa.sequences().get(t));
    MSAPoset result = new MSAPoset(sequences);
    for (Column c : msa.columns())
    {
      Column newC = new Column();
      for (SequenceId t : inter)
        if (c.points.containsKey(t))
          newC.points.put(t,c.points.get(t));
      if (!result.tryAdding(newC))
        throw new RuntimeException();
    }
    return result;
  }
  
  public static MSAPoset union(Collection<MSAPoset> collection)
  {
    Map<SequenceId,String> sequences = new HashMap<SequenceId, String>();
    for (MSAPoset msa : collection)
      for (SequenceId key : msa.sequences().keySet())
      {
        String current = sequences.get(key);
        if (current == null)
          sequences.put(key, msa.sequences().get(key));
        else if (!current.equals(msa.sequences().get(key)))
          throw new RuntimeException();
      }
    MSAPoset result = new MSAPoset(sequences);
    for (MSAPoset msa : collection)
      for (Column c : msa.columns())
        if (!result.tryAdding(c))
          return null;
    return result;
  }
  
  /**
   *  try to add all the links in this column
   *  Returns if this was succesful
   *  
   *  Warning: if it fails, the set of columns will be the same as before calling,
   *  but the linearization might be different (but still guaranteed to be consistent)
   */
  public boolean tryAdding(Column c)
  {
    return tryAdding(c.points);
  }
  public boolean tryAdding(Map<SequenceId,Integer> points)
  {
    for (Edge e : spanningEdges(points))
      if (!tryAdding(e))
        return false;
    return true;
  }
  /**
   *  Warning: if it fails, the set of columns will be the same as before calling,
   *  but the linearization might be different (but still guaranteed to be consistent)
   */
  public boolean tryAdding(Edge alignmentLink) { return _tryAdding(alignmentLink, true); }
  /**
   *  Warning: the set of columns will be the same as before calling,
   *  but the linearization might be different (but still guaranteed to be consistent)
   */
  public boolean isValidAddition(Edge alignmentLink)   { return _tryAdding(alignmentLink, false); }
  
  // we call 'arcs' the links in the Hesse diagram, and 'alignmentLink' the links in the multi alignment
  private boolean _tryAdding(Edge alignmentLink, boolean commitChanges)
  {
    Column mergeTo = getColumn(alignmentLink,true), mergeFrom = getColumn(alignmentLink,false);
    if (mergeTo == null || mergeFrom == null)
      throw new RuntimeException();
    if (!mergeTo.disjoint(mergeFrom)) return false;
    
    boolean success = true;
    if (linearizationEnabled)
      for (Pair<Column,Column> currentArc : arcs(mergeFrom))
      {
        // peek to see if we can add it
        final Pair<Column, Column> copyArc = copyArc(currentArc, mergeFrom, mergeTo); 
        if (!(success = TopoSort.onlineTopologicalSort(
            poset, 
            linearizedColumns, 
            copyArc.getLeft(),copyArc.getRight())))
          break;
      }
    if (success && commitChanges)
    {
      mergeTo.points.putAll(mergeFrom.points);
      for (SequenceId t : mergeFrom.points.keySet())
      {
        final int strIdx = mergeFrom.points.get(t);
        final Column[] currentMap = columnMaps.get(t);
        if (currentMap[strIdx] != mergeFrom)
          throw new RuntimeException();
        currentMap[strIdx] = mergeTo;
      }
      Double currentLoc = linearizedColumns.get(mergeFrom);
      linearizedColumns.remove(mergeFrom);
      if (linearizationEnabled)
        linearizedLocations.remove(currentLoc);
    }
    return success;
  }
  
  private static Pair<Column, Column> copyArc(Pair<Column, Column> currentArc,
      Column mergeFrom, Column mergeTo)
  {
    if (currentArc.getLeft() == mergeFrom && currentArc.getRight() == mergeFrom)
      throw new RuntimeException();
    else if (currentArc.getLeft() == mergeFrom)
      return Pair.of(mergeTo, currentArc.getRight());
    else if (currentArc.getRight() == mergeFrom)
      return Pair.of(currentArc.getLeft(), mergeTo);
    else 
      throw new RuntimeException();
  }
  
  public class ImplicitPoset implements PartialOrder<Column>, Serializable 
  {
    private static final long serialVersionUID = 1L;

    @Override
    public Set<Column> next(Column n)
    {
      Set<Column> result = Sets.newLinkedHashSet();
      for (final SequenceId t : n.points.keySet())
      {
        final int curPos = n.points.get(t);
        if (curPos + 1 < sequences().get(t).length())
        {
          // next nhb:
          Column next = columnMaps.get(t)[curPos+1];
          result.add(next);
        }
      }
      return result;
    }

    @Override
    public Set<Column> prev(Column n)
    {
      Set<Column> result = Sets.newLinkedHashSet();
      for (final SequenceId t : n.points.keySet())
      {
        final int curPos = n.points.get(t);
        if (curPos > 0)
        {
          // prev nhb:
          Column prev = columnMaps.get(t)[curPos-1];
          result.add(prev);
        }
      }
      return result;
    }

    @Override
    public Set<Column> nodes()
    {
      return linearizedColumns.keySet();
    }
    
  }
  
  private List<Pair<Column,Column>> arcs(Column mergeFrom)
  {
    List<Pair<Column,Column>> result = Lists.newArrayList();
    for (final SequenceId t : mergeFrom.points.keySet())
    {
      final int curPos = mergeFrom.points.get(t);
      if (curPos > 0)
      {
        // prev nhb:
        Column prev = columnMaps.get(t)[curPos-1];
        result.add(Pair.of(prev, mergeFrom));
      }
      if (curPos + 1 < sequences().get(t).length())
      {
        // next nhb:
        Column next = columnMaps.get(t)[curPos+1];
        result.add(Pair.of(mergeFrom, next));
      }
    }
    return result;
    
  }

  private Column getColumn(Edge alignmentLink, boolean b)
  {
    if (b) return columnMaps.get(alignmentLink.lang1())[alignmentLink.index1()];
    else   return columnMaps.get(alignmentLink.lang2())[alignmentLink.index2()];
  }
  
//  public MultiAlignment toMultiAlignmentObject()
//  {
//    MultiAlignment result = new MultiAlignment(sequences);
//    for (Column c : linearizedColumns.keySet()) //poset.nodes())
//    {
//      for (Edge e : c.spanningEdges())
//        result.addAlign(e.lang1(), e.index1(), e.lang2(), e.index2());
//    }
//    return result;
//  }
  
  private Set<SequenceId> languagesAtOtherEndOfArc(Column reference, Pair<Column,Column> arc)
  {
         if (arc.getLeft() == reference)
      return arc.getRight().points.keySet();
    else if (arc.getRight() == reference)
      return arc.getLeft().points.keySet();
    else
      throw new RuntimeException();
  }

  
  public Column column(SequenceId lang, int index)
  {
    return columnMaps.get(lang)[index];
  }
  
  public Set<Map<SequenceId,Integer>> points()
  {
    Set<Map<SequenceId,Integer>> result = Sets.newLinkedHashSet();
    for (Column c : linearizedColumns.keySet())
      result.add(c.points);
    return result;
  }
  
  public Set<Column> relevantColumns(Set<SequenceId> langs)
  {
    Set<Column> result = Sets.newLinkedHashSet();
    for (Column c : linearizedColumns.keySet())
      if (BriefCollections.intersects(c.points.keySet(), langs))
        result.add(c);
    return result;
  }
  
  public static List<Edge> spanningEdges(Map<SequenceId,Integer> points)
  {
    List<Edge> result = Lists.newArrayList();
    if (points.size() < 2)
      return result;
    List<SequenceId> langs = Lists.newArrayList(points.keySet());
    final SequenceId baseLang = langs.get(0);
    final int basePos = points.get(baseLang);
    for (int i = 1; i < langs.size(); i++)
    {
      final SequenceId otherLang = langs.get(i);
      result.add(new Edge(basePos, points.get(otherLang), baseLang, otherLang));
    }
    return result;
  }
  
  public boolean isFull(Column c)
  {
    return c.points.keySet().equals(sequences.keySet());
  }

  public static class Column implements Serializable
  {
    private static final long serialVersionUID = 1L;
    private Map<SequenceId,Integer> points;
    
    /**
     * a map from language to a POSITION in the string
     * NOT an indexed character!
     * @return
     */
    public Map<SequenceId,Integer> getPoints() 
    { 
      return Collections.unmodifiableMap(points); 
    }
    
    public Column() 
    {
      points = new HashMap<SequenceId, Integer>(2, 0.75f);
    }
    public Column(SequenceId lang, int i)
    {
      points = new HashMap<SequenceId, Integer>(2, 0.75f);
      points.put(lang,i);
    }
    
    public List<Edge> spanningEdges()
    {
      return MSAPoset.spanningEdges(this.points);
    }
    
    private boolean disjoint(Column c2)
    {
      return Collections.disjoint(this.points.keySet(), c2.points.keySet());
    }
    // note: should NOT add a different .equal, .hashcode
    @Override public boolean equals(Object obj) { return super.equals(obj); }
    @Override public int hashCode() { return super.hashCode(); }
    public String toString() { return points.toString(); }
  }
  
  public Collection<Column> columns()
  {
    return Collections.unmodifiableCollection(linearizedColumns.keySet());
  }
  
  public List<Column> linearizedColumns()
  {
    if (!linearizationEnabled)
      throw new RuntimeException();
    Column [] all = new Column[linearizedLocations.size()];
    Map<Double,Integer> conversion = Maps.newLinkedHashMap();
    int cur = 0;
    for (Double d : linearizedLocations)
      conversion.put(d, cur++);
    for (Column c : linearizedColumns.keySet())
    {
      Double key = linearizedColumns.get(c);
      all[conversion.get(key)] = c;
    }
    return Arrays.asList(all);
  }
  
  public double getIdentityStatistic() { return basicStat(true); }
  public double getAlignedStatistics() { return basicStat(false);}
  
  public double getMeanSequenceLength() 
  {
    SummaryStatistics result = new SummaryStatistics();
    for (String seq : sequences().values())
      result.addValue(seq.length());
    return result.getMean();
  }
  
  private double basicStat(boolean isIdentity)
  {
    SummaryStatistics stat = new SummaryStatistics();
    for (SequenceId t1 : taxa())
      for (SequenceId t2 : taxa())
        if (!t1.equals(t2))
          stat.addValue(basicStat(t1,t2,isIdentity));
    return stat.getMean();
  }
  
  private double basicStat(SequenceId t1, SequenceId t2, boolean isIdentity) // o.w., just check aligned
  {
    SummaryStatistics stat = new SummaryStatistics();
    for (int i = 0; i < sequences().get(t1).length(); i++)
    {
      Column c = column(t1, i);
      if(!c.points.containsKey(t2))
        stat.addValue(0.0);
      else if (!isIdentity)
        stat.addValue(1.0);
      else
        stat.addValue( charAt(c, t1) == charAt(c, t2) ? 1.0 : 0.0);
    }
    return stat.getMean();
  }

  @Override
  public String toString() { return toString(null); }
  public String toString(Set<SequenceId> restriction)
  {
    StringBuilder [] builders = createPaddedStrings(restriction);
    StringBuilder result = new StringBuilder();
    List<SequenceId> languagePrintOrder = printOrder(restriction);
    for (int i = 0; i < languagePrintOrder.size(); i++)
    {
      result.append(builders[i]);
      result.append("|");
      result.append(languagePrintOrder.get(i));
      result.append("\n");
    }
    return result.toString();
  }
  
  
  private List<SequenceId> printOrder(Set<SequenceId> restriction)
  {
    List<SequenceId> _printOrder = new ArrayList<SequenceId>(
        restriction == null ? 
          sequences.keySet() :
          Sets.intersection(restriction,sequences.keySet()));
    Collections.sort(_printOrder);
    return _printOrder;
  }
  
  private StringBuilder [] createPaddedStrings(Set<SequenceId> restriction)
  {
    List<SequenceId> printOrder = printOrder(restriction);
    Map<SequenceId,Integer> languagePrintOrder = invert(printOrder);
    StringBuilder [] builders = new StringBuilder[printOrder.size()];
    
    for (SequenceId lang : languagePrintOrder.keySet())
    {
      final int row = languagePrintOrder.get(lang);
      StringBuilder current = new StringBuilder();
      builders[row] = current;

    }
    for (Column c : linearizedColumns())
      if (restriction == null || BriefCollections.intersects(restriction, c.points.keySet()))
      {
        for (SequenceId lang : languagePrintOrder.keySet())
        {
          String currentChar = c.points.keySet().contains(lang)   ?
              "" + sequences.get(lang).charAt(c.points.get(lang)) :
              "-";
          builders[languagePrintOrder.get(lang)].append(currentChar);
        }
      }
    return builders;
  }
  
  public static <T> Map<T,Integer> invert(List<T> list)
  {
    Map<T,Integer> result = new HashMap<T,Integer>();
    for (int i =0; i < list.size(); i++)
      result.put(list.get(i),i);
    return result;
  }
  public static <K,V> Map<V,Set<K>> invert(Map<K,V> map)
  {
    Map<V,Set<K>> result = new HashMap<V,Set<K>>();
    for (K key : map.keySet())
      BriefMaps.getOrPutSet(result, map.get(key)).add(key);
    return result;
  }
  
  public static double columnRecall(MSAPoset gold, MSAPoset guess)
  {
    Set<Map<SequenceId,Integer>> goldColumns = gold.points(),
      guessColumns = guess.points();
    double num = 0.0, denom = 0.0;
    for (Map<SequenceId,Integer> goldColumn : goldColumns)
      if (goldColumn.size() > 1) // in bali references, this will have the effect of only counting core blocks
      {
        denom++;
        if (guessColumns.contains(goldColumn))
          num++;
      }
    return num/denom;
  }
  
  public static double edgeRecall(MSAPoset gold, MSAPoset guess)
  {
    double num = 0.0, denom = 0.0;
    for (Edge e : gold.edges())
    {
      denom++;
      if (guess.containsEdge(e))
        num++;
    }
    return num/denom;
  }
  
  public static double edgePrecision(MSAPoset gold, MSAPoset guess)
  {
    return edgeRecall(guess, gold);
  }
  
  public static double edgeF1(MSAPoset gold, MSAPoset guess)
  {
    return f1Score(edgePrecision(gold, guess), edgeRecall(gold, guess));
  }
  
  public static double f1Score(double precision, double recall)
  {
    if (precision + recall == 0.0) return 0.0;
    return 2 * (precision * recall) / (precision + recall);
  }
  
  /**
   * keeps only capitalized links, then capitalize everything
   * @param msa
   * @return
   */
  public static MSAPoset processBenchmarkReference(MSAPoset msa)
  {  
    MSAPoset result = new MSAPoset(msa.sequences);
    result.disableLinearization();
    for (Column c : msa.linearizedColumns.keySet())
    {
      Map<SequenceId,Integer> processed = Maps.newLinkedHashMap();
      for (SequenceId item : c.points.keySet())
        if (Character.isUpperCase(msa.charAt(c, item)))
          processed.put(item, c.points.get(item));
        result.tryAdding(processed);
    }
    for (SequenceId lang : result.taxa())
      result.sequences.put(lang, result.sequences.get(lang).toUpperCase());
    result.enableLinearization();
    return result;
  }
//  /**
//   * keeps only capitalized links, then capitalize everything
//   * @param msa
//   * @return
//   * @deprecated
//   */
//  public static MSAPoset _processBenchmarkReference(MSAPoset msa)
//  {  
//    System.out.println("Warning! Using legacy implementation!");
//    MSAPoset result = new MSAPoset(msa.sequences);
//    for (Edge e : msa.edges())
//      if (Character.isUpperCase(msa.charAt(e, true)) &&
//          Character.isUpperCase(msa.charAt(e, false)))
//        result.tryAdding(e);
//    for (Taxon lang : result.taxa())
//      result.sequences.put(lang, result.sequences.get(lang).toUpperCase());
//    return result;
//  }
  
//  public static class SaveMSAPoset implements Runnable
//  {
//    @Option public File path = null;
//    public static void main(String [] args)
//    {
//      IO.runLight(args, new SaveMSAPoset());
//    }
//
//    @Override
//    public void run()
//    {
//      int n = Integer.MAX_VALUE;
//      MSAPoset msa = 
//        PreprocessGutellData.randomDataSet(new File("/Users/bouchard/Documents/data/gutell/16S.3.alnfasta"), 1, n, new Random(1)).get(0);
//      File outFile = new File(path.getAbsolutePath() + ".bin");
//      
//      System.gc();System.gc();System.gc();System.gc();System.gc();
//      System.out.println("" + n + "\t" + Runtime.getRuntime().totalMemory()/1024.0/1024.0 + "");
//      
//      MSAPoset.save(msa, outFile);
//      
//      System.out.println("Done");
//    }
//  }
//  

//  public static void main(String [] args)
//  {
//    
////    MSAPoset read = parseAlnOrMsfFormats(new File(args[0]));
////    
////    System.out.println(read);
////    
////    if (true) return;
////    
////    {
//////      PreprocessG(new File("/Users/bouchard/Downloads/5S.3.alnfasta"));
//////      System.out.println(msa);
//////      if (true) return;
////    }
//    
////    {
////      MSAPoset msa = parseFASTA(new File("/Users/bouchard/w/evolvere/data/bench1.0/bali2dna/ref/1aab_ref1"));
//////      System.out.println(msa);
//////      msa = keepOnlyEdgesBetweenCapitalizedSymbols(msa);
//////      System.out.println(msa);
//////      msa = capitalize(msa);
//////      System.out.println(msa);
//////      if (true) return; 
////    }
////    
//    
////    {
////      BalibaseCorpusOptions baliopt = new BalibaseCorpusOptions();
////      baliopt.referenceAlignmentsPath.clear();
////      File path = new File("data/BAliBASE/ref1/test1/");
////      baliopt.referenceAlignmentsPath.add(path.getAbsolutePath());
////      path = new File("data/BAliBASE/ref1/test2/");
////      baliopt.referenceAlignmentsPath.add(path.getAbsolutePath());
////      path = new File("data/BAliBASE/ref1/test3/");
////      baliopt.referenceAlignmentsPath.add(path.getAbsolutePath());
////      path = new File("data/BAliBASE/ref2/test/");
////      baliopt.referenceAlignmentsPath.add(path.getAbsolutePath());
////      path = new File("data/BAliBASE/ref3/test/");
////      baliopt.referenceAlignmentsPath.add(path.getAbsolutePath());
////      path = new File("data/BAliBASE/ref4/test/");
////      baliopt.referenceAlignmentsPath.add(path.getAbsolutePath());
////      path = new File("data/BAliBASE/ref5/test/");
////      baliopt.referenceAlignmentsPath.add(path.getAbsolutePath());
//////      for (String arg :args)
//////        baliopt.referenceAlignmentsPath.add(arg);
////      final BalibaseCorpus bc = new BalibaseCorpus(baliopt);
////      Random rand = new Random(1);
////      for (CognateId id : bc.intersectedIds())
////      {
////        MSAPoset gold = MSAPoset.coreBlocksFromMultiAlignmentObject(bc.getMultiAlignment(id));
////        System.out.println(id);
////        System.out.println("Gold:\n" + gold);
////        // construct a random align
////        MSAPoset randomGuess = new MSAPoset(gold.sequences());
////        for (int i = 0; i < 100; i++)
////        {
////          Taxon l1 = Sampling.randomElt(gold.taxa(), rand),
////                   l2 = Sampling.randomElt(gold.taxa(), rand);
////          int i1= rand.nextInt(gold.sequences().get(l1).length()),
////              i2= rand.nextInt(gold.sequences().get(l2).length());
////          randomGuess.tryAdding(new Edge(i1,i2,l1,l2));
////          randomGuess.tryAdding(Sampling.randomElt(gold.edges(), rand));
////        }
////        System.out.println("Random guess:\n" + randomGuess); double oldv, newv;
////        System.out.println("SP (old way): " + (oldv=gold.toMultiAlignmentObject().sumOfPairsScore(randomGuess.toMultiAlignmentObject())));
////        System.out.println("SP (new way): " + (newv=edgeRecall(gold, randomGuess)));
////        if (!MathUtils.close(oldv,newv))
////          throw new RuntimeException();
////      }
////    }
//    
//    
////    System.out.println(restore(new File(args[0])));
//    int maxL = 9;
//    Random rand = new Random(1);
//    int nLangs = 10;
//    
////    int maxL = 3;
////    Random rand = new Random(1);
////    int nLangs = 3;
//    
//    Map<Taxon,String> seqs= Maps.newLinkedHashMap();
//    for (int i =0 ; i < nLangs; i++)
//    {
//      Taxon l = new Taxon("l" + i);
//      String word = "";
//      for (int w = 0; w < maxL; w++)
//        word += w;
//      seqs.put(l, word);
//    }
//    MSAPoset msa = new MSAPoset(seqs);
//    System.out.println(msa);
//    List<Taxon> list = new ArrayList<Taxon>(seqs.keySet());
//    int id = 0;
//    for (int i = 0; i < 1000000; i ++)
//    {
////      System.out.println(msa.linearizedLocations);
////      new MSAPoset(msa);
////      if (i==4)
////        System.out.println("debug");
//      Taxon l1 = list.get(rand.nextInt(list.size())),
//               l2 = list.get(rand.nextInt(list.size()));
//      int i1 = rand.nextInt(seqs.get(l1).length()),
//          i2 = rand.nextInt(seqs.get(l2).length());
//      Edge e = new Edge(i1,i2,l1,l2);
//      if (i % 10000 == 0)
//        System.out.println(i+" Trying to add edge: " + e);
//      
//      MultiAlignment bu = msa.toMultiAlignmentObject();
//      boolean test = msa.isValidAddition(e);
//      if (!bu.equals(msa.toMultiAlignmentObject()))
//        throw new RuntimeException();
//      
//      boolean success = msa.tryAdding(e);
//      if (_testArcs(msa) > 0)
//        throw new RuntimeException("Sanity failed after insert. ID:" + id);
//      id++;
//      if (success != test)
//        throw new RuntimeException();
//    if (i % 10000 == 0)
//    {
//      System.out.println("Success: " + success);
//
//        System.out.println("New align:\n" + msa);
//      System.out.println();
//    }
//      
//        
//      if (i % 2 == 0)
//      {
//        List<Column> cols = Lists.newArrayList();
//        for (Column c : msa.linearizedColumns.keySet())
//          if (c.points.size() > 1)
//            cols.add(c);
//        if (cols.size() > 1)
//        {
//          Column c = cols.get(rand.nextInt(cols.size()));
//          Set<Taxon> toKeep = Sets.newLinkedHashSet(); 
//          Queue<Taxon> available = new LinkedList<Taxon>(c.points.keySet());
//          toKeep.add(available.poll());
//          available.poll();
//          for (Taxon lang : available)
//            if (rand.nextBoolean())
//              toKeep.add(lang);
//          if (i % 10000 == 0)
//            System.out.println("Splitting(" + c.points + "," + toKeep+")");
//          msa.split(c, toKeep);
//          if (_testArcs(msa) > 0)
//            throw new RuntimeException("Sanity failed after split");
//          if (i % 10000 == 0)
//          {
//          System.out.println("New align:\n" + msa);
//          System.out.println();
//          }
//        }
//      }
//    }
//  }
  
  public void toFASTA(File f)
  {
    PrintWriter out = BriefIO.output(f); //IOUtils.openOutHard(f);
    for (SequenceId t : taxa)
    {
      final String curSeq = sequences.get(t);
      out.append(">" + t + "\n");
      for (Column c : linearizedColumns())
        if (c.points.containsKey(t))
          out.append("" + curSeq.charAt(c.points.get(t)));
        else
          out.append("-");
      out.append("\n");
    }
    out.close();
  }
  

  public static MSAPoset parseFASTA(File f)
  {
    // gather data
    Map<SequenceId,List<Integer>> alignData = Maps.newLinkedHashMap();
    Map<SequenceId,StringBuilder> stringData = Maps.newLinkedHashMap();
    Counter<SequenceId> highestIndex = new Counter<SequenceId>();
    SequenceId currentTaxon = null;
    for (String line : BriefIO.readLines(f))
      if (line.matches("^\\s*$") || line.matches("^[;].*"))
        ;
      else if (line.matches("[>].*"))
      {
        currentTaxon = new SequenceId( BriefStrings.firstGroupFromFirstMatch("[>](.*)",line));
        if (alignData.containsKey(currentTaxon))
          throw new RuntimeException("Duplicated taxon name:" + currentTaxon);
        alignData.put(currentTaxon, new ArrayList<Integer>());
        stringData.put(currentTaxon, new StringBuilder());
      }
      else if (line.matches("[a-zA-Z.-]*"))
      {
        if (currentTaxon == null)
          throw new RuntimeException("Sequences should be preceded by the taxon name using a line " +
          		"of the form \">[name]\"");
        for (char c : line.toCharArray())
          if (c == '.' || c == '-')
            alignData.get(currentTaxon).add(null);
          else
          {
            List<Integer> currentAlignData = alignData.get(currentTaxon);
            currentAlignData.add((int) (highestIndex.getCount(currentTaxon)));
            highestIndex.incrementCount(currentTaxon, 1.0);
            stringData.get(currentTaxon).append(c);
          }
      }
      else
        throw new RuntimeException("Invalid line:" + line);
    // construct the align
    Map<SequenceId,String> strings = Maps.newLinkedHashMap();
    for (SequenceId lang :  stringData.keySet())
      strings.put(lang, stringData.get(lang).toString());
    MSAPoset result = new MSAPoset(strings);
    result.disableLinearization();
    int len = -1;
    for (SequenceId lang : alignData.keySet())
    {
      if (len == -1)
        len = alignData.get(lang).size();
      else if (len != alignData.get(lang).size())
        throw new RuntimeException("Invalid alignment spec: " +
        		"all gap-padded seqns should have the len");
//      if (len != stringData.get(lang).length())
//        throw new RuntimeException();
    }
    
//    System.out.println("Got here2");
    
    
    List<SequenceId> langs = Lists.newArrayList(alignData.keySet());
    for (int p = 0; p < len; p++)
    {
      Map<SequenceId,Integer> points = Maps.newLinkedHashMap();
      for (int l = 0; l < langs.size(); l++)
      {
        SequenceId l1 = langs.get(l);
        if (alignData.get(l1).get(p) != null)
          points.put(l1, alignData.get(l1).get(p));
      }
      result.tryAdding(points);
    }
    result.enableLinearization();
    return result;
  }
  

  public PartialOrder<Column> getPoset()
  {
    return poset; // an implicit representation
  }

  public static boolean deepEquals(MSAPoset msa1, MSAPoset msa2)
  {
    if (!msa1.sequences().equals(msa2.sequences())) return false;
    return Sets.newLinkedHashSet(msa1.edges()).equals(Sets.newLinkedHashSet(msa2.edges()));
  }

}
