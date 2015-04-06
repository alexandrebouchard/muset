package muset;


import java.io.File;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import muset.Alphabet.Letter;
import muset.util.Edge;
import muset.util.ROCPoint;
import muset.util.SupportFunctions;
import muset.util.TopoSort;
import muset.util.TopoSort.PartialOrder;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.CombinatoricsUtils;

import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefStrings;
import briefj.collections.Counter;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;


/**
 * A modifiable multiple sequence alignment (MSA).
 * 
 * Uses the algorithm described in 
 * http://bioinformatics.oxfordjournals.org/content/23/2/e24.full
 * to efficiently check that proposed additions
 * of links still yield a valid MSA.
 * 
 * @author bouchard
 *
 */
public class MSAPoset implements Serializable
{
  private static final long serialVersionUID = 1L;

  private final PartialOrder<Column> poset;
  private TreeSet<Double> linearizedLocations = new TreeSet<Double>();
  Map<Column,Double> linearizedColumns = Maps.newLinkedHashMap();
  private final Map<SequenceId,Column[]> columnMaps = Maps.newLinkedHashMap();
  final Map<SequenceId,Sequence> sequences;
  private final List<SequenceId> allSequenceIds;
  private boolean linearizationEnabled = true;
  
  /**
   * Create a new empty alignment over the given sequences.
   * 
   * @param sequences
   */
  public MSAPoset(Map<SequenceId,Sequence> sequences)
  {
    this.allSequenceIds = Lists.newArrayList(sequences.keySet());
    Collections.sort(allSequenceIds);
    this.sequences = sequences;
    this.poset = new ImplicitPoset(); 
    Counter<Column> fraction = new Counter<Column>();
    
    for (SequenceId sequenceId : allSequenceIds)
    {
      Column[] currentMap = new Column[sequences.get(sequenceId).length()];
      columnMaps.put(sequenceId, currentMap);
      for (int i =0 ; i < sequences.get(sequenceId).length(); i++)
      {
        Column currentCol = new Column(sequenceId, i);
        currentMap[i] = (currentCol);
        fraction.setCount(currentCol, -((double) i) / (1.0 + ((double) sequences.get(sequenceId).length())));
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
  
  void disableLinearization()
  {
    this.linearizationEnabled = false;
    this.linearizedLocations = null;
  }
  
  void enableLinearization()
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
  
  /**
   * 
   * @return All the sequence ids contained in this MSA
   */
  public final List<SequenceId> allSequenceIds() { return Collections.unmodifiableList(allSequenceIds); }
  
  /**
   * Performs a deep copy.
   * 
   * @param base
   */
  public MSAPoset(MSAPoset base)
  {
    this(base.sequences);
    for (Column c : base.linearizedColumns.keySet())
      if (!this.tryAdding(c))
        throw new RuntimeException();
  }
  
  /**
   * Create an alignment by greedily adding, if possible, the links (edges) provided 
   * by the counter, in order of highest to lowest score. The efficient poset representation
   * is used to ensure that this yields a valid alignment (if adding one would break the
   * MSA properties, then it is just skipped).
   * 
   * This is called max recall as no attention is spent on precision under this strategy.
   * 
   * @param sequences
   * @param edgeCounter
   * @return
   */
  public static MSAPoset maxRecallMSA(Map<SequenceId,Sequence> sequences, Counter<Edge> edgeCounter)
  {
    MSAPoset result = new MSAPoset(sequences);
    for (Edge e : edgeCounter)
      result.tryAdding(e);
    return result;
  }
  
  /**
   * Creates a full ROC path of alignment quality scores by adding 
   * alignment links (edges) in order of score (provided by the counter).
   * 
   * @see maxRecallMSA, this is essentially a more instrumented version.
   * 
   * @param sequences
   * @param edgeCounter
   * @param ref
   * @param nPoints
   * @return
   */
  public static List<ROCPoint> ROC(
      Map<SequenceId,Sequence> sequences, 
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
              Metrics.edgePrecision(ref, currentGuess),
              Metrics.edgeRecall(ref, currentGuess),
              edgeCounter.getCount(e)));
      }
    return result;
  }
  
  /**
   * The number of sequences in the MSA (e.g. for pairwise, this is 2).
   * @return
   */
  public int nSequences() { return sequences.size(); }
  
  /**
   * The total number of edges, obtained by multiplying the number of MSA
   * columns by the binomial coefficient of the number of non-gap symbol of
   * each column. 
   * 
   * @return
   */
  public int nEdges()
  {
    int result = 0;
    for (Column c : linearizedColumns.keySet())
    {
      int nPts = c.points.size();
      result += nPts < 2 ? 0 : CombinatoricsUtils.binomialCoefficient(nPts, 2); 
    }
    return result;
  }
  
  /**
   * The list of all edges. 
   * 
   * @see nEdges()
   * 
   * @return
   */
  public Collection<Edge> edges()
  {
    List<Edge> result = Lists.newArrayList();
    for (Column c : linearizedColumns.keySet())
      for (int l1i = 0; l1i < allSequenceIds.size(); l1i++)
      {
        final SequenceId l1 = allSequenceIds.get(l1i);
        if (c.points.containsKey(l1))
          for (int l2i = l1i+1; l2i < allSequenceIds.size(); l2i++)
          {
            final SequenceId l2 = allSequenceIds.get(l2i);
            if (!l1.equals(l2) && c.points.containsKey(l2))
              result.add(new Edge(c.points.get(l1), c.points.get(l2), l1, l2));
          }
      }
    return result;
  }

  /**
   * The original sequences (without gap symbols).
   * 
   * @return
   */
  public Map<SequenceId,Sequence> sequences() 
  {
    return Collections.unmodifiableMap(sequences);
  }
  
  /**
   * Access a character in a column.
   * 
   * @param c
   * @param sequenceId
   * @return
   */
  public Letter charAt(Column c, SequenceId sequenceId)
  {
    return sequences.get(sequenceId).letterAt(c.points.get(sequenceId));
  }
 
  private static boolean isValidSplit(Column c, Set<SequenceId> keepInCurrent)
  {
    if (keepInCurrent.size() == 0 || keepInCurrent.size() == c.points.size())
      return false; // there must be at least one point in each cc
    if (!c.points.keySet().containsAll(keepInCurrent))
      return false; // malformed request!
    return true;
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
  
  /**
   * Break (split) a column into two columns.
   * 
   * Modifies this object in place.
   * 
   * @param c The current column to be split
   * @param keepInCurrent Specifies one of the two new columns (the other 
   *  is the complement relative to the current column).
   */
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
  
  /**
   * 
   * @param e
   * @return
   */
  public boolean containsEdge(Edge e)
  {
    return getColumn(e,true) == getColumn(e,false);
  }
  
  /**
   * Restrict (project) this alignment to the given subset.
   * Discards all alignment links where both end points are not
   * in the restriction.
   * 
   * Modifies this object in place.
   * 
   * @param msa
   * @param sequenceIdsRestriction
   * @return
   */
  public static MSAPoset restrict(MSAPoset msa, Set<SequenceId> sequenceIdsRestriction)
  {
    Map<SequenceId,Sequence> sequences = new HashMap<SequenceId, Sequence>();
    Set<SequenceId> inter = Sets.intersection(sequenceIdsRestriction, msa.sequences().keySet());
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
  
  /**
   * Creates a new alignment where the sequences and edges are given by the union
   * of the provided MSAs.
   * 
   * @param collection
   * @return
   */
  public static MSAPoset union(Collection<MSAPoset> collection)
  {
    Map<SequenceId,Sequence> sequences = new HashMap<SequenceId, Sequence>();
    for (MSAPoset msa : collection)
      for (SequenceId key : msa.sequences().keySet())
      {
        Sequence current = sequences.get(key);
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
   *  Attempts to add all the links in the input column.
   *  Returns if this was succesful
   *  
   *  Warning: if it fails, the set of columns will be the same as before calling,
   *  but the linearization might be different (but still guaranteed to be consistent).
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
  
  private class ImplicitPoset implements PartialOrder<Column>, Serializable 
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
    if (b) return columnMaps.get(alignmentLink.sequenceId1())[alignmentLink.index1()];
    else   return columnMaps.get(alignmentLink.sequenceId2())[alignmentLink.index2()];
  }
  
  /**
   * 
   * @param sequenceId
   * @param index
   * @return The column in which the given character index in the given sequenceId belongs to
   */
  public Column column(SequenceId sequenceId, int index)
  {
    return columnMaps.get(sequenceId)[index];
  }
  
  Set<Map<SequenceId,Integer>> points()
  {
    Set<Map<SequenceId,Integer>> result = Sets.newLinkedHashSet();
    for (Column c : linearizedColumns.keySet())
      result.add(c.points);
    return result;
  }

  
  private static List<Edge> spanningEdges(Map<SequenceId,Integer> points)
  {
    List<Edge> result = Lists.newArrayList();
    if (points.size() < 2)
      return result;
    List<SequenceId> sequenceIds = Lists.newArrayList(points.keySet());
    final SequenceId baseLang = sequenceIds.get(0);
    final int basePos = points.get(baseLang);
    for (int i = 1; i < sequenceIds.size(); i++)
    {
      final SequenceId otherLang = sequenceIds.get(i);
      result.add(new Edge(basePos, points.get(otherLang), baseLang, otherLang));
    }
    return result;
  }

  /**
   * A column is an equivalence class of (sequenceId, position)'s.
   * 
   * The interpretation is that all members of this equivalence class are 
   * hypothesized to have a common ancestor.
   * 
   * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
   *
   */
  public static class Column implements Serializable
  {
    private static final long serialVersionUID = 1L;
    Map<SequenceId,Integer> points;
    
    /**
     * a map from sequenceId to a POSITION in the string
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
    public Column(SequenceId sequenceId, int i)
    {
      points = new HashMap<SequenceId, Integer>(2, 0.75f);
      points.put(sequenceId,i);
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
  
  /**
   * 
   * @return A set of columns in no particular order.
   */
  public Collection<Column> columns()
  {
    return Collections.unmodifiableCollection(linearizedColumns.keySet());
  }
  
  /**
   * A list of columns in a MSA-linearized order (i.e. topological order with
   * respect to the order of all sequences
   * @return
   */
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
  
  
  
  public double basicStat(SequenceId t1, SequenceId t2, boolean isIdentity) // o.w., just check aligned
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
        stat.addValue( charAt(c, t1).equals(charAt(c, t2)) ? 1.0 : 0.0);
    }
    return stat.getMean();
  }

  @Override
  public String toString() { return toString(null); }
  public String toString(Set<SequenceId> restriction)
  {
    StringBuilder [] builders = createPaddedStrings(restriction);
    StringBuilder result = new StringBuilder();
    List<SequenceId> sequenceIdPrintOrder = printOrder(restriction);
    for (int i = 0; i < sequenceIdPrintOrder.size(); i++)
    {
      result.append(builders[i]);
      result.append("|");
      result.append(sequenceIdPrintOrder.get(i));
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
    Map<SequenceId,Integer> sequenceIdPrintOrder = SupportFunctions.invert(printOrder);
    StringBuilder [] builders = new StringBuilder[printOrder.size()];
    
    int padLength = -1;
    for (SequenceId sequenceId : sequenceIdPrintOrder.keySet())
    {
      padLength = Math.max(padLength, sequences.get(sequenceId).alphabet.getMaxLetterStringLength());
      final int row = sequenceIdPrintOrder.get(sequenceId);
      StringBuilder current = new StringBuilder();
      builders[row] = current;
    }
    final String gap = "-" + Strings.repeat(" ", Math.max(0,padLength - 1));
    for (Column c : linearizedColumns())
      if (restriction == null || BriefCollections.intersects(restriction, c.points.keySet()))
      {
        for (SequenceId sequenceId : sequenceIdPrintOrder.keySet())
        {
          String currentChar = c.points.keySet().contains(sequenceId)   ?
              "" + Alphabet.toPaddedString(sequences.get(sequenceId).letterAt(c.points.get(sequenceId)), padLength) :
              gap;
          builders[sequenceIdPrintOrder.get(sequenceId)].append(currentChar);
        }
      }
    return builders;
  }
  
  /**
   * Save the alignment in fasta format.
   * 
   * Uses '-' for gaps.
   * 
   * @param f
   */
  public void toFASTA(File f)
  {
    PrintWriter out = BriefIO.output(f); 
    StringBuilder [] builders = createPaddedStrings(null);
    List<SequenceId> sequenceIdPrintOrder = printOrder(null);
    for (int i = 0; i < sequenceIdPrintOrder.size(); i++)
    {
      out.println(">" + sequenceIdPrintOrder.get(i));
      out.println(builders[i]);
    }
    out.close();
  }
  
  /**
   * Parse an alignment in fasta format.
   * 
   * @param f
   * @return
   */
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
          throw new RuntimeException("Duplicated sequenceId:" + currentTaxon);
        alignData.put(currentTaxon, new ArrayList<Integer>());
        stringData.put(currentTaxon, new StringBuilder());
      }
      else if (line.matches("[a-zA-Z.-]*"))
      {
        if (currentTaxon == null)
          throw new RuntimeException("Sequences should be preceded by the sequenceId using a line " +
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
    Map<SequenceId,Sequence> strings = Maps.newLinkedHashMap();
    Alphabet alphabet = new Alphabet();
    for (SequenceId sequenceId :  stringData.keySet())
      strings.put(sequenceId, Sequence.buildSimpleSequence(alphabet, stringData.get(sequenceId).toString()));
    MSAPoset result = new MSAPoset(strings);
    result.disableLinearization();
    int len = -1;
    for (SequenceId sequenceId : alignData.keySet())
    {
      if (len == -1)
        len = alignData.get(sequenceId).size();
      else if (len != alignData.get(sequenceId).size())
        throw new RuntimeException("Invalid alignment spec: " +
        		"all gap-padded seqns should have the len");
    }
    
    List<SequenceId> sequenceIds = Lists.newArrayList(alignData.keySet());
    for (int p = 0; p < len; p++)
    {
      Map<SequenceId,Integer> points = Maps.newLinkedHashMap();
      for (int l = 0; l < sequenceIds.size(); l++)
      {
        SequenceId l1 = sequenceIds.get(l);
        if (alignData.get(l1).get(p) != null)
          points.put(l1, alignData.get(l1).get(p));
      }
      result.tryAdding(points);
    }
    result.enableLinearization();
    return result;
  }

  public static boolean deepEquals(MSAPoset msa1, MSAPoset msa2)
  {
    if (!msa1.sequences().equals(msa2.sequences())) return false;
    return Sets.newLinkedHashSet(msa1.edges()).equals(Sets.newLinkedHashSet(msa2.edges()));
  }

}
