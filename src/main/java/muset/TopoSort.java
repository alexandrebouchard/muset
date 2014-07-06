package muset;
import java.io.*;
import java.util.*;

import org.apache.commons.lang3.tuple.Pair;




public class TopoSort
{
  public static interface PartialOrder<T>
  {
    public Set<T> next(T n);
    public Set<T> prev(T n);
    public Set<T> nodes();
  }
  
  /* modify current topo sort to accomodate an extra edge x -> y */
  /**
   * L is a total order space, e.g. double or integer
   */
  public static <T,L extends Comparable<L>> boolean onlineTopologicalSort(final PartialOrder<T> po, final Map<T,L> n2i, final T x, final T y)
  {
    if (x.equals(y)) return false;
    
    final L lb = n2i.get(y), ub = n2i.get(x);
    
    if (ub.compareTo(lb) < 0) //(ub < lb) 
      return true;
    
    final Set<T> 
      fwd = new HashSet<T>(4, 0.75f),
      bwd = new HashSet<T>(4, 0.75f);
    
    Deque<T> deque = new ArrayDeque<T>();
    deque.add(y);
    while (!deque.isEmpty())
    {
      final T node = deque.pop();
      fwd.add(node);
      for (final T next : po.next(node))
      {
        final L nextPos = n2i.get(next);
        if (nextPos.compareTo(ub) == 0) return false;
        if (!fwd.contains(next) && nextPos.compareTo(ub) < 0) //nextPos < ub)
          deque.push(next);
      }
    }
   
    deque = new ArrayDeque<T>();
    deque.add(x);
    while (!deque.isEmpty())
    {
      final T node = deque.pop();
      bwd.add(node);
      for (final T prev : po.prev(node))
      {
        final L prevPos = n2i.get(prev);
        if (!bwd.contains(prev) && prevPos.compareTo(lb) > 0) //prevPos > lb)
          deque.push(prev);
      }
    }
    
    final List<Pair<T,L>> 
      spfwd = sortedPairedList(fwd, n2i),
      spbwd = sortedPairedList(bwd, n2i);
    // merged contains the indices of the positions in which thing will be inserted
    final List<L> merged = new ArrayList<L>(spfwd.size() + spbwd.size());
    for (final Pair<T,L> p : spfwd)
      merged.add(p.getRight());
    for (final Pair<T,L> p : spbwd)
      merged.add(p.getRight());
    Collections.sort(merged);
    int currentIdxIdx = 0;
    for (final Pair<T,L> cur : spbwd)
      n2i.put(cur.getLeft(), merged.get(currentIdxIdx++));
    for (final Pair<T,L> cur : spfwd)
      n2i.put(cur.getLeft(), merged.get(currentIdxIdx++));
    
    return true;
    
//    OnlineSorter<T> sorter = new OnlineSorter<T>(po,x,y,n2i);
//    return sorter.doIt();
  }
  
  private static <T, L extends Comparable<L>> List<Pair<T,L>> sortedPairedList(Set<T> set, Map<T,L> n2i)
  {
    List<Pair<T,L>> list = new ArrayList<Pair<T,L>>(set.size());
    for (T elt : set)
      list.add(Pair.of(elt, n2i.get(elt)));
    Collections.sort(list, new Comparator<Pair<T,L>>() {
      @Override public int compare(Pair<T, L> arg0, Pair<T, L> arg1)
      {
        return arg0.getRight().compareTo(arg1.getRight());
      }
    });
    return list;
  }
  
  private static final class OnlineSorter<T>
  {
    private final PartialOrder<T> po;
    private final T x, y;
    private final Map<T,Integer> n2i;
    private final int lb, ub;
    
    public OnlineSorter(PartialOrder<T> po, T x, T y, Map<T, Integer> n2i)
    {
      this.po = po;
      this.x = x;
      this.y = y;
      this.n2i = n2i;
      this.lb = n2i.get(y);
      this.ub = n2i.get(x);
    }

    private final Set<T> 
      fwd = new HashSet<T>(),
      bwd = new HashSet<T>();
    
    private boolean doIt()
    {
      if (ub < lb) 
        return true;
      //
      if (!fwd2(y)) return false;
      bwd2(x);
      reassign();
      return true;
    }
    
    private void reassign()
    {
      List<Pair<T,Integer>> 
        spfwd = sortedPairedList(fwd, n2i),
        spbwd = sortedPairedList(bwd, n2i);
      // merged contains the indices of the positions in which thing will be inserted
      List<Integer> merged = new ArrayList<Integer>();
      for (Pair<T,Integer> p : spfwd)
        merged.add(p.getRight());
      for (Pair<T,Integer> p : spbwd)
        merged.add(p.getRight());
      Collections.sort(merged);
      int currentIdxIdx = 0;
      for (Pair<T,Integer> cur : spbwd)
        n2i.put(cur.getLeft(), merged.get(currentIdxIdx++));
      for (Pair<T,Integer> cur : spfwd)
        n2i.put(cur.getLeft(), merged.get(currentIdxIdx++));
    }



//    private boolean fwd(T node)
//    {
//      fwd.add(node);
//      for (T next : po.next(node))
//      {
//        final int nextPos = n2i.get(next);
//        if (nextPos == ub) return false;
//        if (!fwd.contains(next) && nextPos < ub)
//          if (!fwd(next))
//            return false;
//      }
//      return true;
//    }
    
    private boolean fwd2(T node)
    {
      Deque<T> deque = new ArrayDeque<T>();
      deque.add(node);
      while (!deque.isEmpty())
      {
        node = deque.pop();
        fwd.add(node);
        for (T next : po.next(node))
        {
          final int nextPos = n2i.get(next);
          if (nextPos == ub) return false;
          if (!fwd.contains(next) && nextPos < ub)
            deque.push(next);
        }
      }
      return true;
    }
    
    private void bwd2(T node)
    {
      Deque<T> deque = new ArrayDeque<T>();
      deque.add(node);
      while (!deque.isEmpty())
      {
        node = deque.pop();
        bwd.add(node);
        for (T prev : po.prev(node))
        {
          final int prevPos = n2i.get(prev);
          if (!bwd.contains(prev) && prevPos > lb)
            deque.push(prev);
        }
      }
    }
    
//    private void bwd(T node)
//    {
//      bwd.add(node);
//      for (T prev : po.prev(node))
//      {
//        final int prevPos = n2i.get(prev);
//        if (!bwd.contains(prev) && prevPos > lb)
//          bwd(prev);
//      }
//    }
  }
  
  
  public static <T> List<T> topologicalSort(PartialOrder<T> po)
  {
    Set<T> remaining = new HashSet<T>(po.nodes());
    List<T> result = new ArrayList<T>(remaining.size());
    while (!remaining.isEmpty())
    {
      // find a min
      Set<T> visited = new HashSet<T>();
      T elt = remaining.iterator().next();
      Set<T> temp = null;
      while (!(temp = po.prev(elt)).isEmpty())
      {
        elt = temp.iterator().next();
        if (visited.contains(elt)) return null;
        visited.add(elt);
      }
      // dfs from that min, adding stuff to the list at the same time
      dfs(po, elt,remaining,result);

    }
    Collections.reverse(result);
    if (!isCorrectTopoSort(po, result))
      return null;
    return result;
  }
  
  private static <T> void dfs(PartialOrder<T> po, T current, Set<T> remaining, List<T> result)
  {
    remaining.remove(current);
    for (T child : po.next(current))
      if (remaining.contains(child))
        dfs(po, child, remaining, result);
    result.add(current);
  }
  
  private static <T> void dfs2(PartialOrder<T> po, T _current, Set<T> remaining, List<T> result)
  {
    Deque<T> stack = new ArrayDeque<T>();
    stack.push(_current);
    while (!stack.isEmpty())
    {
      T current = stack.pop();
      if (!remaining.contains(current))
      {
        remaining.remove(current);
        result.add(current);
      }
    }
    Collections.reverse(result);
  }
  
  public static final class HashOrder<T> implements PartialOrder<T>, Serializable
  {
    private static final long serialVersionUID = 1L;
    private final Map<T,Set<T>> 
      nexts = new HashMap<T,Set<T>>(),
      prevs = new HashMap<T,Set<T>>();
    public HashOrder() {}
    public HashOrder(Set<T> initialObjects)
    {
      for (T o : initialObjects)
      {
        nexts.put(o, new HashSet<T>());
        prevs.put(o, new HashSet<T>());
      }
    }
    @Override public Set<T> next(T n)
    {
      return /*Collections.unmodifiableSet(*/nexts.get(n);//);
    }

    @Override public Set<T> nodes()
    {
      return /*Collections.unmodifiableSet(*/nexts.keySet();//);
    }

    @Override public Set<T> prev(T n)
    {
      return /*Collections.unmodifiableSet(*/prevs.get(n);//);
    }
    /* if needed, add elt */
    public void add(T x)
    {
      if (nexts.containsKey(x)) return;
      nexts.put(x, new HashSet<T>());
      prevs.put(x, new HashSet<T>());
    }
    public void add(T x, T y) // add x \prec y
    {
      add(x); 
      add(y);
      nexts.get(x).add(y);
      prevs.get(y).add(x);
    }
    public void add(Pair<T,T> p) { add(p.getLeft(),p.getRight()); }
    public void remove(T x, T y)
    {
      if (!nexts.containsKey(x) || !nexts.containsKey(y))
        throw new RuntimeException();
      nexts.get(x).remove(y);
      prevs.get(y).remove(x);
    }
    public void remove(Pair<T,T> p) { remove(p.getLeft(),p.getRight()); }
  }

//  public static <T> Pair<T,T> randomEdge(Random rand, PartialOrder<T> superOrder, HashOrder<T> subOrder)
//  {
//    List<T> nodes = new ArrayList<T>(superOrder.nodes());
//    Collections.shuffle(nodes, rand);
//    for (T node : nodes)
//    {
//      Set<T> inSub = subOrder.next(node);
//      for (T nxt : superOrder.next(node))
//        if (!inSub.contains(nxt))
//          return Pair.makePair(node, nxt);
//    }
//    return null;
//  }
  
  public static class SetOrder implements PartialOrder<Set<Integer>>
  {
    public final int size;
    public final boolean cyclic;
    public SetOrder(int size) { this.size = size; cyclic = false;  }
    public SetOrder(int size, boolean introduceCycle) {this.size = size; this.cyclic = introduceCycle; }
    @Override
    public Set<Set<Integer>> next(Set<Integer> n)
    {
      Set<Set<Integer>> result = new HashSet<Set<Integer>>();
      for (int i =0 ; i < size; i++)
        if (!n.contains(i))
        {
          Set<Integer> cur = new HashSet<Integer>(n.size()+1);
          cur.addAll(n);
          cur.add(i);
          result.add(cur);
        }
      if (cyclic && result.size() ==0) result.add(new HashSet<Integer>());
      return result;
    }

    @Override
    public Set<Set<Integer>> nodes()
    {
      Set<Set<Integer>> result = new HashSet<Set<Integer>>();
      nodes(new HashSet<Integer>(), result);
      return result;
    }

    private void nodes(Set<Integer> cur, Set<Set<Integer>> result)
    {
      result.add(cur);
      for (Set<Integer> successor : next(cur))
        if (!result.contains(successor))
          nodes(successor, result);
    }

    @Override
    public Set<Set<Integer>> prev(Set<Integer> n)
    {
      Set<Set<Integer>> result = new HashSet<Set<Integer>>();
      for (int i =0 ; i < size; i++)
        if (n.contains(i))
        {
          Set<Integer> cur = new HashSet<Integer>(n.size());
          cur.addAll(n);
          cur.remove(i);
          result.add(cur);
        }
      return result;
    }
    
  }
  
  public static <T> boolean isCorrectTopoSort(PartialOrder<T> po, List<T> proposed)
  {
    Map<T,Integer> inverted= MSAPoset.invert(proposed);
    for (T elt : proposed)
    {
      int current = inverted.get(elt);
      for (T other : po.next(elt))
        if (inverted.get(other) <= current)
          return false;
    }
    return true;
  }
  
  public static void main(String[]args)
  {
    SetOrder so = new SetOrder(5);
    System.out.println(so.nodes());
    System.out.println(topologicalSort(so));
    System.out.println("Testing online version...");
    HashOrder<Set<Integer>> hashOrder = new HashOrder<Set<Integer>>(so.nodes());
    List<Set<Integer>> all = new ArrayList<Set<Integer>>(so.nodes());
    Random rand = new Random(1);
    Map<Set<Integer>, Integer> current = MSAPoset.invert(topologicalSort(hashOrder));
    for (int i = 0; i < 1000; i++)
    {
      if (i == 42)
        System.out.println("going to crash");
      System.out.println(i + " Current po: " + hashOrder.nexts);
      Set<Integer> 
        x = all.get(rand.nextInt(all.size())),
        y = all.get(rand.nextInt(all.size()));
      System.out.println("Current edge considered: " + x + " -> "+y);
      boolean online = onlineTopologicalSort(hashOrder, current, x, y);
      System.out.println("Online: " + online);
      List<Set<Integer>> converted = convert(current);
      System.out.println("Online order: " + converted);
      
      hashOrder.add(x,y);
      
      // if so, check the ordering is right  
      if (online && !isCorrectTopoSort(hashOrder, converted))
        throw new RuntimeException();
      
      // check if it's true...
      if ((topologicalSort(hashOrder) == null ? false : true) != online)
        throw new RuntimeException();
      
      // only keep if is in so
      if (!so.next(x).contains(y))
      {
        System.out.println("rejected");
        hashOrder.remove(x,y);
      }
      else
        System.out.println("kept");
      
      System.out.println("----");
    }
  }

  private static <T> List<T> convert(Map<T, Integer> current)
  {
    List<T> result = new ArrayList<T>(current.size());
    for (int i =0; i < current.size();i++)
      result.add(null);
    for (T o : current.keySet())
      result.set(current.get(o), o);
    return result;
  }
}
