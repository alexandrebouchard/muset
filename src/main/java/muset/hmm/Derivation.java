package muset.hmm;


import java.io.Serializable;
import java.util.Arrays;




/**
 * Derivation tells you how phonemes in the current word were derived from the 
 * ancestral word
 * @author Alexandre Bouchard
 */
public final class Derivation implements Serializable
{
  private static final long serialVersionUID = -8207543758229656446L;
  private final int [] ancestors;
  private final String ancestorWord, currentWord;
  public static final int INSERTED = -1;

  public Derivation(final int[] ancestors, final String ancestorWord, 
      final String currentWord)
  {
    this.ancestors = ancestors;
    this.ancestorWord = ancestorWord;
    this.currentWord = currentWord;
    if (!isConsitent()) 
      throw new RuntimeException("Inconsistent derivation");
  }
  @Override
  public int hashCode()
  {
    final int PRIME = 31;
    int result = 1;
    result = PRIME * result + ancestorWord.hashCode();
    result = PRIME * result + Arrays.hashCode(ancestors);
    result = PRIME * result + currentWord.hashCode();
    return result;
  }
  @Override
  public boolean equals(Object o)
  {
    if (!(o instanceof Derivation)) return false;
    Derivation od = (Derivation) o;
    return od.ancestorWord.equals(ancestorWord) && od.currentWord.equals(currentWord)
      && Arrays.equals(od.ancestors, ancestors);
  }
  public String alignmentToString() 
  { 
    return Arrays.toString(ancestors).replaceAll("" + INSERTED, "X"); 
  }
  public boolean isConsitent()
  {
    return ancestorWord != null && currentWord != null && 
      ancestors.length == currentWord.length();
  }
  public String getAncestorWord() { return ancestorWord; }
  public String getCurrentWord() { return currentWord; }

  /**
   * returns the index of the ancestor character c_a of the provided character c_b,
   * provided c_a mutated to c_b
   * 
   * Throws a runtime exception if the index provided was actually inserted
   * 
   * @param phonemeIndex
   * @return
   */
  public int ancestor(int phonemeIndex)
  {
    int result = ancestors[phonemeIndex];
    if (result == INSERTED) throw new RuntimeException("Called ancestor(int) with inserted index");
    return result;
  }
  public boolean hasAncestor(int phonemeIndex)
  {
    if (phonemeIndex >= ancestors.length)
      throw new RuntimeException("Out of bound phoneme index in hasAncestor()");
    return  (ancestors[phonemeIndex] != INSERTED);
  }
  /**
   * returns the derivation from the ancestral word to the current word  
   * (pretending the time is reversed...)
   * @return
   */
  public Derivation invert()
  {
    int [] reversed = new int[ancestorWord.length()];
    for (int i = 0; i < reversed.length; i++) reversed[i] = INSERTED;
    for (int i = 0; i < ancestors.length; i++)
      if (hasAncestor(i))
        reversed[ancestor(i)] = i;
    // order of the last two arguments reversed on purpose!
    Derivation result = new Derivation(reversed, currentWord, ancestorWord);
    return result; 
  }
  public boolean isMonotonic()
  {
    int last = 0;
    for (int i = 0; i < currentWord.length(); i++)
      if (hasAncestor(i) && ancestor(i) < last) return false;
      else if (hasAncestor(i)) last = ancestor(i);
    return true;
  }
}