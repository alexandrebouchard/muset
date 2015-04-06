package muset;

import java.util.ArrayList;
import java.util.List;

import muset.Alphabet.Letter;


/**
 * A list of Letter's
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class Sequence
{
  private final List<Letter> letters;
  public final Alphabet alphabet;
  
  public Sequence(Alphabet alphabet, List<Letter> letters)
  {
    super();
    this.letters = letters;
    this.alphabet = alphabet;
  }

  public Letter letterAt(int i)
  {
    return letters.get(i);
  }
  
  /**
   * 
   * @param letter
   * @return A new Sequence containing all the letters in this plus an additional one at the end.
   */
  public Sequence append(Letter letter)
  {
    List<Letter> result = new ArrayList<Letter>(letters);
    result.add(letter);
    return new Sequence(alphabet, result);
  }
  
  public int length()
  {
    return letters.size();
  }
  
  @Override
  public String toString()
  {
    return letters.toString();
  }

  @Override
  public int hashCode()
  {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((letters == null) ? 0 : letters.hashCode());
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
    Sequence other = (Sequence) obj;
    if (letters == null)
    {
      if (other.letters != null)
        return false;
    } else if (!letters.equals(other.letters))
      return false;
    return true;
  }

  /**
   * 
   * @param string
   * @return A sequence where each letter is on Character long
   */
  public static Sequence buildSimpleSequence(Alphabet alphabet, String string)
  {
    List<Letter> result = new ArrayList<Letter>(string.length());
    for (char c : string.toCharArray())
      result.add(alphabet.getLetter("" + c));
    return new Sequence(alphabet, result);
  }

  public Sequence subsequence(int fromIndex, int toIndex)
  {
    return new Sequence(alphabet, letters.subList(fromIndex, toIndex));
  }
}
