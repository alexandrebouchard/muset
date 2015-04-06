package muset;


import briefj.Indexer;


/**
 * Maintain a collection of letters. 
 * 
 * A letter is a generalization of character, which can
 * have a string representation longer then a single character.
 * 
 * Example: single nucleotides ('A', 'C', 'G', 'T'), codons 
 * ('AAA', 'AAC', ..), phonemes, etc.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class Alphabet
{
  /**
   * Bijection between {0, 1, 2, ..} to the letters created by this
   * instance
   */
  public final Indexer<Letter> indexer = new Indexer<Letter>();
  
  private int letterStringLength = 0;
  
  /**
   * 
   * @return The maximum length of the string representations of the
   *    letters created by this instance
   */
  public int getMaxLetterStringLength() 
  { 
    return letterStringLength; 
  }
  
  public Alphabet() {}
  
  /**
   * 
   * @param alphabet A copy of the provided alphabet
   */
  public Alphabet(Alphabet alphabet)
  {
    for (Letter letter : alphabet.indexer.objectsList())
      this.getLetter(letter.contents);
  }
  
  public boolean containsLetter(String letterContents)
  {
    Letter temp = new Letter(letterContents);
    return indexer.containsObject(temp);
  }

  /**
   * Return the letter after making sure that the letter is in the alphabet,
   * adding it if needed.
   * 
   * @param contents The string representation of the letter
   * @return
   */
  public Letter getLetter(String contents)
  {
    Letter result = new Letter(contents);
    if (!indexer.containsObject(result))
    {
      indexer.addToIndex(result);
      letterStringLength = Math.max(letterStringLength, contents.length());
    }
    return result;
  }
  
  /**
   * Return the letter or throw a RuntimeException if the letter is not
   * already in the alphabet.
   * @param contents
   * @return
   */
  public Letter getExistingLetter(String contents)
  {
    Letter result = new Letter(contents);
    if (!indexer.containsObject(result))
      throw new RuntimeException("Letter was assumed to exist but did not:" + contents);
    return result;
  }
  
  /**
   * Return a fixed length representation of the letter, used for alignments.
   * @param letter
   * @return
   */
  public String toPaddedString(Letter letter)
  {
    return toPaddedString(letter, letterStringLength);
  }
  
  public static String toPaddedString(Letter letter, int letterStringLength)
  {
    String result = letter.contents;
    while (result.length() < letterStringLength)
      result = result + " ";
    return result;
  }
  
  /**
   * A letter is a generalization of character, which can
   * have a string representation longer then a single character.
   * 
   * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
   *
   */
  public static class Letter
  {
    private final String contents;

    private Letter(String contents)
    {
      this.contents = contents;
    }
    
    @Override
    public String toString()
    {
      return contents;
    }

    @Override
    public int hashCode()
    {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((contents == null) ? 0 : contents.hashCode());
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
      Letter other = (Letter) obj;
      if (contents == null)
      {
        if (other.contents != null)
          return false;
      } else if (!contents.equals(other.contents))
        return false;
      return true;
    }
  }
}
